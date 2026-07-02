library(Spectra)
library(GlycoMsHelper)
# library(devtools)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(pROC)
library(eulerr)
library(patchwork)

# devtools::install_github("FujitaLab-Glycobiology/GlycoMsHelper")


#===================
# Read the MS file
#===================
slc35a2_file_path = '../raw_data/241114_SLC35A2_SeqTypsin.mzML'
slc35a2_file_name = stringr::str_remove(basename(slc35a2_file_path), "\\.mzML$")

slc35a2_mass_spectrum_data = Spectra::Spectra(slc35a2_file_path, source = MsBackendMzR())


#======================
# Check the .mzML file
#======================
GlycoMsHelper::MsFileChecker(slc35a2_mass_spectrum_data)


#======================
# ms file QC
#======================
slc35a2_qc_results = GlycoMsHelper::SpectrumQcFilter(
  ms_data = slc35a2_mass_spectrum_data, 
  filter_method_ms1 = c(peaksCount = 'mean_sd', totIonCurrent = 'quantile_prob', rtime = 'start_end'), 
  threshold_ms1 = list(peaksCount = 2, totIonCurrent = c(0.05, 1), rtime = c(5*60, 50*60)), 
  filter_method_ms2 = c(peaksCount = 'quantile_prob', totIonCurrent = 'quantile_prob', rtime = 'start_end'), 
  threshold_ms2 = list(peaksCount = c(0.1, 1), totIonCurrent = c(0.05, 1), rtime = c(5*60, 50*60)), 
  plot_option = T
)

slc35a2_mass_spectrum_data_filtered = slc35a2_qc_results$filtered_ms_data


#====================
# ms file de-noising 
#====================
ms2_denoising_info = list(
  spline_segmentation_regression = list(
    spar_start = -1.5, spar_end = 1.5, spar_step = 0.02, 
    RMSE_weight = 0.3, CV_weight = 0.3, D2_weight = 0.3, D1_weight = 0.1, 
    use_cv = T, top_n_to_remove = 5, 
    segmentated_non_linear_transform_fun = function(z) z^2+z
  ), 
  spline_regression = list(
    spar_start = -1.5, spar_end = 1.5, spar_step = 0.02, 
    RMSE_weight = 0.3, CV_weight = 0.3, D2_weight = 0.3, D1_weight = 0.1, 
    use_cv = T, top_n_to_remove = 5
  ), 
  segmentation_regression = list(
    segmentated_non_linear_transform_fun = function(z) z^2+z
  ), 
  quantile_prob = 0.05, 
  fixed_value = 30 
)

slc35a2_denoised_results = GlycoMsHelper::MS2SpectrumDenoising(
  ms_data = slc35a2_mass_spectrum_data_filtered, 
  ms2_spectrum_transform_method = 'log2_transform', 
  ms2_denoising_method = 'spline_segmentation_regression', 
  ms2_denoising_detail = ms2_denoising_info
) 

slc35a2_mass_spectrum_data_filtered_denoised = slc35a2_denoised_results$denoised_ms_data

# export(slc35a2_mass_spectrum_data_filtered_denoised,
#        backend = MsBackendMzR(),
#        file = paste0(slc35a2_file_name, '_qc_denoised.mzML'), BPPARAM = SerialParam())



#===========================
# construct the glycan lib 
#===========================
molecular_formula_all = c(
  # monosaccharides
  Hex = 'C6H10O5',
  HexNAc = 'C8H13N1O5',
  dHex = 'C6H10O4',
  Neu5Ac = 'C11H17N1O8',
  HexA = 'C6H8O6',
  Neu5Gc = 'C11H17N1O9',
  Pentose = 'C5H8O4',
  KDN = 'C9H14O8', 
  EtNP = 'C2H6N1O3P1', 
  AHM = 'C6H8O4', 
  # label
  ProA = 'C13H23N3O1', 
  AB = 'C7H10N2', 
  PA = 'C5H8N2', 
  # adduct
  H = 'H1',
  Na = 'Na1',
  K = 'K1',
  Li = 'Li1',
  Mg = 'Mg1'
)

N_glycan_lib = GlycoMsHelper::ConstructGlycanLibrary(
  glycan_type = 'N_glycan', 
  min_charge_state = 1, max_charge_state = 3, 
  derivatization_type = 'ProA', 
  adduct_type = c('H', 'Na', 'K'), 
  min_total_monosaccharides_num = 3, 
  max_total_monosaccharides_num = 22, 
  min_Hex_num = 1,      max_Hex_num = 12, 
  min_HexNAc_num = 2,   max_HexNAc_num = 10, 
  min_dHex_num = 0,     max_dHex_num = 4, 
  min_Neu5Ac_num = 0,   max_Neu5Ac_num = 4, 
  min_HexA_num = 0,     max_HexA_num = 0, 
  min_Neu5Gc_num = 0,   max_Neu5Gc_num = 0, 
  min_Pentose_num = 0,  max_Pentose_num = 0, 
  min_KDN_num = 0,      max_KDN_num = 0
)

# N_glycan_library = N_glycan_lib$glycan_monosaccharides_library
N_glycan_library_iso_info = GlycoMsHelper::GetMonoisoAndIsotopologueRatio(glycan_lib = N_glycan_lib$glycan_monosaccharides_library,
                                                                          molecular_names = colnames(N_glycan_lib$monosaccharides_adduct_num),
                                                                          molecular_formula_list = molecular_formula_all,
                                                                          monosaccharides_names = colnames(N_glycan_lib$monosaccharides_combination),
                                                                          threshold_iso_probalility = 0.01)

N_glycan_library = N_glycan_library_iso_info$glycan_monosaccharides_library_isotopic_info

# N_glycan_library = dplyr::mutate(N_glycan_library_iso_info$glycan_monosaccharides_library_isotopic_info, 
#   glycan_monoisotopic_mz = if_else(
#     theoretical_monoisotopic_isotopologue_abundance_ratio <= 0.85, 
#     theoretical_isotopologue_mz, glycan_monoisotopic_mz
#   ) 
# )
# 
# openxlsx::write.xlsx(N_glycan_library,
#                      file = 'N_glycan_library.xlsx')
N_glycan_library = N_glycan_library %>% 
  dplyr::filter(!(
    (H == 0 & Na == 1 & K == 2) | 
      (H == 0 & Na == 2 & K == 1) |
      (H == 0 & Na == 3 & K == 0) |
      (H == 0 & Na == 0 & K == 3)
  ))





#=================================================================
# Find the spectrum likely to be glycan based on diagnostic ions 
#=================================================================
diagnostic_frags = c(
  HexNAc =              204.08667, 
  HexNAc_ProA =         441.2708, 
  dHex_HexNAc_ProA =    587.3287, 
  HexNAc_HexNAc_ProA =  644.3502, 
  Hex_ProA =            400.2442, 
  Hex_Hex_ProA =        562.2970, 
  dHex_Hex_ProA =       546.3021, 
  Hex =                 163.06007, 
  Bi_HexNAc =           407.16607, 
  Bisecting =           1009.4824, 
  Bisecting_dHex =      1155.5403, 
  Hex_HexNAc_ProA =     603.3236
)

slc35a2_diagnostic_results = GlycoMsHelper::FindSpectrumByDiagnosticFragments(
  ms_data = slc35a2_mass_spectrum_data_filtered_denoised, 
  ms_data_raw = slc35a2_mass_spectrum_data_filtered, 
  diagnostic_frags_list = diagnostic_frags, 
  diagnostic_frags_exp = 'HexNAc & (HexNAc_ProA | dHex_HexNAc_ProA) & !Hex_HexNAc_ProA', 
    # 'HexNAc_ProA | dHex_HexNAc_ProA', # 'HexNAc & (HexNAc_ProA | dHex_HexNAc_ProA) & !Hex_HexNAc_ProA', 
  ppm_val = 70
)


slc35a2_likely_glycan_spectrum_info = slc35a2_diagnostic_results$spectrum_info

# export(slc35a2_diagnostic_results$selected_ms_data,
#        backend = MsBackendMzR(),
#        file = paste0(slc35a2_file_name, '_likely_glycan.mzML'), BPPARAM = SerialParam())


#===========================================
# match the ms2 spectrum info to glycan lib 
#===========================================
slc35a2_spectrum_matching_result = GlycoMsHelper::FindPossibleGlycanComposition(
  spectrum_info = slc35a2_likely_glycan_spectrum_info, 
  glycan_lib = N_glycan_library, 
  max_precursor_mz_ppm = 20, 
  max_possible_candidates_num = 3
)

# openxlsx::write.xlsx(slc35a2_spectrum_matching_result,
#                      file = paste0(slc35a2_file_name, '_ms2_spectrum_matching_result.xlsx'))



#================================================
# find the composition by isotopics distribution
#================================================
slc35a2_glycan_spectrum_composition_info = GlycoMsHelper::ValidateGlycanCompositionByIsotopePattern(
  spectrum_matching_info = slc35a2_spectrum_matching_result, 
  molecular_names = colnames(N_glycan_lib$monosaccharides_adduct_num), 
  molecular_formula_list = molecular_formula_all, 
  ms_data = slc35a2_mass_spectrum_data_filtered, 
  ms1_window_left = 1, 
  ms1_window_right = 2, 
  bin_width = 0.3, 
  threshold_iso_probalility = 0.01
)

# openxlsx::write.xlsx(slc35a2_glycan_spectrum_composition_info,
#                      file = paste0(slc35a2_file_name, '_ms2_spectrum_composition_info.xlsx'))

















slc35a2_ground_truth_df = readr::read_csv('241114_SLC35A2_SeqTypsin_ground_truth_info.csv')

#====================================
# glycan composition quantification 
#====================================
spectra_id_filtered = slc35a2_mass_spectrum_data_filtered[["spectrumId"]] 

rtime_ground_truth = c()
peakscount_ground_truth = c()
tic_ground_truth = c()
precursor_charge_ground_truth = c()
precursor_mz_ground_truth = c()

for (i in slc35a2_ground_truth_df$ms2_spectrum_id) {
  
  spectrum_idx = which(spectra_id_filtered == i)
  
  rtime_ground_truth = c(rtime_ground_truth, slc35a2_mass_spectrum_data_filtered[["rtime"]][spectrum_idx])
  peakscount_ground_truth = c(peakscount_ground_truth, slc35a2_mass_spectrum_data_filtered[["peaksCount"]][spectrum_idx])
  tic_ground_truth = c(tic_ground_truth, slc35a2_mass_spectrum_data_filtered[["totIonCurrent"]][spectrum_idx])
  precursor_charge_ground_truth = c(precursor_charge_ground_truth, slc35a2_mass_spectrum_data_filtered[["precursorCharge"]][spectrum_idx])
  precursor_mz_ground_truth = c(precursor_mz_ground_truth, slc35a2_mass_spectrum_data_filtered[["precursorMz"]][spectrum_idx])
  
}

slc35a2_ground_truth_info = dplyr::mutate(slc35a2_ground_truth_df, 
                                          r_time = rtime_ground_truth, 
                                          peaks_count = peakscount_ground_truth, 
                                          tic = tic_ground_truth, 
                                          precursor_charge = precursor_charge_ground_truth, 
                                          precursor_mz = precursor_mz_ground_truth, 
                                          r_time_min = rtime_ground_truth/60)

h_num = c()
na_num = c()
k_num = c()

for (i in unique(slc35a2_ground_truth_df$glycan_string)) {
  temp_glycan_lib = dplyr::filter(N_glycan_library, glycan_string == i)

  temp_data = dplyr::filter(slc35a2_ground_truth_info, glycan_string == i)
  
  for (pre_mz in temp_data$precursor_mz){
    lib_idx = which.min(abs(pre_mz - temp_glycan_lib$glycan_monoisotopic_mz))
    
    h_num = c(h_num, temp_glycan_lib[lib_idx, ]$H)
    na_num = c(na_num, temp_glycan_lib[lib_idx, ]$Na)
    k_num = c(k_num, temp_glycan_lib[lib_idx, ]$K)
  }
}


slc35a2_ground_truth_info = dplyr::mutate(slc35a2_ground_truth_info, 
                                          H = h_num, Na = na_num, K = k_num) %>% 
  dplyr::mutate(adduct_type = paste0(ifelse(H  > 0, paste0("H", H), ""), 
                                     ifelse(Na > 0, paste0("Na", Na), ""), 
                                     ifelse(K  > 0, paste0("K", K), ""))) %>% 
  dplyr::select(ms2_spectrum_id, glycan_string, r_time_min, dplyr::everything())

# write.csv(slc35a2_ground_truth_info, file = '241114_SLC35A2_SeqTypsin_ground_truth_info.csv')

# function=4 process=0 scan=1054 and function=4 process=0 scan=1055 were manually changed
slc35a2_ground_truth_info = read_csv("241114_SLC35A2_SeqTypsin_ground_truth_info.csv")


tic_all = sum(slc35a2_ground_truth_info$tic)

slc35a2_ground_truth_info_groupby_adduct = slc35a2_ground_truth_info %>%
  dplyr::group_by(glycan_string, adduct_type) %>%
  dplyr::summarise(
    tic_sum = sum(tic, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  dplyr::mutate(tic_relative_abundance = tic_sum*100/tic_all)

slc35a2_ground_truth_info_groupby_precursorcharge = slc35a2_ground_truth_info %>%
  dplyr::group_by(glycan_string, precursor_charge) %>%
  dplyr::summarise(
    tic_sum = sum(tic, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  dplyr::mutate(tic_relative_abundance = tic_sum*100/tic_all)

glycan_order <- slc35a2_ground_truth_info_groupby_adduct %>%
  dplyr::group_by(glycan_string) %>%
  dplyr::summarise(total_tic = sum(tic_sum)) %>%
  dplyr::arrange(total_tic) %>%
  dplyr::pull(glycan_string)

# 
# pdf(file = '241114_SLC35A2_SeqTypsin_quantification_groupby_adduct.pdf',
#     width = 10, height = 7)
custom_adduct_order <- c(
  "H3", "H2K1", "H1K2",
  "H2", "H1K1", "K2", "H1Na1", "Na2" ,
  "H1"
)
ggplot(
  slc35a2_ground_truth_info_groupby_adduct %>%
    mutate(
      glycan_string = factor(glycan_string, levels = glycan_order),
      adduct_type = factor(adduct_type, levels = custom_adduct_order)
    ),
  aes(y = glycan_string, x = tic_relative_abundance, fill = adduct_type)
) +
  geom_col(position = position_stack(reverse = T)) +
  scale_fill_manual(values = c(
    "H1" = "#A3CDA0", "H2" = "#6BAF68", "H3" = "#3D8138",
    "H1K1" = "#9CC1D6", "H2K1" = "#5D98BB", "H1K2" = "#2F7AA9", "K2" = "#0E5F95",
    "H1Na1" = "#E89F9E", "Na2" =  "#D17673"
  )) +
  theme_classic() + theme(legend.position = "none")
# dev.off()
# dev.off()

# 
# pdf(file = '241114_SLC35A2_SeqTypsin_quantification_groupby_precursor_charge.pdf',
#     width = 10, height = 7)
ggplot(
  slc35a2_ground_truth_info_groupby_precursorcharge %>%
    mutate(glycan_string = factor(glycan_string, levels = glycan_order)),
  aes(y = glycan_string, x = tic_relative_abundance, fill = factor(precursor_charge))
) +
  geom_col() +
  scale_fill_manual(values = c(
    '1' = "#DB9195",
    '2' = "#9CC1D6",
    '3' = "#83B881"
  )) +
  theme_classic() + theme(legend.position = "none")
# dev.off()
# dev.off()









#================================
# evaluate performance (MS2 ppm)
#================================
res_list <- vector("list", length(seq(10, 150, 10)))
ppm_vec <- seq(10, 150, 10)
for (idx in seq_along(ppm_vec)) {
  i <- ppm_vec[idx]
  
  slc35a2_diagnostic_results_ms2 = GlycoMsHelper::FindSpectrumByDiagnosticFragments(
    ms_data = slc35a2_mass_spectrum_data_filtered_denoised, ms_data_raw = slc35a2_mass_spectrum_data_filtered, 
    diagnostic_frags_list = diagnostic_frags, diagnostic_frags_exp = 'HexNAc & (HexNAc_ProA | dHex_HexNAc_ProA) & !Hex_HexNAc_ProA', 
    ppm_val = i
  )
  slc35a2_likely_glycan_spectrum_info_ms2 = slc35a2_diagnostic_results_ms2$spectrum_info
  # match the ms2 spectrum info to glycan lib 
  slc35a2_spectrum_matching_result_ms2 = GlycoMsHelper::FindPossibleGlycanComposition(
    spectrum_info = slc35a2_likely_glycan_spectrum_info_ms2, 
    glycan_lib = N_glycan_library, max_precursor_mz_ppm = 20, max_possible_candidates_num = 3
  )
  # find the composition by isotopics distribution
  slc35a2_glycan_spectrum_composition_info_ms2 = GlycoMsHelper::ValidateGlycanCompositionByIsotopePattern(
    spectrum_matching_info = slc35a2_spectrum_matching_result_ms2, molecular_names = colnames(N_glycan_lib$monosaccharides_adduct_num), 
    molecular_formula_list = molecular_formula_all, ms_data = slc35a2_mass_spectrum_data_filtered, 
    ms1_window_left = 1, ms1_window_right = 2, bin_width = 0.3, threshold_iso_probalility = 0.01)

  res_list[[idx]] = slc35a2_glycan_spectrum_composition_info_ms2 %>%
    dplyr::select(glycan_string, ms2_spectrum_id) %>% 
    dplyr::mutate(ppm = i)
}

final_result <- dplyr::bind_rows(res_list)
# write.csv(final_result, file = 'roc_data.csv')


ground_truth <- slc35a2_ground_truth_df %>% 
  dplyr::rename(glycan_truth = glycan_string)

slc35a2_mass_spectrum_data_filtered_ms2 = Spectra::filterMsLevel(slc35a2_mass_spectrum_data_filtered, msLevel. = 2)
ms2_id = slc35a2_mass_spectrum_data_filtered_ms2[["spectrumId"]]


roc_data_all = data.frame(
  spectrum_id = rep(ms2_id, times = length(ppm_vec)),
  ppm_val = rep(ppm_vec, each = length(ms2_id))
) %>% 
  dplyr::left_join(final_result, by = c('spectrum_id' = 'ms2_spectrum_id', 'ppm_val' = 'ppm')) %>% 
  dplyr::left_join(ground_truth, by = c('spectrum_id' = 'ms2_spectrum_id')) %>% 
  dplyr::mutate(
    class = dplyr::case_when(
      !is.na(glycan_string) & !is.na(glycan_truth) & glycan_string == glycan_truth ~ "tp", 
      !is.na(glycan_string) & !is.na(glycan_truth) & glycan_string != glycan_truth ~ "fp", 
      !is.na(glycan_string) & is.na(glycan_truth) ~ "fp",
      is.na(glycan_string) & !is.na(glycan_truth) ~ "fn",
      TRUE ~ "tn"
    )
  )

# F1-ppm plot
f1_score_summary <- roc_data_all %>%
  dplyr::group_by(ppm_val) %>%
  dplyr::summarise(
    tp = sum(class == "tp", na.rm = TRUE),
    fp = sum(class == "fp", na.rm = TRUE),
    fn = sum(class == "fn", na.rm = TRUE),
    tn = sum(class == "tn", na.rm = TRUE),
    f1 = 2*tp / (2*tp + fp + fn), 
    .groups = "drop"
  ) %>% 
  dplyr::add_row(ppm_val = 0, tp = 0, fp = 0, fn = 0, tn = 0, f1 = 0) %>% 
  dplyr::arrange(ppm_val)

best_point <- f1_score_summary %>%
  dplyr::filter(f1 == max(f1)) %>% 
  dplyr::filter(ppm_val == min(ppm_val))


pdf(file = '241114_SLC35A2_SeqTypsin_f1_ppm.pdf', 
    width = 7, height = 7)

ggplot(f1_score_summary, aes(x = ppm_val, y = f1)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  annotate("point", x = best_point$ppm_val, y = best_point$f1, 
           color = "red", size = 3) + 
  annotate("segment", x = best_point$ppm_val, xend = best_point$ppm_val,
           y = 0, yend = best_point$f1,
           linetype = "dashed", color = "red") +
  annotate("text", x = best_point$ppm_val, y = best_point$f1,
           label = best_point$f1,
           vjust = -1, color = "red", size = 4) + 
  scale_x_continuous(
    breaks = sort(unique(c(f1_score_summary$ppm_val,
                           best_point$ppm_val)))
  ) +
  theme_classic()

dev.off()
dev.off()


# precision recall plot
pr_summary <- roc_data_all %>%
  dplyr::group_by(ppm_val) %>%
  dplyr::summarise(
    tp = sum(class == "tp", na.rm = TRUE),
    fp = sum(class == "fp", na.rm = TRUE),
    fn = sum(class == "fn", na.rm = TRUE),
    tn = sum(class == "tn", na.rm = TRUE),
    tpr = tp / (tp + fn),   # sensitivity / recall
    fpr = fp / (fp + tn),   # false positive rate
    precision = tp/(tp+fp), 
    f1 = 2*tp / (2*tp + fp + fn), 
    .groups = "drop"
  ) %>%
  dplyr::arrange(tpr)

P <- pr_summary$tp[1] + pr_summary$fn[1]
N <- pr_summary$tn[1] + pr_summary$fp[1]

pr_extra <- tibble::tibble(
  tpr = c(0, 1),
  precision = c(1, P/(P+N))
)

pr_plot <- bind_rows(
  pr_extra[1,],  # (0,1)
  pr_summary %>% select(tpr, precision),
  pr_extra[2,]   # (1, prevalence)
)

best_point <- pr_summary %>%
  dplyr::filter(f1 == max(f1)) %>% 
  dplyr::filter(ppm_val == min(ppm_val))


pdf(file = '241114_SLC35A2_SeqTypsin_pr_curve.pdf', 
    width = 7, height = 7)

ggplot(pr_plot, aes(tpr, precision)) +
  geom_step(size = 1.2, color = "#2C7BB6") + 
  geom_point(data = best_point,
             aes(tpr, precision),
             color = "red", size = 3) + 
  geom_text(data = best_point, aes(tpr, precision, 
                                   label = paste0("Recall = ", round(tpr, 3), 
                                                  "\nPrecision = ", round(precision, 3), 
                                                  "\nF1 = ", round(f1, 3))), 
            vjust = -1, color = "red", size = 4) + 
  geom_hline(yintercept = P/(P+N), linetype = "dashed") + 
  scale_y_continuous(breaks = sort(unique(c(pretty(pr_plot$precision), P/(P+N))))) + 
  labs(x = "Recall", y = "Precision") +
  theme_classic()

dev.off()
dev.off()


# euler plot (venn)
ground_truth_composition = unique(slc35a2_ground_truth_df$glycan_string)

venn_data = roc_data_all %>% 
  dplyr::filter(ppm_val == 70)
identified_composition = unique(na.omit(venn_data$glycan_string))

fit <- eulerr::euler(list(
  GroundTruth = ground_truth_composition,
  Identified = identified_composition
))


pdf(file = '241114_SLC35A2_SeqTypsin_venn.pdf', 
    width = 7, height = 7)
plot(
  fit,
  fills =c("#E89F9E", "#8FC6E8"), # c("#9FC1CC", "#E8B7AC"), 
  alpha = 0.7,
  quantities = TRUE,
  labels = TRUE
)

dev.off()
dev.off()

# detail info for venn (pie chart)
groundtruth_only_comp <- setdiff(ground_truth_composition, identified_composition)

groundtruth_only_info = venn_data %>% 
  dplyr::filter(glycan_truth %in% groundtruth_only_comp)

identified_as_glycan_wrong_comp = c()
not_identified = c()

for (i in unique(groundtruth_only_info$glycan_truth)) {
  temp_data = dplyr::filter(groundtruth_only_info, glycan_truth == i)
  
  if (any(!is.na(temp_data$glycan_string))) {
    identified_as_glycan_wrong_comp = c(identified_as_glycan_wrong_comp, i)
  } else {
    not_identified = c(not_identified, i)
  }
}

groundtruth_only_detail = data.frame(group = c('identified_as_glycan_wrong_comp_num', 'not_identified_num'), 
                                     number = c(length(identified_as_glycan_wrong_comp), length(not_identified)))


pdf(file = '241114_SLC35A2_SeqTypsin_venn_groundtruth_only_detail.pdf', 
    width = 7, height = 7)

ggplot(groundtruth_only_detail, aes(x = "", y = number, fill = group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") + 
  geom_text(
    aes(label = number),  
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(values = c(
    'identified_as_glycan_wrong_comp_num' = "#E6B0AF", 
    'not_identified_num' = "#CCA4B8"
  )) +
  theme_void()

dev.off()
dev.off()


temp_ground_truth = slc35a2_ground_truth_info %>% 
  dplyr::select(ms2_spectrum_id, r_time_min, precursor_charge, precursor_mz,adduct_type)

temp_glycomshelper_result = slc35a2_glycan_spectrum_composition_info %>% 
  dplyr::select(ms2_spectrum_id, glycan_string, theoretical_monoisotopic_mz, H, Na, K, total_charge) %>% 
  dplyr::mutate(adduct_type = paste0(ifelse(H  > 0, paste0("H", H), ""), 
                                   ifelse(Na > 0, paste0("Na", Na), ""), 
                                   ifelse(K  > 0, paste0("K", K), ""))) %>% 
  dplyr::select(ms2_spectrum_id, glycan_string, adduct_type, theoretical_monoisotopic_mz, total_charge) %>% 
  dplyr::rename_with(~ paste0(.x, "_glycomshelper"))


identified_only_comp <- setdiff(identified_composition, ground_truth_composition)

identified_only_info = venn_data%>% 
  dplyr::filter(glycan_string %in% identified_only_comp) %>% 
  dplyr::left_join(temp_ground_truth, by = c('spectrum_id' = 'ms2_spectrum_id')) %>% 
  dplyr::left_join(temp_glycomshelper_result, by = c('spectrum_id' = 'ms2_spectrum_id_glycomshelper'))


identified_as_glycan_wrong_comp = c()
identified_as_glycan_not_glycan = c()

for (i in unique(identified_only_info$glycan_string)) {
  temp_data = dplyr::filter(identified_only_info, glycan_string == i)
  
  if(any(!is.na(temp_data$glycan_truth))) {
    identified_as_glycan_wrong_comp = c(identified_as_glycan_wrong_comp, i)
  } else {
    identified_as_glycan_not_glycan = c(identified_as_glycan_not_glycan, i)
  }
}

identified_only_detail = data.frame(group = c('identified_as_glycan_wrong_comp_num', 'identified_as_glycan_not_glycan_num'), 
                                     number = c(length(identified_as_glycan_wrong_comp), length(identified_as_glycan_not_glycan)))


pdf(file = '241114_SLC35A2_SeqTypsin_venn_identified_only_detail.pdf', 
    width = 7, height = 7)

ggplot(identified_only_detail, aes(x = "", y = number, fill = group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") + 
  geom_text(
    aes(label = number),  
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(values = c(
    'identified_as_glycan_wrong_comp_num' = "#9CC1D6", 
    'identified_as_glycan_not_glycan_num' = "#70B7C1" 
  )) +
  theme_void()

dev.off()
dev.off()




























#================================
# evaluate performance (MS1 ppm)
#================================
res_list_ms1 <- vector("list", length(seq(10, 150, 10)))
ppm_vec_ms1 <- seq(10, 150, 10)
for (idx in seq_along(ppm_vec_ms1)) {
  i <- ppm_vec_ms1[idx]
  
  slc35a2_likely_glycan_spectrum_info_ms1 = slc35a2_diagnostic_results$spectrum_info
  # match the ms2 spectrum info to glycan lib 
  slc35a2_spectrum_matching_result_ms1 = GlycoMsHelper::FindPossibleGlycanComposition(
    spectrum_info = slc35a2_likely_glycan_spectrum_info_ms1, 
    glycan_lib = N_glycan_library, max_precursor_mz_ppm = i, max_possible_candidates_num = 3
  )
  # find the composition by isotopics distribution
  slc35a2_glycan_spectrum_composition_info_ms1 = GlycoMsHelper::ValidateGlycanCompositionByIsotopePattern(
    spectrum_matching_info = slc35a2_spectrum_matching_result_ms1, molecular_names = colnames(N_glycan_lib$monosaccharides_adduct_num), 
    molecular_formula_list = molecular_formula_all, ms_data = slc35a2_mass_spectrum_data_filtered, 
    ms1_window_left = 1, ms1_window_right = 2, bin_width = 0.3, threshold_iso_probalility = 0.01)
  
  res_list_ms1[[idx]] = slc35a2_glycan_spectrum_composition_info_ms1 %>%
    dplyr::select(glycan_string, ms2_spectrum_id) %>% 
    dplyr::mutate(ppm = i)
}

final_result_ms1 <- dplyr::bind_rows(res_list_ms1)
# write.csv(final_result_ms1, file = 'roc_data.csv')


ground_truth <- slc35a2_ground_truth_df %>% 
  dplyr::rename(glycan_truth = glycan_string)

slc35a2_mass_spectrum_data_filtered_ms2 = Spectra::filterMsLevel(slc35a2_mass_spectrum_data_filtered, msLevel. = 2)
ms2_id = slc35a2_mass_spectrum_data_filtered_ms2[["spectrumId"]]

roc_data_all_ms1 = data.frame(
  spectrum_id = rep(ms2_id, times = length(ppm_vec_ms1)),
  ppm_val = rep(ppm_vec_ms1, each = length(ms2_id))
) %>% 
  dplyr::left_join(final_result_ms1, by = c('spectrum_id' = 'ms2_spectrum_id', 'ppm_val' = 'ppm')) %>% 
  dplyr::left_join(ground_truth, by = c('spectrum_id' = 'ms2_spectrum_id')) %>% 
  dplyr::mutate(
    class = dplyr::case_when(
      !is.na(glycan_string) & !is.na(glycan_truth) & glycan_string == glycan_truth ~ "tp", 
      !is.na(glycan_string) & !is.na(glycan_truth) & glycan_string != glycan_truth ~ "fp", 
      !is.na(glycan_string) & is.na(glycan_truth) ~ "fp",
      is.na(glycan_string) & !is.na(glycan_truth) ~ "fn",
      TRUE ~ "tn"
    )
  )


# F1-ppm plot
f1_score_summary_ms1 <- roc_data_all_ms1 %>%
  dplyr::group_by(ppm_val) %>%
  dplyr::summarise(
    tp = sum(class == "tp", na.rm = TRUE),
    fp = sum(class == "fp", na.rm = TRUE),
    fn = sum(class == "fn", na.rm = TRUE),
    tn = sum(class == "tn", na.rm = TRUE), 
    tpr = tp / (tp + fn),   # sensitivity / recall
    fpr = fp / (fp + tn),   # false positive rate
    precision = tp/(tp+fp), 
    f1 = 2*tp / (2*tp + fp + fn), 
    
    .groups = "drop"
  ) %>% 
  dplyr::add_row(ppm_val = 0, tp = 0, fp = 0, fn = 0, tn = 0, tpr = 0, fpr = 0, precision = 0, f1 = 0) %>% 
  dplyr::arrange(ppm_val)

best_point_ms1 <- f1_score_summary_ms1 %>%
  dplyr::filter(f1 == max(f1)) %>% 
  dplyr::filter(ppm_val == min(ppm_val))


pdf(file = '241114_SLC35A2_SeqTypsin_f1_ppm_ms1.pdf',
    width = 7, height = 7)

ggplot(f1_score_summary_ms1, aes(x = ppm_val, y = f1)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  annotate("point", x = best_point_ms1$ppm_val, y = best_point_ms1$f1, 
           color = "red", size = 3) + 
  annotate("segment", x = best_point_ms1$ppm_val, xend = best_point_ms1$ppm_val,
           y = 0, yend = best_point_ms1$f1,
           linetype = "dashed", color = "red") +
  annotate("text", x = best_point_ms1$ppm_val, y = best_point_ms1$f1,
           label = best_point_ms1$f1,
           vjust = -1, color = "red", size = 4) + 
  scale_x_continuous(
    breaks = sort(unique(c(f1_score_summary_ms1$ppm_val,
                           best_point_ms1$ppm_val)))
  ) +
  theme_classic()

dev.off()
dev.off()

















#=====================================
# glycan lib plot for paper
#=====================================
N_glycan_lib_paper = GlycoMsHelper::ConstructGlycanLibrary(
  glycan_type = 'N_glycan', 
  min_charge_state = 1, max_charge_state = 3, 
  derivatization_type = 'ProA', 
  adduct_type = c('H'), 
  min_total_monosaccharides_num = 3, 
  max_total_monosaccharides_num = 22, 
  min_Hex_num = 1,      max_Hex_num = 12, 
  min_HexNAc_num = 2,   max_HexNAc_num = 10, 
  min_dHex_num = 0,     max_dHex_num = 2, 
  min_Neu5Ac_num = 0,   max_Neu5Ac_num = 4, 
  min_HexA_num = 0,     max_HexA_num = 0, 
  min_Neu5Gc_num = 0,   max_Neu5Gc_num = 0, 
  min_Pentose_num = 0,  max_Pentose_num = 0, 
  min_KDN_num = 0,      max_KDN_num = 0
)

# N_glycan_library = N_glycan_lib$glycan_monosaccharides_library
N_glycan_library_iso_info_paper = GlycoMsHelper::GetMonoisoAndIsotopologueRatio(glycan_lib = N_glycan_lib_paper$glycan_monosaccharides_library,
                                                                          molecular_names = colnames(N_glycan_lib_paper$monosaccharides_adduct_num),
                                                                          molecular_formula_list = molecular_formula_all,
                                                                          monosaccharides_names = colnames(N_glycan_lib_paper$monosaccharides_combination),
                                                                          threshold_iso_probalility = 0.01)

# pdf(file = 'N_glycan_library_paper.pdf', width = 8, height = 4)
# N_glycan_library_iso_info_paper$isotopic_distribution_plot
# dev.off()
# dev.off()











#======================
# de-noising detail 
#======================
ms2_denoising_detail_default = list(spline_segmentation_regression = list(spar_start = -1.5, spar_end = 1.5, spar_step = 0.02,
                                                                          RMSE_weight = 0.3, CV_weight = 0.3, D2_weight = 0.3, D1_weight = 0.1,
                                                                          use_cv = T, top_n_to_remove = 5,
                                                                          segmentated_non_linear_transform_fun = function(z) z^2+z
))

GetDenoiseInfo = function(denoising_detail, denoising_method, ms2_spectrum_transform_method, ms_id, ms_data) {
# GetSplineSegmentationNoise = function(denoising_detail, denoising_method, transform_fun, ms_id, ms_rt) {
  
  GetSplineSparEvaValue = function(spline_fit_result, x_val, y_val, y_predict_val) {
    # add more point to get more precise estimation
    d1 = stats::predict(spline_fit_result,
                        seq(from = x_val[1], to = x_val[length(x_val)], length.out = length(x_val)*2), deriv = 1)$y
    d2 = stats::predict(spline_fit_result,
                        seq(from = x_val[1], to = x_val[length(x_val)], length.out = length(x_val)*2), deriv = 2)$y
    
    # RMSE
    RMSE_val = sqrt(mean((y_val - y_predict_val)^2))
    
    # cv: in the smooth.spline(): leave-one-out cross-validation(LOOCV)
    if (!is.null(spline_fit_result$cv.crit)) {
      loocv = spline_fit_result$cv.crit
    } else {
      loocv = NA
    }
    
    # first derivative
    # first derivative smoothness: total variation
    # (smaller changes in the derivative indicate greater smoothness)
    # because of the data was sorted by the intensity, so theoretically,
    # d1 should always be positive value
    d1_total_var = sum(abs(diff(d1)))
    
    # second derivative
    # rate of change of tangent slope
    dx = diff(seq(from = x_val[1], to = x_val[length(x_val)], length.out = length(x_val)*2))
    dx = c(dx, mean(dx))
    d2_smoothness = sum((d2^2) * dx)
    
    return(list(RMSE = RMSE_val, CV = loocv, D1 = d1_total_var, D2 = d2_smoothness))
  }
  
  
  
  FindSplineSpar = function(x, y,
                            spar_start, spar_end, spar_step,
                            RMSE_weight, CV_weight, D1_weight, D2_weight,
                            use_cv = TRUE, plot = TRUE) {
    stopifnot(length(x) == length(y))
    
    spar_seq = seq(spar_start, spar_end, by = spar_step)
    RMSE_w = RMSE_weight
    CV_w = CV_weight
    D1_w = D1_weight
    D2_w = D2_weight
    
    spar_data_info = data.frame(spar = numeric(), RMSE = numeric(), CV = numeric(), D1 = numeric(), D2 = numeric())
    
    for (i in spar_seq) {
      fit = stats::smooth.spline(x, y, spar = i, cv = use_cv)
      # if (inherits(fit, "try-error")) return(NULL)
      y_predict = stats::predict(fit, x)$y
      
      eva_value = GetSplineSparEvaValue(spline_fit_result = fit,
                                        x_val = x, y_val = y, y_predict_val = y_predict)
      
      new_spar_info = data.frame(spar = i, RMSE = eva_value$RMSE, CV = eva_value$CV,
                                 D1 = eva_value$D1, D2 = eva_value$D2)
      
      spar_data_info = dplyr::bind_rows(spar_data_info, new_spar_info)
    }
    
    # max min normalization
    max_min_norm = function(v) {
      (v - min(v)) / (max(v) - min(v) + 1e-12)
    }
    spar_data_info_norm = data.frame(spar = spar_data_info$spar,
                                     RMSE_norm = max_min_norm(spar_data_info$RMSE),
                                     CV_norm = max_min_norm(spar_data_info$CV),
                                     D1_norm = max_min_norm(spar_data_info$D1),
                                     D2_norm = max_min_norm(spar_data_info$D2))
    
    # evaluation and determine the best spar
    if (all(is.na(spar_data_info_norm$CV_norm))) {
      RMSE_w = RMSE_w + CV_w
      CV_w = 0
      spar_data_info_norm$CV_norm = 0
    }
    spar_data_info_norm = dplyr::mutate(spar_data_info_norm,
                                        evaluation_score = RMSE_norm*RMSE_w + CV_norm*CV_w + D1_norm*D1_w + D2_norm*D2_w)
    
    best_idx = which.min(spar_data_info_norm$evaluation_score)
    best_spar = spar_data_info_norm$spar[best_idx]
    best_fit = stats::smooth.spline(x, y, spar = best_spar, cv = use_cv)
    
    # plotting part
    final_plot = NULL
    
    if (plot) {
      base_plot <- function(data, mapping, title, y_lab) {
        ggplot2::ggplot(data, mapping) +
          ggplot2::geom_point(color = "black") +
          ggplot2::geom_vline(xintercept = best_spar, linetype = "dashed", color = "red") +
          ggplot2::theme_classic() +
          ggplot2::labs(title = title, x = "spar", y = y_lab)
      }
      
      # RMSE
      RMSE_plot = base_plot(data = spar_data_info,
                            mapping = ggplot2::aes(x=spar, y=RMSE),
                            title = "RMSE",
                            y_lab = "RMSE")
      # CV
      if (!all(is.na(spar_data_info$CV))) {
        
        cv_plot = base_plot(data = spar_data_info,
                            mapping = ggplot2::aes(x=spar, y=CV),
                            title = "Leave-one-out cross-validation (LOOCV)",
                            y_lab = "LOOCV")
      } else {
        
        cv_plot = ggplot2::ggplot() +
          ggplot2::annotate("text", x=0.5, y=0.5, label="CV not computed") + ggplot2::theme_classic()
      }
      # D2
      D2_plot = base_plot(data = spar_data_info,
                          mapping = ggplot2::aes(x=spar, y=D2),
                          #title = "D2 smoothness (∫(f'' )^2)",
                          title = "D2 smoothness",
                          y_lab = "D2_smoothness")
      # D1
      D1_plot = base_plot(data = spar_data_info,
                          mapping = ggplot2::aes(x=spar, y=D1),
                          title = "D1 total variation",
                          y_lab = "D1_total_var")
      
      final_plot = (RMSE_plot + cv_plot) / (D2_plot + D1_plot) +
        patchwork::plot_annotation(title = paste("Optimization for best spar:", round(best_spar, 4)),
                                   subtitle = "Red dashed line indicates the selected values")
    }
    
    return(list(
      best_spar = best_spar,
      spar_data_info = spar_data_info,
      spar_data_info_norm = spar_data_info_norm,
      spar_data_info_plot = final_plot
    )
    )
  }
  
  
  # checking the ms2_spectrum_transform_method
  if (is.character(ms2_spectrum_transform_method)) {
    ms2_spectrum_transform_fun = switch(ms2_spectrum_transform_method,
                                        log2_transform = function(z) log2(z+1),
                                        asinh_transform = function(z) asinh(z),
                                        non_transform = function(z) z,
                                        stop("Unknown method"))
  } else if (is.function(ms2_spectrum_transform_method)) {
    ms2_spectrum_transform_fun = ms2_spectrum_transform_method
  } else {
    stop("'ms2_spectrum_transform_method' must be a function or character", call. = FALSE)
  }
  
  ms_id_all = ms_data[["spectrumId"]]
  ms_idx = which(ms_id_all == ms_id)
  
  spectra_data_raw = as.data.frame(Spectra::peaksData(ms_data)[[ms_idx]])
  
  spectra_distinct = dplyr::distinct(spectra_data_raw, intensity, .keep_all = TRUE)
  
  # background noise level detection
  # the intensity is transformed!!!!!
  y = ms2_spectrum_transform_fun(sort(spectra_distinct$intensity))
  n_total = length(y)
  x_idx = 1:n_total
  
  # according to the manual annotation, it is found that segmented regression can only well describe the
  # spectrum that have good fragmentation pattern (ions are completely fragmented),
  # however, for the spectrum that do NOT have good fragmentation pattern, which means the ions are not completely fragmented,
  # signal ions have low intensity, for these 'NOT have good fragmentation pattern' spectrm,
  # using a cubic smoothing spline to do the signal regression
  # use this method to determine the fragmentation is good or not good
  
  top4_peaks = dplyr::arrange(spectra_distinct, dplyr::desc(intensity)) |>
    utils::head(4)
  
  # for some of the isotopic patterns, monoisotopic peak have similar intensity with M+1 peak, or even M+1 peak are higher than monoisotopic
  # if sort from the highest intensity to lower intensity, the difference(abs(diff(intensity_tops))) might be 0.5, 1, 0.5
  # so here sort() is used
  intensity_tops = sort(top4_peaks$mz)
  
  intensity_tops_diff = abs(diff(intensity_tops))
  
  if (length(intensity_tops_diff) >= 2) {
    # calculate cv to see if these m/z have regular pattern
    
    avg_diff <- mean(intensity_tops_diff)
    sd_diff  <- stats::sd(intensity_tops_diff)
    
    diff_cv <- sd_diff / avg_diff
    
    is_complex_pattern <- diff_cv < 0.05 && avg_diff > 0.09 && avg_diff < 1.11
    
  } else if (length(intensity_tops_diff) == 1) {
    is_complex_pattern <- all(abs(intensity_tops_diff - 1) < 0.1) ||
      all(abs(intensity_tops_diff - 0.5) < 0.05) ||
      all(abs(intensity_tops_diff - 0.33) < 0.03) ||
      all(abs(intensity_tops_diff - 0.25) < 0.02)
  } else {
    is_complex_pattern <- FALSE
  }
  
  # for some of the isotopic patterns, especially for the non-proton ions,
  # fragmentation is not complete, the precursor isotopes have the highest intensity,
  # and the intensity of these precursor isotopes sometimes are similar
  # so if do not remove the top 5 peaks, the minimum slope point might be among these top peaks
  if (is_complex_pattern) {
    y_ss <- y[1:max(1, (n_total - denoising_detail[[denoising_method]]$top_n_to_remove))]
    x_ss <- 1:length(y_ss)
    
    method_for_regression = 'Spline_regression'
    
    # determine the best spar parameter for the stats::smooth.spline()
    best_spar_info = FindSplineSpar(x_ss, y_ss,
                                    spar_start = denoising_detail[[denoising_method]]$spar_start,
                                    spar_end = denoising_detail[[denoising_method]]$spar_end,
                                    spar_step = denoising_detail[[denoising_method]]$spar_step,
                                    RMSE_weight = denoising_detail[[denoising_method]]$RMSE_weight,
                                    CV_weight = denoising_detail[[denoising_method]]$CV_weight,
                                    D2_weight = denoising_detail[[denoising_method]]$D2_weight,
                                    D1_weight = denoising_detail[[denoising_method]]$D1_weight,
                                    use_cv = denoising_detail[[denoising_method]]$use_cv,
                                    plot = T)
    
    try(
      {
        fit_ss = stats::smooth.spline(x_idx, y, spar = best_spar_info$best_spar, cv = denoising_detail[[denoising_method]]$use_cv)
        
        # regression spline and first-order derivative (slope)
        slope = stats::predict(fit_ss, x_ss, deriv = 1)$y
        
        threshold_idx = which.min(slope)
        
        # check the threshold index
        if (!is.null(threshold_idx) && threshold_idx > 0 && threshold_idx < n_total) {
          # intensity threshold
          threshold_intensity = sort(spectra_distinct$intensity)[threshold_idx]
        } else {
          stop(sprintf("No threshold for denoising was detected in %s spectrum.", ms_id),
               call. = FALSE)
        }
        
        fit_ss_eva_value = GetSplineSparEvaValue(
          spline_fit_result = fit_ss,
          x_val = x_idx,
          y_val = y,
          y_predict_val = stats::predict(fit_ss, x_idx)$y
        )
        
        new_regression = data.frame(ms2_spectrum_id = ms_id, 
                                    threshold_value = threshold_intensity,
                                    denoising_method = method_for_regression,
                                    spar = best_spar_info$best_spar,
                                    RMSE = fit_ss_eva_value$RMSE, LOOCV = fit_ss_eva_value$CV, D1 = fit_ss_eva_value$D1, D2 = fit_ss_eva_value$D2,
                                    
                                    r_squared_linear = 0, adj_r_squared_linear = 0, residual_standard_error_linear = 0,
                                    slope_linear = 0, intercept_linear = 0,
                                    r_squared_non_linear = 0, adj_r_squared_non_linear = 0, residual_standard_error_non_linear = 0,
                                    slope_non_linear = 0, intercept_non_linear = 0)
        
        
        x_ss_seq <- seq(min(x_idx), max(x_idx), length.out = 1000)
        y_ss_predict_plot <- predict(fit_ss, x_ss_seq)
        
        ss_plot_df <- data.frame(x = y_ss_predict_plot$x, y = y_ss_predict_plot$y) %>%
          dplyr::mutate(part = dplyr::if_else(x <= threshold_idx, "before", "after"))
        
        
        de_noise_plot = ggplot2::ggplot() + 
          ggplot2::geom_point(data = data.frame(x = x_idx, y = y), ggplot2::aes(x = x, y = y), color = "black", size = 1) + 
          ggplot2::geom_line(data = ss_plot_df[ss_plot_df$part == "before", ], 
                             ggplot2::aes(x = x, y = y), color = "#A8B1C1", linewidth = 2) + 
          ggplot2::geom_line(data = ss_plot_df[ss_plot_df$part == "after", ], 
                             ggplot2::aes(x = x, y = y), color = "#8DBFDF", linewidth = 2) + 
          ggplot2::geom_vline(xintercept = threshold_idx, linetype = "dashed", color = "#9FC97F", linewidth = 2) + 
          ggplot2::annotate("text", x = min(x_idx) + 0.35 * diff(range(x_idx)), y = max(y) * 0.78,
                            label = ms_id,
                            size = 8
          ) +   
          ggplot2::labs(x = "Index", y = "Transformed intensity") +   
          ggplot2::theme_classic(base_size = 18) +
          ggplot2::theme(legend.position = "none", axis.text = element_text(color = "black"), axis.title = element_text(color = "black"))
          
     

        
        
      }, silent = F)
    
  } else {
    method_for_regression = 'Segmentation_regresion'
    
    non_linear_func_expression = denoising_detail[[denoising_method]]$segmentated_non_linear_transform_fun
    
    if (is.character(non_linear_func_expression)) {
      non_linear_transform_fun = switch(non_linear_func_expression,
                                        square_transform = function(z) z^2,
                                        exponential_transform = function(z) 2^z,
                                        stop("Unknown method"), call. = FALSE)
    } else if (is.function(non_linear_func_expression)) {
      non_linear_transform_fun = non_linear_func_expression
    } else {
      stop("'segmentated_non_linear_transform_fun' must be a function or character", call. = FALSE)
    }
    
    try({
      # segmented fitting
      fit_lm_init <- stats::lm(y ~ x_idx)
      sum_fit_lm_init = summary(fit_lm_init)
      
      # find break point
      # possible modification: psi = list(x_idx = c(30, 60)), use 2 break point and choose the smaller one 
      fit_lm_seg = try(segmented::segmented(fit_lm_init, seg.Z = ~x_idx, psi = list(x_idx = length(y)/2)))
      if (inherits(fit_lm_seg, "try-error")) {
        
        stop("'segmented::segmented' try-error", call. = FALSE)
      }
      
      break_point = fit_lm_seg$psi[,"Est."]
      threshold_idx = round(break_point)
      
      # linear and non_linear fitting
      y_linear = y[1:threshold_idx]
      x_linear = x_idx[1:threshold_idx]
      
      fit_lm_linear = stats::lm(y_linear ~ x_linear)
      # fit_lm_linear = stats::lm(y[1:threshold_idx] ~ x_idx[1:threshold_idx])
      
      sum_fit_lm_linear = summary(fit_lm_linear)
      
      x_non_linear = x_idx[(threshold_idx+1):length(y)]
      x_idx_non_linear = non_linear_transform_fun(x_non_linear)
      y_non_linear = y[(threshold_idx+1):length(y)]
      
      fit_lm_non_linear = stats::lm(y_non_linear ~ x_idx_non_linear)
      # fit_lm_non_linear = stats::lm(log(y_non_linear) ~ x2)
      sum_fit_lm_non_linear = summary(fit_lm_non_linear)
      
      # summary(fit_lm_seg)
      
      # check the threshold index
      if (!is.null(threshold_idx) && threshold_idx > 0 && threshold_idx < n_total) {
        # intensity threshold
        threshold_intensity = sort(spectra_distinct$intensity)[threshold_idx]
      } else {
        stop(sprintf("No threshold for denoising was detected in %s spectrum.", ms_id),
             call. = FALSE)
      }
      
      new_regression = data.frame(ms2_spectrum_id = ms_id, 
                                  
                                  threshold_value = threshold_intensity,
                                  
                                  denoising_method = method_for_regression,
                                  spar = 0, RMSE = 0, LOOCV = 0, D1 = 0, D2 = 0,
                                  
                                  r_squared_linear = sum_fit_lm_linear$r.squared, adj_r_squared_linear = sum_fit_lm_linear$adj.r.squared, residual_standard_error_linear = sum_fit_lm_linear$sigma,
                                  slope_linear = fit_lm_linear$coefficients[[2]], intercept_linear = fit_lm_linear$coefficients[[1]],
                                  r_squared_non_linear = sum_fit_lm_non_linear$r.squared,
                                  adj_r_squared_non_linear = sum_fit_lm_non_linear$adj.r.squared,
                                  residual_standard_error_non_linear = sum_fit_lm_non_linear$sigma,
                                  slope_non_linear = fit_lm_non_linear$coefficients[[2]],
                                  #slope_non_linear_2 = fit_lm_non_linear$coefficients["x2_non_linear"],
                                  intercept_non_linear = fit_lm_non_linear$coefficients[[1]])
      
      
      seg_plot_df <- data.frame(x_idx = x_idx, y = y, part = ifelse(x_idx <= threshold_idx, "linear_part", "nonlinear_part"))
      
      # fit_lm_init df
      init_df <- data.frame(x_idx = x_idx)
      init_df$y_pred <- predict(fit_lm_init, newdata = init_df)
      
      # fit_lm_linear df (before breakpoint)
      linear_df <- data.frame(x_linear = x_idx[1:threshold_idx])
      linear_df$y_pred <- predict(fit_lm_linear, newdata = linear_df)
      
      # fit_lm_non_linear df (after breakpoint)
      nonlinear_df <- data.frame(x_idx = x_idx[(threshold_idx + 1):length(y)])
      nonlinear_transformed_df = data.frame(x_idx_non_linear = non_linear_transform_fun(nonlinear_df$x_idx))
      
      nonlinear_df$y_pred <- stats::predict(fit_lm_non_linear, newdata = nonlinear_transformed_df)
      
      
      de_noise_plot = ggplot2::ggplot(seg_plot_df, ggplot2::aes(x = x_idx, y = y)) + 
        ggplot2::geom_point(color = "black", size = 1, alpha = 1) + 
        # scale_color_manual(values = c(linear_part = "#A8B1C1", nonlinear_part = "black")) +  
        
        # fit_lm_init
        ggplot2::geom_line(data = init_df, ggplot2::aes(x = x_idx, y = y_pred), color = "#D7E7EC", linewidth = 2, alpha = 0.8) + 
        
        # fit_lm_linear (before breakpoint)
        ggplot2::geom_line(data = linear_df, ggplot2::aes(x = x_linear, y = y_pred), color = "#A8B1C1", linewidth = 2, alpha = 0.95) + 
        
        # fit_lm_non_linear (after breakpoint)
        ggplot2::geom_line(data = nonlinear_df, ggplot2::aes(x = x_idx, y = y_pred), color = "#D48C8C", linewidth = 2, alpha = 0.95) + 
        
        # breakpoint dashed line
        ggplot2::geom_vline(xintercept = break_point, color = "#9FC97F", linewidth = 2, linetype = "dashed" ) + 
        
        ggplot2::annotate("text", x = min(x_idx) + 0.35 * diff(range(x_idx)), y = max(y) * 0.78,
                 label = ms_id,
                 size = 8
        ) +   
        ggplot2::labs(x = "Index", y = "Transformed intensity") +   
        ggplot2::theme_classic(base_size = 18) +
        ggplot2::theme(legend.position = "none", axis.text = element_text(color = "black"), axis.title = element_text(color = "black"))
      
      
    }, silent = F)
    
  }
  
  
  print(is_complex_pattern)

  top_peaks = spectra_data_raw %>%
    dplyr::arrange(dplyr::desc(intensity)) %>%
    dplyr::slice(1:15)
  
  if (is_complex_pattern) {
    #break_point = stats::quantile(spectra_data_raw$intensity, 9/10)
    break_point = sort(spectra_data_raw$intensity, decreasing = T)[10]
    
    # original plot
    spectra_data_raw_plot_raw = ggplot2::ggplot(spectra_data_raw, ggplot2::aes(x = mz, y = intensity)) + 
      ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0), linewidth = 0.5) + 
      ggplot2::geom_text(data = top_peaks, ggplot2::aes(label = round(mz, 3)), 
                         size = 3, vjust = 0, color = 'black') + 
      ggplot2::scale_y_continuous(expand = c(0, 0)) + 
      ggplot2::scale_x_continuous(expand = c(0.01, 0), limits = c(150, 2000)) + 
      ggplot2::labs(x = "m/z", y = "Intensity") +
      ggplot2::theme_classic()
    
    p1_raw <- spectra_data_raw_plot_raw +
      ggplot2::coord_cartesian(ylim = c(0, break_point))
    
    p2_raw <- spectra_data_raw_plot_raw +
      ggplot2::coord_cartesian(ylim = c(break_point, 1.05*max(spectra_data_raw$intensity)))
    
    spectra_data_raw_plot = p2_raw / p1_raw
    
    
    # denoised plot
    spectra_data_denoised_plot_raw = ggplot2::ggplot(spectra_data_raw, ggplot2::aes(x = mz, y = intensity)) + 
      ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0, color = intensity > threshold_intensity), 
                            linewidth = 0.5) + 
      ggplot2::geom_text(data = top_peaks, ggplot2::aes(label = round(mz, 3)), 
                         size = 3, vjust = 0, color = 'black') + 
      ggplot2::scale_color_manual(values = c("TRUE" = "black", "FALSE" = "#A8B1C1")) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) + 
      ggplot2::scale_x_continuous(expand = c(0.01, 0), limits = c(150, 2000)) + 
      ggplot2::labs(x = "m/z", y = "Intensity") +
      ggplot2::theme_classic() + 
      ggplot2::theme(legend.position = "none")
    
    p1_denoised <- spectra_data_denoised_plot_raw +
      ggplot2::coord_cartesian(ylim = c(0, break_point))
    
    p2_denoised <- spectra_data_denoised_plot_raw +
      ggplot2::coord_cartesian(ylim = c(break_point, 1.05*max(spectra_data_raw$intensity)))
    
    spectra_data_denoised_plot = p2_denoised / p1_denoised
    
  } else {
    # original plot
    spectra_data_raw_plot = ggplot2::ggplot(spectra_data_raw, ggplot2::aes(x = mz, y = intensity)) + 
      ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0), linewidth = 0.5) + 
      ggplot2::geom_text(data = top_peaks, ggplot2::aes(label = round(mz, 3)), 
                         size = 3, vjust = 0, color = 'black') + 
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05*max(spectra_data_raw$intensity))) + 
      ggplot2::scale_x_continuous(expand = c(0.01, 0), limits = c(150, 2000)) + 
      ggplot2::labs(x = "m/z", y = "Intensity") +
      ggplot2::theme_classic()
    
    # denoised plot
    spectra_data_denoised_plot = ggplot2::ggplot(spectra_data_raw, ggplot2::aes(x = mz, y = intensity)) + 
      ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0, color = intensity > threshold_intensity), 
                            linewidth = 0.5) + 
      ggplot2::geom_text(data = top_peaks, ggplot2::aes(label = round(mz, 3)), 
                         size = 3, vjust = 0, color = 'black') + 
      ggplot2::scale_color_manual(values = c("TRUE" = "black", "FALSE" = "#A8B1C1")) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1.05*max(spectra_data_raw$intensity))) + 
      ggplot2::scale_x_continuous(expand = c(0.01, 0), limits = c(150, 2000)) + 
      ggplot2::labs(x = "m/z", y = "Intensity") +
      ggplot2::theme_classic() + 
      ggplot2::theme(legend.position = "none")
  }
  
  
  
  if (is_complex_pattern) {
    return(
      list(
        threshold_value = threshold_intensity,
        row_regression_info = new_regression, 
        de_noise_info_plot = de_noise_plot, 
        mass_spectra_plot = spectra_data_raw_plot, 
        mass_spectra_denoised_plot = spectra_data_denoised_plot, 
        spline_spar_info_plot = best_spar_info$spar_data_info_plot
      )
    )
  } else {
    return(
      list(
        threshold_value = threshold_intensity,
        row_regression_info = new_regression, 
        de_noise_info_plot = de_noise_plot, 
        mass_spectra_plot = spectra_data_raw_plot, 
        mass_spectra_denoised_plot = spectra_data_denoised_plot
        #spline_spar_info_plot = best_spar_info$spar_data_info_plot
      )
    )
  }

  
}








# Hex3HexNAc5dHex1 2H
Hex3HexNAc5dHex1_2H = GetDenoiseInfo(denoising_detail = ms2_denoising_detail_default, 
               denoising_method = 'spline_segmentation_regression', 
               ms2_spectrum_transform_method = 'log2_transform', 
               ms_id = "function=2 process=0 scan=844", 
               ms_data = slc35a2_mass_spectrum_data_filtered) 


pdf(file = 'Hex3HexNAc5dHex1_2H_mass_spectra_denoised_plot.pdf', width = 18, height = 6)
Hex3HexNAc5dHex1_2H$mass_spectra_denoised_plot
dev.off()
dev.off()


pdf(file = 'Hex3HexNAc5dHex1_2Hde_noise_info_plot.pdf', width = 6, height = 6)
Hex3HexNAc5dHex1_2H$de_noise_info_plot
dev.off()
dev.off()

# Hex3HexNAc5dHex1 H+K
Hex3HexNAc5dHex1_H_K = GetDenoiseInfo(denoising_detail = ms2_denoising_detail_default, 
                      denoising_method = 'spline_segmentation_regression', 
                      ms2_spectrum_transform_method = 'log2_transform', 
                      ms_id = "function=2 process=0 scan=824", 
                      ms_data = slc35a2_mass_spectrum_data_filtered) 

pdf(file = 'Hex3HexNAc5dHex1_H_K_mass_spectra_denoised_plot.pdf', width = 18, height = 8)
Hex3HexNAc5dHex1_H_K$mass_spectra_denoised_plot
dev.off()
dev.off()


pdf(file = 'Hex3HexNAc5dHex1_H_K_de_noise_info_plot.pdf', width = 6, height = 6)
Hex3HexNAc5dHex1_H_K$de_noise_info_plot
dev.off()
dev.off()

pdf(file = 'Hex3HexNAc5dHex1_H_K_spline_spar_info_plot.pdf', , width = 8, height = 8)
Hex3HexNAc5dHex1_H_K$spline_spar_info_plot
dev.off()
dev.off()










#=========================================
# ms2 similairty cosine similarity score
#=========================================
ms2_spectrum_similarity_info = GlycoMsHelper::GetMS2SpectrumSimilarityScore(ms_data = slc35a2_mass_spectrum_data_filtered, 
                                                                            spectrum_matching_result = slc35a2_glycan_spectrum_composition_info, 
                                                                            glycan_composition_str = 'Hex3HexNAc5dHex1', 
                                                                            adduct_type = c(H = 2, Na = 0, K = 0), 
                                                                            bin_width = 1, 
                                                                            ms2_range_start = 100, 
                                                                            ms2_range_end = 2200) 

pdf(file = 'Hex3HexNAc5dHex1_2H_ms2_cosine_similarity_heatmap.pdf')
ms2_spectrum_similarity_info$similarity_score_heatmap
dev.off()
dev.off()


















spectra_data_ms1 = Spectra::filterMsLevel(ms_data, 1)
spectra_data_ms2 = Spectra::filterMsLevel(ms_data, 2)


all_ms2_peaks_data = Spectra::peaksData(spectra_data_ms2)
all_ms2_spectra_id = spectra_data_ms2[['spectrumId']]
all_ms2_spectra_rt = spectra_data_ms2[['rtime']]




current_ms2_peaks_data = as.data.frame(all_ms2_peaks_data[[i]])
current_ms2_peaks_data_distinct = dplyr::distinct(current_ms2_peaks_data, intensity, .keep_all = TRUE)

current_ms2_spectra_id = all_ms2_spectra_id[i]
current_ms2_spectra_rt = all_ms2_spectra_rt[i]

if (ms2_denoising_method == 'spline_segmentation_regression') {
  
  denoising_info = GetSplineSegmentationNoise(denoising_detail = ms2_denoising_detail,
                                              denoising_method = ms2_denoising_method,
                                              transform_fun = ms2_spectrum_transform_fun,
                                              spectra_distinct = current_ms2_peaks_data_distinct,
                                              ms_id = current_ms2_spectra_id,
                                              ms_rt = current_ms2_spectra_rt
  )
  
  thres_val = denoising_info$threshold_value
  new_row_spectrum_info = denoising_info$row_regression_info
  
} 















































stack <- c(
  "#9FC5DC",  # 很浅
  "#5F9FC4",  # 中浅
  "#2E7FB2",  # 中
  "#06639E"   

  "#A6D3A3",  # 浅绿
  "#6EB66A",  # 中绿
  "#3D8539"   # 深绿（原色）



  "#E59A97",  # 浅红
  "#CE5B58"   # 原色





  "#F0A3A1",  # 浅珊瑚红
  "#D97976"   # 中等红（比 #CE5B58 更浅）















test = slc35a2_ground_truth_info_groupby_adduct %>% 
  mutate(glycan_string_new = factor(glycan_string, levels = glycan_order))





ggplot(glycan_data,
       aes(x = glycan_string, y = abundance, fill = adduct_type)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





glycan_abundance = list()

for (i in unique(slc35a2_ground_truth_info$glycan_string)) {
  temp_data = dplyr::filter(slc35a2_ground_truth_info, glycan_string == i)
  glycan_abundance[[i]] = sum(temp_data$tic)
}

glycan_abundance_df = data.frame(
  glycan_string = names(glycan_abundance),
  sum_tic = unlist(glycan_abundance)
)

rownames(glycan_abundance_df) <- NULL

tic_all = sum(glycan_abundance_df$sum_tic)

glycan_abundance_df = dplyr::mutate(glycan_abundance_df, relative_abundance = sum_tic/tic_all)





ggplot(
  glycan_abundance_df,
  aes(y = reorder(glycan_string, relative_abundance),
      x = relative_abundance)
) +
  geom_col() +
  theme_classic()


























"#9CC1D6"

'#9ABDD1'


cyan_palette =   c("#C6DCDC", "#9FC4C4", "#7FA8A8", "#5E8787")
blue_palette =   c("#D2DCEB", "#AFC1DC", "#8FA8C6", "#6D87A8")
purple_palette = c("#E1D3E0", "#C9AFC7", "#B58FB2", "#946F92")
rose_palette =   c("#F0D2D2", "#E5AFAF", "#D48C8C", "#B46B6B")
beige_palette =  c("#F3E6D7", "#E6D0B8", "#D9BFA0", "#C7A888")

sage_palette =   c("#E3EDDE", "#C9DDBF", "#B2CEA4", "#93B486")
teal_palette =   c("#D7E7EC", "#B6D2DB", "#9FC1CC", "#7EA8B5")






identified_only_info = venn_data %>% 
  dplyr::filter(glycan_string %in% identified_only_comp)




groundtruth_only_info_pie_chart = list(unique())



df <- data.frame(
  group = c("TP", "FP", "FN", "TN"),
  value = c(479, 68, 458, 6863)
)

ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void()





length(groundtruth_only_data)
length(identified_only_data)







































cyan_palette =   c("#C6DCDC", "#9FC4C4", "#7FA8A8", "#5E8787")
blue_palette =   c("#D2DCEB", "#AFC1DC", "#8FA8C6", "#6D87A8")
purple_palette = c("#E1D3E0", "#C9AFC7", "#B58FB2", "#946F92")
rose_palette =   c("#F0D2D2", "#E5AFAF", "#D48C8C", "#B46B6B")
beige_palette =  c("#F3E6D7", "#E6D0B8", "#D9BFA0", "#C7A888")

sage_palette =   c("#E3EDDE", "#C9DDBF", "#B2CEA4", "#93B486")
teal_palette =   c("#D7E7EC", "#B6D2DB", "#9FC1CC", "#7EA8B5")
mauve_palette =  c("#E8D9E5", "#D3B7CF", "#C29DBE", "#A87FA4")
salmon_palette = c("#F4D9D2", "#E8B7AC", "#DB9A8F", "#C57C71")
orange_palette = c("#F6E0C8", "#EDC08F", "#E1A35C", "#C5823F")



df <- data.frame(
  group = c("TP", "FP", "FN", "TN"),
  value = c(479, 68, 458, 6863)
)

ggplot(df, aes(x = "", y = value, fill = group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  theme_void()















ms2_spectrum_similarity_info = GetMS2SpectrumSimilarityScore(ms_data = slc35a2_mass_spectrum_data_filtered, 
                                                             spectrum_matching_result = slc35a2_glycan_spectrum_composition_info, 
                                                             glycan_composition_str = 'Hex3HexNAc8dHex1', 
                                                             adduct_type = c(H = 2, Na = 0, K = 0), 
                                                             bin_width = 0.3, 
                                                             ms2_range_start = 100, 
                                                             ms2_range_end = 2200) 












testx = c(1, 2, 3, 4, 5, 6, 7, 8)
testy = c(1, 4.5, 10, 19, 30, 40, 60, 70)

testx_non_linear = testx^2
fit_test = stats::lm(testy ~ testx_non_linear)

pred_y <- stats::predict(fit_test, testx)$y


stats::predict(fit_test, testx_non_linear, deriv = 1)$y

pred_y$y





testx = c(1, 2, 3, 4, 5, 6, 7, 8)
testy = c(1, 4.5, 10, 19, 30, 40, 60, 70)

testx_non_linear = testx^2
fit_test = stats::lm(testy ~ testx_non_linear)

pred_y <- stats::predict(
  fit_test,
  newdata = data.frame(testx_non_linear = testx^2)
)

pred_y



ms2_id_gt = slc35a2_ground_truth_info$ms2_spectrum_id

test = slc35a2_denoised_results$denoising_regression_info %>% 
  dplyr::filter(denoising_method == 'Segmentation_regresion') %>% 
  filter(ms2_spectrum_id %in% ms2_id_gt)



ggplot(test, aes(x = slope_non_linear, y = slope_non_linear)) +
  geom_jitter(width = 0.2, size = 2) +
  theme_classic()





ms2_spectrum_similarity_info = GlycoMsHelper::GetMS2SpectrumSimilarityScore(ms_data = slc35a2_mass_spectrum_data_filtered, 
                                                             spectrum_matching_result = slc35a2_glycan_spectrum_composition_info, 
                                                             glycan_composition_str = 'Hex3HexNAc4dHex1', 
                                                             adduct_type = c(H = 2, Na = 0, K = 0), 
                                                             bin_width = 0.3, 
                                                             ms2_range_start = 100, 
                                                             ms2_range_end = 2200) 

col_fun = circlize::colorRamp2(
  seq(0, 1, length.out = 10),
  c(
    '#F8F8F8',
    '#F3EFD9',
    '#ECE4C9',
    '#E4D6B6',
    '#DBC6A4',
    '#D1B390',
    '#C69D82',
    '#B88476',
    '#A96C6E',
    '#955667'
  )
)





youden_df <- f1_score_summary_ms1 %>%
  mutate(youden_index = tpr - fpr)

# 找最佳ppm
best_point <- youden_df %>%
  filter(youden_index == max(youden_index, na.rm = TRUE))

# 作图
ggplot(youden_df, aes(x = ppm_val, y = youden_index)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  # 标出最佳点
  geom_point(data = best_point, color = "red", size = 3) +
  geom_vline(data = best_point, aes(xintercept = ppm_val), 
             linetype = "dashed") +
  
  labs(
    x = "ppm tolerance",
    y = "Youden index (TPR - FPR)",
    title = "Youden Index vs ppm tolerance"
  ) +
  theme_classic()




test <- f1_score_summary_ms1 %>%
  mutate(
    mcc = (tp * tn - fp * fn) /
      sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
  )

plot_df <- test %>%
  filter(!is.na(mcc), !is.nan(mcc))

# 找最大 MCC 点
best_point <- plot_df %>%
  filter(mcc == max(mcc, na.rm = TRUE))

ggplot(plot_df, aes(x = ppm_val, y = mcc)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  # 标出最优点
  geom_point(data = best_point, color = "red", size = 3) +
  geom_vline(data = best_point,
             aes(xintercept = ppm_val),
             linetype = "dashed") +
  
  labs(
    x = "ppm tolerance",
    y = "MCC (Matthews Correlation Coefficient)",
    title = "MCC vs ppm tolerance"
  ) +
  theme_classic()






