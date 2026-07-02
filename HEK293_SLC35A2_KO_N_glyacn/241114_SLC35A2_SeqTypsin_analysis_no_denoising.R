library(Spectra)
library(GlycoMsHelper)
# library(devtools)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(pROC)
library(eulerr)
library(patchwork)
library(readr)

# devtools::install_github("FujitaLab-Glycobiology/GlycoMsHelper")

setwd('D:/Paper/Manual_scripts/Research_paper/GlycoMSHelper/r_data')



#=================================================================
# Find the spectrum likely to be glycan based on diagnostic ions 
#=================================================================
slc35a2_diagnostic_results_no_denoising = GlycoMsHelper::FindSpectrumByDiagnosticFragments(
  ms_data = slc35a2_mass_spectrum_data_filtered, 
  ms_data_raw = slc35a2_mass_spectrum_data_filtered, 
  diagnostic_frags_list = diagnostic_frags, 
  diagnostic_frags_exp = 'HexNAc & (HexNAc_ProA | dHex_HexNAc_ProA) & !Hex_HexNAc_ProA', 
  # 'HexNAc_ProA | dHex_HexNAc_ProA', # 'HexNAc & (HexNAc_ProA | dHex_HexNAc_ProA) & !Hex_HexNAc_ProA', 
  ppm_val = 80
)


slc35a2_likely_glycan_spectrum_info_no_denoising = slc35a2_diagnostic_results_no_denoising$spectrum_info

# export(slc35a2_diagnostic_results_no_denoising$selected_ms_data,
#        backend = MsBackendMzR(),
#        file = paste0(slc35a2_file_name, '_likely_glycan.mzML'), BPPARAM = SerialParam())


#===========================================
# match the ms2 spectrum info to glycan lib 
#===========================================
slc35a2_spectrum_matching_result_no_denoising = GlycoMsHelper::FindPossibleGlycanComposition(
  spectrum_info = slc35a2_likely_glycan_spectrum_info_no_denoising, 
  glycan_lib = N_glycan_library, 
  max_precursor_mz_ppm = 20, 
  max_possible_candidates_num = 3
)

# openxlsx::write.xlsx(slc35a2_spectrum_matching_result,
#                      file = paste0(slc35a2_file_name, '_ms2_spectrum_matching_result.xlsx'))



#================================================
# find the composition by isotopics distribution
#================================================
slc35a2_glycan_spectrum_composition_info_no_denoising = GlycoMsHelper::ValidateGlycanCompositionByIsotopePattern(
  spectrum_matching_info = slc35a2_spectrum_matching_result_no_denoising, 
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



slc35a2_glycan_spectrum_composition_info_no_denoising_sum = slc35a2_glycan_spectrum_composition_info_no_denoising %>% 
  dplyr::mutate(adduct_type = paste0(ifelse(H  > 0, paste0("H", H), ""), 
                                     ifelse(K  > 0, paste0("K", K), ""), 
                                     ifelse(Na  > 0, paste0("Na", Na), ""))) %>% 
  dplyr::select(Hex, HexNAc, dHex, Neu5Ac, HexA, Neu5Gc, adduct_type, total_charge, 
                glycan_string, ms2_spectrum_id, ms2_retention_time, 
                ion_formula, theoretical_monoisotopic_mz, 
                ms1_spectrum_id, ms2_precursor_mz, ms2_total_ion_current) %>% 
  group_by(glycan_string, adduct_type) %>%
  summarise(
    across(c(Hex, HexNAc, dHex, Neu5Ac, HexA, Neu5Gc, total_charge, 
             ion_formula, theoretical_monoisotopic_mz), first),
    ms2_spectrum_ids = paste(ms2_spectrum_id, collapse = ", "),
    ms2_precursor_mzs = paste(ms2_precursor_mz, collapse = ", "), 
    ms2_retention_times = paste(ms2_retention_time, collapse = ", "), 
    ms1_spectrum_ids = paste(unique(ms1_spectrum_id), collapse = ", "), 
    ms2_tic_sum = sum(ms2_total_ion_current), 
    n_spectra = n(),
    .groups = "drop"
  )

write.csv(slc35a2_glycan_spectrum_composition_info_no_denoising_sum, file = 'no_denoise_241114_SLC35A2_SeqTypsin_glycan_composition_summary_GlycoMsHelper.csv')





#================================
# evaluate performance (MS2 ppm)
#================================
res_list_no_denoising <- vector("list", length(seq(10, 150, 10)))
ppm_vec <- seq(10, 150, 10)
for (idx in seq_along(ppm_vec)) {
  i <- ppm_vec[idx]
  
  slc35a2_diagnostic_results_ms2 = GlycoMsHelper::FindSpectrumByDiagnosticFragments(
    ms_data = slc35a2_mass_spectrum_data_filtered, ms_data_raw = slc35a2_mass_spectrum_data_filtered, 
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
  
  res_list_no_denoising[[idx]] = slc35a2_glycan_spectrum_composition_info_ms2 %>%
    dplyr::select(glycan_string, ms2_spectrum_id) %>% 
    dplyr::mutate(ppm = i)
}

final_result_no_denoising <- dplyr::bind_rows(res_list_no_denoising)
# write.csv(final_result, file = 'roc_data.csv')














# ground_truth <- slc35a2_ground_truth_df %>% 
#   dplyr::rename(glycan_truth = glycan_string)
# 
# slc35a2_mass_spectrum_data_filtered_ms2 = Spectra::filterMsLevel(slc35a2_mass_spectrum_data_filtered, msLevel. = 2)
# ms2_id = slc35a2_mass_spectrum_data_filtered_ms2[["spectrumId"]]
# 

roc_data_all_no_denoising = data.frame(
  spectrum_id = rep(ms2_id, times = length(ppm_vec)),
  ppm_val = rep(ppm_vec, each = length(ms2_id))
) %>% 
  dplyr::left_join(final_result_no_denoising, by = c('spectrum_id' = 'ms2_spectrum_id', 'ppm_val' = 'ppm')) %>% 
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
f1_score_summary_no_denoising <- roc_data_all_no_denoising %>%
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

best_point_no_denoising <- f1_score_summary_no_denoising %>%
  dplyr::filter(f1 == max(f1)) %>% 
  dplyr::filter(ppm_val == min(ppm_val))


pdf(file = '241114_SLC35A2_SeqTypsin_f1_ppm_no_denoising.pdf',
    width = 7, height = 7)

ggplot(f1_score_summary_no_denoising, aes(x = ppm_val, y = f1)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  annotate("point", x = best_point_no_denoising$ppm_val, y = best_point_no_denoising$f1, 
           color = "red", size = 3) + 
  annotate("segment", x = best_point_no_denoising$ppm_val, xend = best_point_no_denoising$ppm_val,
           y = 0, yend = best_point_no_denoising$f1,
           linetype = "dashed", color = "red") +
  annotate("text", x = best_point_no_denoising$ppm_val, y = best_point_no_denoising$f1,
           label = best_point_no_denoising$f1,
           vjust = -1, color = "red", size = 4) + 
  scale_x_continuous(
    breaks = sort(unique(c(f1_score_summary_no_denoising$ppm_val,
                           best_point_no_denoising$ppm_val)))
  ) +
  theme_classic()

dev.off()
dev.off()


# precision recall plot
pr_summary_no_denoising <- roc_data_all_no_denoising %>%
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

# P_no_denoising <- pr_summary_no_denoising$tp[1] + pr_summary_no_denoising$fn[1]
# N_no_denoising <- pr_summary_no_denoising$tn[1] + pr_summary_no_denoising$fp[1]

pr_extra_no_denoising <- tibble::tibble(
  tpr = c(0, 1),
  precision = c(1, dim(ground_truth)[1]/length(ms2_id))
)

pr_plot_no_denoising <- bind_rows(
  pr_extra_no_denoising[1,],  # (0,1)
  pr_summary_no_denoising %>% select(tpr, precision),
  pr_extra_no_denoising[2,]   # (1, prevalence)
)

best_point_no_denoising <- pr_summary_no_denoising %>%
  dplyr::filter(f1 == max(f1)) %>% 
  dplyr::filter(ppm_val == min(ppm_val))


# pdf(file = '241114_SLC35A2_SeqTypsin_pr_curve.pdf', 
#     width = 7, height = 7)
# 
ggplot(pr_plot_no_denoising, aes(tpr, precision)) +
  geom_step(size = 1.2, color = "#2C7BB6") +
  geom_point(data = best_point_no_denoising,
             aes(tpr, precision),
             color = "red", size = 3) +
  geom_text(data = best_point_no_denoising, aes(tpr, precision,
                                   label = paste0("Recall = ", round(tpr, 3),
                                                  "\nPrecision = ", round(precision, 3),
                                                  "\nF1 = ", round(f1, 3))),
            vjust = -1, color = "red", size = 4) +
  geom_hline(yintercept = dim(ground_truth)[1]/length(ms2_id), linetype = "dashed") +
  scale_y_continuous(breaks = sort(unique(c(pretty(pr_plot_no_denoising$precision), dim(ground_truth)[1]/length(ms2_id))))) +
  labs(x = "Recall", y = "Precision") +
  theme_classic()
# 
# dev.off()
# dev.off()







# combined precision recall plot (denoising and no denoising)
pr_plot_no_denoising$group <- "No Denoising"
pr_plot$group <- "Denoising"

pr_combined <- rbind(pr_plot_no_denoising, pr_plot)

best_point_no_denoising$group <- "No Denoising"
best_point$group <- "Denoising"

best_points_combined <- rbind(best_point_no_denoising, best_point)



pdf(file = '241114_SLC35A2_SeqTypsin_pr_curve_combined.pdf',
    width = 9, height = 7)

ggplot(pr_combined, aes(tpr, precision, color = group)) +
  geom_step(size = 1.2) +
  geom_point(data = best_points_combined,
             aes(tpr, precision), 
             color = 'red', 
             size = 3) +
  geom_text(data = best_points_combined,
            aes(tpr, precision, color = group,
                label = paste0("Recall = ", round(tpr, 3),
                               "\nPrecision = ", round(precision, 3),
                               "\nF1 = ", round(f1, 3))),
            vjust = -1, size = 4, show.legend = FALSE) +
  geom_hline(yintercept = dim(ground_truth)[1]/length(ms2_id), linetype = "dashed") + 
  scale_color_manual(values = c("No Denoising" = "#4575B4", "Denoising" = "#CC79A7")) +
  scale_y_continuous(breaks = sort(unique(c(pretty(pr_combined$precision), dim(ground_truth)[1]/length(ms2_id))))) +
  labs(x = "Recall", y = "Precision", color = "Group") +
  theme_classic()

dev.off()
dev.off()




# FP-FN plot
pr_summary_no_denoising$group <- "no_denoising"
pr_summary$group <- "denoising"

combined <- rbind(pr_summary_no_denoising, pr_summary)

combined <- combined %>% 
  dplyr::arrange(group, ppm_val)

best_points_FP_FN_plot <- combined %>%
  mutate(total_error = fp + fn) %>%
  group_by(group) %>%
  arrange(total_error, ppm_val, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

x_axis_min <- min(combined$fp, na.rm = TRUE)
y_axis_min <- min(combined$fn, na.rm = TRUE)

x_gap <- diff(range(combined$fp, na.rm = TRUE)) * 0.02
y_gap <- diff(range(combined$fn, na.rm = TRUE)) * 0.04

pdf(file = '241114_SLC35A2_SeqTypsin_FP-FN_plot.pdf',
    width = 8.5, height = 7)

ggplot(combined, aes(x = fp, y = fn, color = group)) +
  geom_path(size = 1.2) +
  geom_point(size = 2) +
  geom_text(aes(label = ppm_val), vjust = -0.8, size = 3, show.legend = FALSE) + 
  geom_segment(
    data = best_points_FP_FN_plot,
    aes(x = fp, xend = fp, y = y_axis_min, yend = fn - y_gap, color = group),
    linetype = "dashed",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  geom_segment(
    data = best_points_FP_FN_plot,
    aes(x = x_axis_min, xend = fp - x_gap, y = fn, yend = fn, color = group),
    linetype = "dashed",
    linewidth = 0.6,
    show.legend = FALSE
  ) + 
  geom_text(
    data = best_points_FP_FN_plot,
    aes(
      x = fp,
      y = y_axis_min,
      label = fp,
      color = group
    ),
    vjust = 1.6,
    size = 3.5,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = best_points_FP_FN_plot,
    aes(
      x = x_axis_min,
      y = fn,
      label = fn,
      color = group
    ),
    hjust = 1.2,
    size = 3.5,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) + 
  geom_point(
    data = best_points_FP_FN_plot,
    aes(x = fp, y = fn), color = 'red',
    size = 3,
    stroke = 1.2,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("no_denoising" = "#0072B2", "denoising" = "#CC79A7")) +
  labs(x = "FP", y = "FN", color = "Group") +
  theme_classic()

dev.off()
dev.off()
























# euler plot (venn)
ground_truth_composition = unique(slc35a2_ground_truth_df$glycan_string)

venn_data_no_denoising = roc_data_all_no_denoising %>% 
  dplyr::filter(ppm_val == 80)
identified_composition_no_denoising = unique(na.omit(venn_data_no_denoising$glycan_string))

fit_no_denoising <- eulerr::euler(list(
  GroundTruth = ground_truth_composition,
  Identified = identified_composition_no_denoising
))


pdf(file = '241114_SLC35A2_SeqTypsin_venn_no_denoising.pdf', 
    width = 7, height = 7)
plot(
  fit_no_denoising,
  fills =c("#E89F9E", "#8FC6E8"), # c("#9FC1CC", "#E8B7AC"), 
  alpha = 0.7,
  quantities = TRUE,
  labels = TRUE
)

dev.off()
dev.off()

# detail info for venn (pie chart)
groundtruth_only_comp_no_denoising <- setdiff(ground_truth_composition, identified_composition_no_denoising)

groundtruth_only_info_no_denoising = venn_data_no_denoising %>% 
  dplyr::filter(glycan_truth %in% groundtruth_only_comp_no_denoising)

identified_as_glycan_wrong_comp_no_denoising = c()
not_identified_no_denoising = c()

for (i in unique(groundtruth_only_info_no_denoising$glycan_truth)) {
  temp_data = dplyr::filter(groundtruth_only_info_no_denoising, glycan_truth == i)
  
  if (any(!is.na(temp_data$glycan_string))) {
    identified_as_glycan_wrong_comp_no_denoising = c(identified_as_glycan_wrong_comp_no_denoising, i)
  } else {
    not_identified_no_denoising = c(not_identified_no_denoising, i)
  }
}

groundtruth_only_detail_no_denoising = data.frame(group = c('identified_as_glycan_wrong_comp_num', 'not_identified_num'), 
                                     number = c(length(identified_as_glycan_wrong_comp_no_denoising), length(not_identified_no_denoising)))


pdf(file = '241114_SLC35A2_SeqTypsin_venn_groundtruth_only_detail_no_denoisng.pdf', 
    width = 7, height = 7)

ggplot(groundtruth_only_detail_no_denoising, aes(x = "", y = number, fill = group)) +
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














slc35a2_ground_truth_info = readr::read_csv("241114_SLC35A2_SeqTypsin_ground_truth_info.csv")



temp_ground_truth = slc35a2_ground_truth_info %>% 
  dplyr::select(ms2_spectrum_id, r_time_min, precursor_charge, precursor_mz,adduct_type)







temp_glycomshelper_result_no_denoising = slc35a2_glycan_spectrum_composition_info_no_denoising %>% 
  dplyr::select(ms2_spectrum_id, glycan_string, theoretical_monoisotopic_mz, H, Na, K, total_charge) %>% 
  dplyr::mutate(adduct_type = paste0(ifelse(H  > 0, paste0("H", H), ""), 
                                     ifelse(Na > 0, paste0("Na", Na), ""), 
                                     ifelse(K  > 0, paste0("K", K), ""))) %>% 
  dplyr::select(ms2_spectrum_id, glycan_string, adduct_type, theoretical_monoisotopic_mz, total_charge) %>% 
  dplyr::rename_with(~ paste0(.x, "_glycomshelper"))


identified_only_comp_no_denoising <- setdiff(identified_composition_no_denoising, ground_truth_composition)

identified_only_info_no_denoising = venn_data_no_denoising%>% 
  dplyr::filter(glycan_string %in% identified_only_comp_no_denoising) %>% 
  dplyr::left_join(temp_ground_truth, by = c('spectrum_id' = 'ms2_spectrum_id')) %>% 
  dplyr::left_join(temp_glycomshelper_result_no_denoising, by = c('spectrum_id' = 'ms2_spectrum_id_glycomshelper'))


identified_as_glycan_wrong_comp_no_denoising = c()
identified_as_glycan_not_glycan_no_denoising = c()

for (i in unique(identified_only_info_no_denoising$glycan_string)) {
  temp_data = dplyr::filter(identified_only_info_no_denoising, glycan_string == i)
  
  if(any(!is.na(temp_data$glycan_truth))) {
    identified_as_glycan_wrong_comp_no_denoising = c(identified_as_glycan_wrong_comp_no_denoising, i)
  } else {
    identified_as_glycan_not_glycan_no_denoising = c(identified_as_glycan_not_glycan_no_denoising, i)
  }
}

identified_only_detail_no_denoising = data.frame(group = c('identified_as_glycan_wrong_comp_num', 'identified_as_glycan_not_glycan_num'), 
                                    number = c(length(identified_as_glycan_wrong_comp_no_denoising), length(identified_as_glycan_not_glycan_no_denoising)))


pdf(file = '241114_SLC35A2_SeqTypsin_venn_identified_only_detail_no_denoising.pdf', 
    width = 7, height = 7)

ggplot(identified_only_detail_no_denoising, aes(x = "", y = number, fill = group)) +
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
# Hex3HexNAc2 H1K1
Hex3HexNAc2_H1K1 = GetDenoiseInfo(denoising_detail = ms2_denoising_detail_default, 
                                  denoising_method = 'spline_segmentation_regression', 
                                  ms2_spectrum_transform_method = 'log2_transform', 
                                  ms_id = "function=2 process=0 scan=536", 
                                  ms_data = slc35a2_mass_spectrum_data_filtered)


pdf(file = 'Hex3HexNAc2_H1K1_mass_spectra_denoised_plot.pdf', width = 15, height = 6)
Hex3HexNAc2_H1K1$mass_spectra_denoised_plot
dev.off()
dev.off()


# pdf(file = 'Hex3HexNAc5dHex1_2Hde_noise_info_plot.pdf', width = 6, height = 6)
Hex3HexNAc2_H1K1$de_noise_info_plot
dev.off()
dev.off()


id_all = slc35a2_mass_spectrum_data[["spectrumId"]]
ms_idx = which(id_all == "function=2 process=0 scan=536")
spectra_data_raw = as.data.frame(Spectra::peaksData(slc35a2_mass_spectrum_data)[[ms_idx]])















