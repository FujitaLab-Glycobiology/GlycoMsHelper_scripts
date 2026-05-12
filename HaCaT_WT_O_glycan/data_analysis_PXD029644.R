library(Spectra)
library(GlycoMsHelper)
# library(devtools)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(pROC)
library(eulerr)
library(patchwork)


setwd('D:/Paper/Manual_scripts/Research_paper/GlycoMSHelper/PXD029644')

#========================================
# construct the 2AB-labeled O-glycan lib 
#========================================
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
  AB = 'C7H10N2O1', 
  PA = 'C5H8N2', 
  # adduct
  H = 'H1',
  Na = 'Na1',
  K = 'K1',
  Li = 'Li1',
  Mg = 'Mg1'
)

O_glycan_lib = GlycoMsHelper::ConstructGlycanLibrary(
  glycan_type = 'O_glycan', 
  min_charge_state = 1, 
  max_charge_state = 1, 
  derivatization_type = 'AB', 
  adduct_type = c('H'), #adduct_type = c('H', 'Na', 'K'), 
  
  monosaccharides_additional_customized_list = c(Pentose = 132.0423), 
  min_num_monosaccharides_additional_customized_list = c(Pentose = 0), 
  max_num_monosaccharides_additional_customized_list = c(Pentose = 2),
  
  min_total_monosaccharides_num = 1, 
  max_total_monosaccharides_num = 10, 
  min_Hex_num = 0,      max_Hex_num = 5, 
  min_HexNAc_num = 0,   max_HexNAc_num = 5, 
  min_dHex_num = 0,     max_dHex_num = 3, 
  min_Neu5Ac_num = 0,   max_Neu5Ac_num = 3, 
  min_HexA_num = 0,     max_HexA_num = 2, 
  min_Neu5Gc_num = 0,   max_Neu5Gc_num = 0, 
  min_Pentose_num = 0,  max_Pentose_num = 0, 
  min_KDN_num = 0,      max_KDN_num = 0
)

O_glycan_library_iso_info = GlycoMsHelper::GetMonoisoAndIsotopologueRatio(glycan_lib = O_glycan_lib$glycan_monosaccharides_library, 
                                                                          molecular_names = colnames(O_glycan_lib$monosaccharides_adduct_num),
                                                                          molecular_formula_list = molecular_formula_all, 
                                                                          monosaccharides_names = colnames(O_glycan_lib$monosaccharides_combination), 
                                                                          threshold_iso_probalility = 0.01) 

O_glycan_library = O_glycan_library_iso_info$glycan_monosaccharides_library_isotopic_info

#write.csv(O_glycan_library, file = 'O_glycan_library.csv')


#=========================
# define diagnostic frags
#=========================
diagnostic_frags = c(
  Hex_2AB =         301.1394, 
  HexNAc_2AB =            342.1660, 
  dHex_2AB = 285.1445, 
  
  HexNAc =              204.08667, 
  'Neu5Ac-H2O' = 274.0921, 
  Neu5Ac = 292.1027, 
  HexNAc_Hex = 366.1395, 
  HexNAc_Hex_dHex = 512.1974, 
  'HexNAc_Hex_Neu5Ac-H2O' = 657.2349, 
  HexNAc2_Hex2 = 731.2717, 
  
  Hex =                 163.06007, 
  Bi_HexNAc =           407.16607
)






#===================
# Read the .mzml file
#===================
mzml_files <- list.files("./raw_data/", pattern = "\\.mzML$", 
                         full.names = T, ignore.case = TRUE)

#==================
# initialization
#==================
all_results = list()

for (mass_file in mzml_files) {
  
  file_name = stringr::str_remove(basename(mass_file), "\\.mzML$")
  
  mass_spectrum_data = Spectra::Spectra(mass_file, source = MsBackendMzR())
  
  #======================
  # Check the .mzML file
  #======================
  GlycoMsHelper::MsFileChecker(mass_spectrum_data)
  
  #======================
  # ms file QC
  #======================
  qc_results = GlycoMsHelper::SpectrumQcFilter(
    ms_data = mass_spectrum_data, 
    filter_method_ms1 = c(peaksCount = 'quantile_prob', totIonCurrent = 'quantile_prob', rtime = 'start_end'), 
    threshold_ms1 = list(peaksCount = c(0.05, 1), totIonCurrent = c(0.05, 1), rtime = c(19*60, 36*60)), 
    filter_method_ms2 = c(peaksCount = 'quantile_prob', totIonCurrent = 'quantile_prob', rtime = 'start_end'), 
    threshold_ms2 = list(peaksCount = c(0.05, 1), totIonCurrent = c(0.05, 1), rtime = c(20*60, 35*60)), 
    plot_option = T
  )
  
  mass_spectrum_data_filtered = qc_results$filtered_ms_data
  
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
  
  denoised_results = GlycoMsHelper::MS2SpectrumDenoising(
    ms_data = mass_spectrum_data_filtered, 
    ms2_spectrum_transform_method = 'log2_transform', 
    ms2_denoising_method = 'spline_regression', 
    ms2_denoising_detail = ms2_denoising_info
  ) 
  
  mass_spectrum_data_filtered_denoised = denoised_results$denoised_ms_data
  
  
  #=================================================================
  # Find the spectrum likely to be glycan based on diagnostic ions 
  #=================================================================
  diagnostic_results = GlycoMsHelper::FindSpectrumByDiagnosticFragments(
    ms_data = mass_spectrum_data_filtered_denoised, 
    ms_data_raw = mass_spectrum_data_filtered, 
    diagnostic_frags_list = diagnostic_frags, 
    diagnostic_frags_exp = 'HexNAc_2AB | Hex_2AB | dHex_2AB', 
    ppm_val = 50
  )
  
  # export(diagnostic_results$selected_ms_data,
  #        backend = MsBackendMzR(),
  #        file = "your_ms_file_path_and_name.mzML", BPPARAM = SerialParam())
  
  
  likely_glycan_spectrum_info = diagnostic_results$spectrum_info

  
  #================================================
  # Match the likely glycan spectrum to glycan lib
  #================================================
  likely_glycan_spectrum_matching_result = GlycoMsHelper::FindPossibleGlycanComposition(
    spectrum_info = likely_glycan_spectrum_info, 
    glycan_lib = O_glycan_library, 
    max_precursor_mz_ppm = 50, 
    max_possible_candidates_num = 3
  )

  #=======================================================================
  # Find the candidate glycan composition based on isotopics distribution
  #=======================================================================
  final_glycan_spectrum_matching_result = GlycoMsHelper::ValidateGlycanCompositionByIsotopePattern(
    spectrum_matching_info = likely_glycan_spectrum_matching_result, 
    molecular_names = colnames(O_glycan_lib$monosaccharides_adduct_num), 
    molecular_formula_list = molecular_formula_all, 
    ms_data = mass_spectrum_data_filtered, 
    ms1_window_left = 1, 
    ms1_window_right = 2, 
    bin_width = 0.05, 
    threshold_iso_probalility = 0.01
  )
  
  #==============
  # save results 
  #==============
  all_results[[file_name]] = list(
    spectrum_info                    = likely_glycan_spectrum_info,
    spectrum_matching_result         = likely_glycan_spectrum_matching_result,
    final_spectrum_matching_result   = final_glycan_spectrum_matching_result
  )
  
}



#====================
# paper compositions
#====================
table_s2_HaCaT_WT_composition = c(
  'HexNAc1', 
  'Hex1dHex1', 
  'Hex1HexNAc1', 
  'Pentose2Hex1', 
  'Hex2dHex1', 
  'Hex1HexNAc1dHex1', 
  'Hex2HexNAc1', 
  'Hex1HexNAc1Neu5Ac1', 
  'Hex2HexNAc1dHex1', 
  'Hex2HexNAc2', 
  'Hex2HexNAc1Neu5Ac1', 
  'Hex2HexNAc2dHex1', 
  'Hex2HexNAc2Neu5Ac1', 
  'Hex3HexNAc3Neu5Ac1', 
  'Hex3HexNAc3Neu5Ac2'
  )

table_s5_HEK293_WT_composition = c(
  'HexNAc1', 
  'Hex1dHex1', 
  'Hex1HexNAc1', 
  'Pentose2Hex1', 
  'Hex1Neu5Ac1', 
  'Hex2HexNAc1', 
  'Hex1HexNAc1Neu5Ac1', 
  'Hex2HexNAc2', 
  'Hex1HexNAc1Neu5Ac2', 
  'Hex2HexNAc2Neu5Ac1',
  'Hex2HexNAc2Neu5Ac2'
)

table_s6_keratinocyte_N_TERT_1_WT_50_composition = c(
  'HexNAc1', 
  'Pentose1Hex1', 
  'Hex1dHex1', 
  'Hex1HexNAc1', 
  'HexNAc2', 
  'Pentose2Hex1', 
  'Hex1Neu5Ac1',
  'Hex1HexNAc1dHex1', 
  'Hex2HexNAc1', 
  'Hex1HexNAc2', 
  'Hex1HexNAc1Neu5Ac1', 
  'Hex2HexNAc2', 
  'Hex2HexNAc1Neu5Ac1', 
  'Hex2HexNAc2dHex1', 
  'Hex1HexNAc1Neu5Ac2', 
  'Hex2HexNAc2Neu5Ac1', 
  'Hex3HexNAc3', 
  'Hex2HexNAc2dHex1Neu5Ac1', 
  'Hex2HexNAc2Neu5Ac2', 
  'Hex3HexNAc3Neu5Ac1', 
  'Hex3HexNAc3Neu5Ac2'
  )	
  
  











get_comp = function(file_names, final_matching_result) {
  all_ideintified_comp = c()
  for (name in file_names) {
    result <- final_matching_result[[name]]$final_spectrum_matching_result
    ideintified_comp = result[['glycan_string']]
    all_ideintified_comp = c(all_ideintified_comp, ideintified_comp)
  }
  return(all_ideintified_comp)
}


HaCaT_WT_day1_file_names = c(
  '210818_Lumos__LC1200__nanoFlex_NH21036_01_02_Ogly', 
  '210818_Lumos__LC1200__nanoFlex_NH21036_01_03_Ogly', 
  '210818_Lumos__LC1200__nanoFlex_NH21036_01_04_Ogly', 
  '210818_Lumos__LC1200__nanoFlex_NH21036_01_06_Ogly', 
  '210818_Lumos__LC1200__nanoFlex_NH21036_01_07_Ogly', 
  '210818_Lumos__LC1200__nanoFlex_NH21036_01_08_Ogly'
  )

HaCaT_WT_day2_file_names = c(
  '210818_Lumos__LC1200__nanoFlex_NH21036_02_02_Ogly', 
  '210818_Lumos__LC1200__nanoFlex_NH21036_02_03_Ogly', 
  '210818_Lumos__LC1200__nanoFlex_NH21036_02_04_Ogly', 
  '210818_Lumos__LC1200__nanoFlex_NH21036_02_06_Ogly', 
  '210818_Lumos__LC1200__nanoFlex_NH21036_02_07_Ogly', 
  '210818_Lumos__LC1200__nanoFlex_NH21036_02_08_Ogly'
  )

HaCaT_WT_target_file_names = c('211021_Lumos__LC1200__nanoFlex_NH21039_17_Ogly_targ')

HaCaT_WT_all_file_names = c(HaCaT_WT_day1_file_names, HaCaT_WT_day2_file_names, HaCaT_WT_target_file_names)


HaCaT_WT_all_ideintified_comp = get_comp(
  file_names = HaCaT_WT_all_file_names, 
  final_matching_result = all_results
  )

fit <- eulerr::euler(list(
  paper_comp = table_s2_HaCaT_WT_composition,
  identified = unique(HaCaT_WT_all_ideintified_comp)
))



pdf(file = 'HaCaT_WT_venn.pdf',
    width = 7, height = 7)
plot(
  fit,
  fills =c("#C2A5CF", "#A6D7B8"), # c("#9FC1CC", "#E8B7AC"), 
  alpha = 0.7,
  quantities = TRUE,
  labels = TRUE
)
dev.off()
dev.off()

HaCaT_WT_paper_only = setdiff(table_s2_HaCaT_WT_composition, unique(HaCaT_WT_all_ideintified_comp))

HaCaT_WT_glycomshelper_only = setdiff(unique(HaCaT_WT_all_ideintified_comp), table_s2_HaCaT_WT_composition)




not_identified_info <- read.xlsx('HaCaT_WT_glycomshelper_only_info.xlsx')

glycomshelper_only_info = data.frame(group = c('not_glycan_num', 'new_glycan_num'), 
                                     number = c(
                                       dim(dplyr::filter(not_identified_info, HaCaT_WT_glycomshelper_only_reason == 'not_glycan'))[1], 
                                       dim(dplyr::filter(not_identified_info, HaCaT_WT_glycomshelper_only_reason == 'glycan'))[1]
                                       )
                                     )


pdf(file = 'HaCaT_WT_glycomshelper_only_venn.pdf',
    width = 7, height = 7)
ggplot(glycomshelper_only_info, aes(x = "", y = number, fill = group)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") + 
  geom_text(
    aes(label = number),  
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(values = c(
    'not_glycan_num' = "#79C17C", 
    'new_glycan_num' = "#2AA268"#, 
    #'not_sure_glycan_num' = '#46C664'
  )) +
  theme_void()

dev.off()
dev.off()




#temp_lib = dplyr::filter(O_glycan_library, glycan_string %in% HaCaT_WT_glycomshelper_only)


c("#E64B35", "#4DBBD5", "#00A087", "#3C5488")
c("#4E79A7", "#F28E2B", "#59A14F", "#E15759")
c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")
c("#7BAFD4", "#C7849B", "#8FBC8F", "#D4A96A")
c("#0072B2", "#E69F00", "#009E73", "#CC79A7")

c("#E89F9E", "#8FC6E8", "#B9A3D6", "#9BC8A6")
c("#E89F9E", "#8FC6E8", "#C2A5CF", "#A6D7B8")


info_list = list()
for (i in HaCaT_WT_all_file_names) {
  final_mat = all_results[[i]]$final_spectrum_matching_result
  detail_info = dplyr::filter(final_mat, glycan_string == 'Hex3dHex1')
  info_list[[i]] = detail_info
}

info_list
  
info_list$`210818_Lumos__LC1200__nanoFlex_NH21036_01_03_Ogly`$ms2_spectrum_id
  























HEK293_WT_file_names = c(
  '210927_Lumos__LC1200__nanoFlex_NH21038_01_Ogly',
  '210928_Lumos__LC1200__nanoFlex_NH21038_22_Ogly',
  '210928_Lumos__LC1200__nanoFlex_NH21038_25_Ogly',
  '210916_Lumos__LC1200__nanoFlex_NH21038_25_Ogly_targ'
  
)

HEK293_WT_ideintified_comp = c()

for (name in HEK293_WT_file_names) {
  
  result <- all_results[[name]]$final_spectrum_matching_result
  
  ideintified_comp = result[['glycan_string']]
  
  HEK293_WT_ideintified_comp = c(HEK293_WT_ideintified_comp, ideintified_comp)
  
}



fit <- eulerr::euler(list(
  paper_comp = table_s5_HEK293_WT_composition,
  identified = unique(HEK293_WT_ideintified_comp)
))


plot(
  fit,
  fills =c("#E89F9E", "#8FC6E8"), # c("#9FC1CC", "#E8B7AC"), 
  alpha = 0.7,
  quantities = TRUE,
  labels = TRUE
)


identify_only = setdiff(unique(HEK293_WT_ideintified_comp), table_s5_HEK293_WT_composition)

result <- all_results[['210916_Lumos__LC1200__nanoFlex_NH21038_25_Ogly_targ']]$final_spectrum_matching_result

test = dplyr::filter(result, glycan_string %in% identify_only) 


all_id = test$ms2_spectrum_id



test = dplyr::filter(all_results[['210916_Lumos__LC1200__nanoFlex_NH21038_25_Ogly_targ']]$spectrum_matching_result, 
                     ms2_spectrum_id %in% all_id) %>% 
  dplyr::select(glycan_string, 
                ms2_retention_time, 
                ms2_spectrum_id, 
                glycan_monoisotopic_mz, 
                total_charge) %>% 
  dplyr::mutate(ms2_retention_time_min = ms2_retention_time/60)




















#======================
# de-noising detail
#======================
ms2_denoising_detail_default = list(spline_segmentation_regression = list(spar_start = -1.5, spar_end = 1.5, spar_step = 0.02,
                                                                          RMSE_weight = 0.3, CV_weight = 0.3, D2_weight = 0.3, D1_weight = 0.1,
                                                                          use_cv = T, top_n_to_remove = 5,
                                                                          segmentated_non_linear_transform_fun = function(z) z^2+z
))

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
  
  if (denoising_method == 'spline_regression') {
    is_complex_pattern <- TRUE
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

  # if (is_complex_pattern) {
  #   #break_point = stats::quantile(spectra_data_raw$intensity, 9/10)
  #   break_point = sort(spectra_data_raw$intensity, decreasing = T)[10]
  # 
  #   # original plot
  #   spectra_data_raw_plot_raw = ggplot2::ggplot(spectra_data_raw, ggplot2::aes(x = mz, y = intensity)) +
  #     ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0), linewidth = 0.5) +
  #     ggplot2::geom_text(data = top_peaks, ggplot2::aes(label = round(mz, 3)),
  #                        size = 3, vjust = 0, color = 'black') +
  #     ggplot2::scale_y_continuous(expand = c(0, 0)) +
  #     ggplot2::scale_x_continuous(expand = c(0.01, 0), limits = c(150, 2000)) +
  #     ggplot2::labs(x = "m/z", y = "Intensity") +
  #     ggplot2::theme_classic()
  # 
  #   p1_raw <- spectra_data_raw_plot_raw +
  #     ggplot2::coord_cartesian(ylim = c(0, break_point))
  # 
  #   p2_raw <- spectra_data_raw_plot_raw +
  #     ggplot2::coord_cartesian(ylim = c(break_point, 1.05*max(spectra_data_raw$intensity)))
  # 
  #   spectra_data_raw_plot = p2_raw / p1_raw
  # 
  # 
  #   # denoised plot
  #   spectra_data_denoised_plot_raw = ggplot2::ggplot(spectra_data_raw, ggplot2::aes(x = mz, y = intensity)) +
  #     ggplot2::geom_segment(ggplot2::aes(xend = mz, yend = 0, color = intensity > threshold_intensity),
  #                           linewidth = 0.5) +
  #     ggplot2::geom_text(data = top_peaks, ggplot2::aes(label = round(mz, 3)),
  #                        size = 3, vjust = 0, color = 'black') +
  #     ggplot2::scale_color_manual(values = c("TRUE" = "black", "FALSE" = "#A8B1C1")) +
  #     ggplot2::scale_y_continuous(expand = c(0, 0)) +
  #     ggplot2::scale_x_continuous(expand = c(0.01, 0), limits = c(150, 2000)) +
  #     ggplot2::labs(x = "m/z", y = "Intensity") +
  #     ggplot2::theme_classic() +
  #     ggplot2::theme(legend.position = "none")
  # 
  #   p1_denoised <- spectra_data_denoised_plot_raw +
  #     ggplot2::coord_cartesian(ylim = c(0, break_point))
  # 
  #   p2_denoised <- spectra_data_denoised_plot_raw +
  #     ggplot2::coord_cartesian(ylim = c(break_point, 1.05*max(spectra_data_raw$intensity)))
  # 
  #   spectra_data_denoised_plot = p2_denoised / p1_denoised
  # 
  # } else {
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
  # }



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
#





# denoise_info = GetDenoiseInfo(
#   denoising_detail = ms2_denoising_info, 
#   denoising_method = 'spline_regression', 
#   ms2_spectrum_transform_method = 'log2_transform', 
#   ms_id = 'controllerType=0 controllerNumber=1 scan=10298', 
#   ms_data = Spectra::Spectra("./raw_data/210818_Lumos__LC1200__nanoFlex_NH21036_01_02_Ogly.mzML", source = MsBackendMzR()))
# 
# denoise_info$mass_spectra_denoised_plot
# denoise_info$de_noise_info_plot
# denoise_info$threshold_value
# 
# 


denoise_info_10886 = GetDenoiseInfo(
  denoising_detail = ms2_denoising_info, 
  denoising_method = 'spline_regression', 
  ms2_spectrum_transform_method = 'log2_transform', 
  ms_id = 'controllerType=0 controllerNumber=1 scan=10886', 
  ms_data = Spectra::Spectra("./raw_data/210818_Lumos__LC1200__nanoFlex_NH21036_01_02_Ogly.mzML", source = MsBackendMzR()))


pdf(file = '210818_Lumos__LC1200__nanoFlex_NH21036_01_02_Ogly_Hex1HexNAc1Neu5Ac1_mass_spectra_denoised_plot.pdf', width = 18, height = 6)
denoise_info_10886$mass_spectra_denoised_plot
dev.off()
dev.off()

pdf(file = '210818_Lumos__LC1200__nanoFlex_NH21036_01_02_Ogly_Hex1HexNAc1Neu5Ac1_de_noise_info_plot.pdf', width = 6, height = 6)
denoise_info_10886$de_noise_info_plot
dev.off()
dev.off()




denoise_info_10471 = GetDenoiseInfo(
  denoising_detail = ms2_denoising_info, 
  denoising_method = 'spline_regression', 
  ms2_spectrum_transform_method = 'log2_transform', 
  ms_id = 'controllerType=0 controllerNumber=1 scan=10471', 
  ms_data = Spectra::Spectra("./raw_data/210818_Lumos__LC1200__nanoFlex_NH21036_01_03_Ogly.mzML", source = MsBackendMzR()))

pdf(file = '210818_Lumos__LC1200__nanoFlex_NH21036_01_03_Ogly_Hex1HexNAc2_mass_spectra_denoised_plot.pdf', width = 20, height = 6)
denoise_info_10471$mass_spectra_denoised_plot
dev.off()
dev.off()





denoise_info_9226 = GetDenoiseInfo(
  denoising_detail = ms2_denoising_info, 
  denoising_method = 'spline_regression', 
  ms2_spectrum_transform_method = 'log2_transform', 
  ms_id = 'controllerType=0 controllerNumber=1 scan=9226', 
  ms_data = Spectra::Spectra("./raw_data/210818_Lumos__LC1200__nanoFlex_NH21036_02_06_Ogly.mzML", source = MsBackendMzR()))


pdf(file = '210818_Lumos__LC1200__nanoFlex_NH21036_02_06_Ogly_Hex3_mass_spectra_denoised_plot.pdf', width = 20, height = 6)
denoise_info_9226$mass_spectra_denoised_plot
dev.off()
dev.off()






denoise_info_9876 = GetDenoiseInfo(
  denoising_detail = ms2_denoising_info, 
  denoising_method = 'spline_regression', 
  ms2_spectrum_transform_method = 'log2_transform', 
  ms_id = 'controllerType=0 controllerNumber=1 scan=9876', 
  ms_data = Spectra::Spectra("./raw_data/210818_Lumos__LC1200__nanoFlex_NH21036_01_02_Ogly.mzML", source = MsBackendMzR()))

pdf(file = '210818_Lumos__LC1200__nanoFlex_NH21036_01_02_Ogly_not_sure_comp_mass_spectra_denoised_plot_scan=9876.pdf', width = 20, height = 6)
denoise_info_9876$mass_spectra_denoised_plot
dev.off()
dev.off()



denoise_info_10518 = GetDenoiseInfo(
  denoising_detail = ms2_denoising_info, 
  denoising_method = 'spline_regression', 
  ms2_spectrum_transform_method = 'log2_transform', 
  ms_id = 'controllerType=0 controllerNumber=1 scan=10518', 
  ms_data = Spectra::Spectra("./raw_data/210818_Lumos__LC1200__nanoFlex_NH21036_01_03_Ogly.mzML", source = MsBackendMzR()))

pdf(file = '210818_Lumos__LC1200__nanoFlex_NH21036_01_03_Ogly_not_sure_comp_mass_spectra_denoised_plot_scan=10518.pdf', width = 20, height = 6)
denoise_info_10518$mass_spectra_denoised_plot
dev.off()
dev.off()






ms_data = Spectra::Spectra("./raw_data/210818_Lumos__LC1200__nanoFlex_NH21036_01_03_Ogly.mzML", source = MsBackendMzR())
ms_id_all = ms_data[["spectrumId"]]
ms_idx = which(ms_id_all == 'controllerType=0 controllerNumber=1 scan=10471')

spectra_data_raw = as.data.frame(Spectra::peaksData(ms_data)[[ms_idx]])

spectra_distinct = dplyr::distinct(spectra_data_raw, intensity, .keep_all = TRUE)

















HaCaT_WT_all_file_names

#setdiff(table_s2_HaCaT_WT_composition, unique(HaCaT_WT_all_ideintified_comp))

identify_only = setdiff(unique(HaCaT_WT_all_ideintified_comp), table_s2_HaCaT_WT_composition)
identify_only

# new glycan and right assignment
"Hex1HexNAc2"   
"Hex3"   

# likely to be glycan but composition assigned wrong
"Pentose1Hex2HexNAc1dHex1" 
"Hex3dHex1"    




  

for (file_name in HaCaT_WT_all_file_names) {
  temp_result = all_results[[file_name]]$final_spectrum_matching_result
  
  temp_result_comp = dplyr::filter(temp_result, glycan_string == "Hex3dHex1" ) 
  ms2_id = temp_result_comp$ms2_spectrum_id
  
  if (dim(temp_result_comp)[1] >= 1) {
    print(file_name)
    print(ms2_id)
  }
  
}



temp_data = Spectra::Spectra('./raw_data/210818_Lumos__LC1200__nanoFlex_NH21036_02_06_Ogly.mzML', source = MsBackendMzR())


which(temp_data[["spectrumId"]] == 'controllerType=0 controllerNumber=1 scan=9226') 


all_peaks_data = Spectra::peaksData(temp_data)

specific_data = all_peaks_data[[9226]]
