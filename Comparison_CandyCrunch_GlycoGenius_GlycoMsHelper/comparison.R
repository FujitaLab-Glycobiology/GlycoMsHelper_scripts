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

setwd('D:/Paper/Manual_scripts/Research_paper/GlycoMSHelper/Comparison_candycrunch_glycogenius')


ground_truth_info = read_csv("../r_data/241114_SLC35A2_SeqTypsin_ground_truth_info.csv")
ground_truth_compositions = unique(ground_truth_info$glycan_string)


glycomshelper_info = read.xlsx('../r_data/241114_SLC35A2_SeqTypsin_ms2_spectrum_composition_info.xlsx')
glycomshelper_compositions = unique(glycomshelper_info$glycan_string)

  
candycrunch_compositions = c('Hex3HexNAc2dHex1', 'Hex3HexNAc3', 'Hex3HexNAc3dHex1', 'Hex5HexNAc3', 'Hex5HexNAc4dHex1')


convert_glycogenius <- function(x) {
  extract <- function(pattern, str) {
    m <- regmatches(str, regexpr(pattern, str))
    if (length(m) == 0 || nchar(m) == 0) return(0)
    as.integer(sub("[A-Z]", "", m))
  }
  
  hex    <- extract("H[0-9]+", x)
  hexnac <- extract("N[0-9]+", x)
  fuc    <- extract("F[0-9]+", x)
  neuac  <- extract("S[0-9]+", x)
  
  result <- ""
  if (hex    > 0) result <- paste0(result, "Hex",    hex)
  if (hexnac > 0) result <- paste0(result, "HexNAc", hexnac)
  if (fuc    > 0) result <- paste0(result, "dHex",   fuc)
  if (neuac  > 0) result <- paste0(result, "Neu5Ac", neuac)
  return(result)
}
glycogenius_raw = read_csv('../GlycoGenius/gg_denoise_test/260617_161516_glycan_abundance_table_compositions.csv')
glycogenius_info = glycogenius_raw[2:dim(glycogenius_raw)[1], ]
glycogenius_info = dplyr::mutate(glycogenius_info, sample_comp = sapply(Sample, convert_glycogenius))
glycogenius_compositions = unique(glycogenius_info$sample_comp)





fit <- euler(list(
  ground_truth = ground_truth_compositions,
  glycomshelper = glycomshelper_compositions, 
  glycogenius = glycogenius_compositions, 
  candycrunch = candycrunch_compositions), shape = "ellipse", control = list(extraopt = TRUE))


# 
# fit <- eulerr::euler(list(
#   ground_truth = ground_truth_compositions,
#   glycomshelper = glycomshelper_compositions,
#   candycrunch = candycrunch_compositions,
#   glycogenius = glycogenius_compositions
# ))


setdiff(glycogenius_compositions, ground_truth_compositions)

setdiff(glycomshelper_compositions, ground_truth_compositions)


intersect(setdiff(glycogenius_compositions, ground_truth_compositions), 
          setdiff(glycomshelper_compositions, ground_truth_compositions)
          )

setdiff(intersect(glycomshelper_compositions, ground_truth_compositions), 
        intersect(glycogenius_compositions, ground_truth_compositions)
          
)



pdf(file = 'compositions_venn.pdf',
    width = 7, height = 10)
plot(
  fit,
  fills =c("#F4A3A8", "#64B9E8", "#F1C76E", "#A6D7B8"), 
  alpha = 0.7,
  quantities = TRUE,
  labels = TRUE
)
dev.off()
dev.off()







