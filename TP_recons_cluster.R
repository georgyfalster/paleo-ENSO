## ---------------------------------------------------------------------------
##
## Script name: hierarchical cluster analysis of different paleo-ENSO reconstructions
##
## Purpose of script: attempt to replicate Fig. 3 of https://doi.org/10.1002/wcc.70036 
##
## Script author: Georgina Falster
##
## Date create/updated: 2026-05-29
##
## Email: georgina.falster@adelaide.edu.au
##
## Citation:

## ----------------------------------------------------------------------------
##
## Reconstruction metadata distilled as follows:
## ## Author_Year_ProxyType_Target_TargetSource_ProxyDist_Nested_Method_TemporalAveraging_ProxySensitivity_ProxyPretreatment
##
## MP = multi-proxy = 4 or more different archive types
## CoralMP = 3 or more proxy types from coral
##
## target given as simply 'ENSO' if the reconstruction wasn't trained explicitly on an index
##
## nested here means proxy depth changes through time. So, this is applicable to DA
## 
## ----------------------------------------------------------------------------

# =============================================================================
# set display options
# =============================================================================

options(scipen = 10, digits = 4)

# =============================================================================
# load packages
# =============================================================================

library(magrittr)
library(tidyverse)
library(patchwork)
library(ncdf4)
library(dendextend)

# change as necessary

fpath = "TP_recons/"

# =============================================================================
# get all the data from the various repos
# =============================================================================

# -----------------------------------------------------------------------------
# first deal with records from Zenodo
# -----------------------------------------------------------------------------
# ---------------------
# Falster et al 2023
# ---------------------

# download zipped recons to your computer
download.file(url='https://zenodo.org/records/10951789/files/Falster_PWCreconstructions.zip?download=1',
              destfile = paste0(fpath, "falster2023.zip"), method='curl')

# unzip
unzip(paste0(fpath, "falster2023.zip"),
      exdir = fpath)

# load the data
Falster_2023_MP_dSLP_Multi_Global_Nested_Multi_calyr_MoistTemp_raw <- read_csv(paste0(fpath, "recons/Falster2023_PWC_reconstruction_full_ensemble_median_and_95pct_range.csv")) %>%
  select(val = `Full ensemble median`, year)

# ---------------------
# Zhu et al 2022
# ---------------------

# download zipped recons to your computer
download.file(url='https://zenodo.org/records/5893781/files/fzhu2e/paper-volcENSO-v1.1.zip?download=1',
              destfile = paste0(fpath, "zhu2022.zip"), method='curl')

# unzip
unzip(paste0(fpath, "zhu2022.zip"),
      exdir = fpath)

{
  dat <- nc_open(paste0(fpath, "fzhu2e-paper-volcENSO-19ef135/recons/recon_Corals.nc"))
  Zhu_2022_CoralMP_Nino34_HadCRUT_GlobalTropics_Nested_DA_DecJanFeb_Temp_various <- data.frame(year = ncvar_get(dat, "year"), val = ncvar_get(dat, "nino34")[,3])
  nc_close(dat)
  
  dat <- nc_open(paste0(fpath, "fzhu2e-paper-volcENSO-19ef135/recons/recon_Corals_Li13b6.nc"))
  Zhu_2022_TreeWidthCoralMP_Nino34_HadCRUT_Global_Nested_DA_DecJanFeb_MoistTemp_various <- data.frame(year = ncvar_get(dat, "year"), val = ncvar_get(dat, "nino34")[,3])
  nc_close(dat)
  
  dat <- nc_open(paste0(fpath, "fzhu2e-paper-volcENSO-19ef135/recons/recon_Li13b6.nc"))
  Zhu_2022_TreeWidth_Nino34_HadCRUT_NH_Nested_DA_DecJanFeb_MoistTemp_various <- data.frame(year = ncvar_get(dat, "year"), val = ncvar_get(dat, "nino34")[,3])
  nc_close(dat)
  
  dat <- nc_open(paste0(fpath, "fzhu2e-paper-volcENSO-19ef135/recons/recon_Ocn2kCorals_Li13b6.nc"))
  Zhu_2022_TreeWidthCoralMP_Nino34_HadCRUT_Global_Nested_DA_DecJanFeb_MoistTemp_variousOcn2k <- data.frame(year = ncvar_get(dat, "year"), val = ncvar_get(dat, "nino34")[,3])
  nc_close(dat)
}

rm(dat)

# -----------------------------------------------------------------------------
# # now all the recons on NOAA
# -----------------------------------------------------------------------------

# https://www.ncdc.noaa.gov/paleo/study/29050 
Datwyler_2019_MP_Nino34_ERSST_Global_Nested_PCA_calyr_MoistTemp_raw <- read.table(
  url("https://www1.ncdc.noaa.gov/pub/data/paleo/reconstructions/datwyler2020/datwyler2020enso.txt"),
  header = TRUE) %>%
  select(year = age_CE, val = ENSOrec)

# https://www.ncei.noaa.gov/access/paleo-search/study/6250
DArrigo_2005_TreeWidth_Nino34_Kaplan_SWUSMexico_Nested_PCR_DecJanFeb_Moist_InstrumentalPersistence <- read.table(
  url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/nino3_recon-noaa.txt"),
  header = TRUE) %>% 
  select(year = age, val = nino3)

# https://www.ncei.noaa.gov/access/paleo-search/study/8409
Braganza_2009_TreeWidthIce_ENSO_NA_PacificBasin_NonNested_PCA_calyr_MoistTemp_threeYrLowPass <- read.table(
  url("https://www.ncei.noaa.gov/pub/data/paleo/contributions_by_author/braganza2009/braganza2009enso-noaa.txt"),
  header = TRUE) %>% 
  select(year = age, val = ensor5)

Braganza_2009_TreeWidthIceCoralMP_ENSO_NA_PacificBasin_NonNested_PCA_calyr_MoistTemp_threeYrLowPass <- read.table(
  url("https://www.ncei.noaa.gov/pub/data/paleo/contributions_by_author/braganza2009/braganza2009enso-noaa.txt"),
  header = TRUE) %>% 
  select(year = age, val = ensor8) %>%
  na.omit()

# https://www.ncei.noaa.gov/access/paleo-search/study/13684
EmileGeay_2013_MP_Nino34_ERSST_GlobalTropics_Nested_RegEM_calyr_MoistTemp_raw <- read.table(
  url("https://www.ncei.noaa.gov/pub/data/paleo/contributions_by_author/emile-geay2012/emile-geay2013-noaa.txt"),header = TRUE) %>% 
  select(year = age, val = ERSSTv3tot)

EmileGeay_2013_MP_Nino34_HadISST_GlobalTropics_Nested_RegEM_calyr_MoistTemp_raw <- read.table(
  url("https://www.ncei.noaa.gov/pub/data/paleo/contributions_by_author/emile-geay2012/emile-geay2013-noaa.txt"),header = TRUE) %>% 
  select(year = age, val = HadSST2itot)

EmileGeay_2013_MP_Nino34_Kaplan_GlobalTropics_Nested_RegEM_calyr_MoistTemp_raw <- read.table(
  url("https://www.ncei.noaa.gov/pub/data/paleo/contributions_by_author/emile-geay2012/emile-geay2013-noaa.txt"),header = TRUE) %>% 
  select(year = age, val = Kaplantot)

# https://www.ncdc.noaa.gov/paleo/study/26270
Freund_2019_CoralSrCorald18O_NinoWarmPool_HadISST_tropicalIndoPacific_Nested_PCR_seasonal_Temp_detrended <- read.table(url("https://www1.ncdc.noaa.gov/pub/data/paleo/reconstructions/freund2019/freund2019enso.txt"),header = TRUE) %>%
  filter(Season == 12) %>%
  select(year = age_CE, val = NWP)

Freund_2019_CoralSrCorald18O_NinoColdTongue_HadISST_tropicalIndoPacific_Nested_PCR_seasonal_Temp_detrended <- read.table(url("https://www1.ncdc.noaa.gov/pub/data/paleo/reconstructions/freund2019/freund2019enso.txt"),header = TRUE) %>%
  filter(Season == 12) %>%
  select(year = age_CE, val = NCT)

# this transformation is minimally explained in the paper but I think is accurate
# Ideally these shouldn't  be included as they are more or less just the same data as the NCT and NWP recons
solve_N3_N4 <- function(Nct, Nwp) {
  
  a <- 2/5
  
  # assume a = 2/5
  denom <- 1 - a^2
  
  N3_a <- (Nct + a * Nwp) / denom
  N4_a <- (Nwp + a * Nct) / denom
  
  use_a <- (N3_a * N4_a) > 0
  
  N3 <- ifelse(use_a, N3_a, Nct)
  N4 <- ifelse(use_a, N4_a, Nwp)
  
  return(data.frame(N3 = N3, N4 = N4))
}

{
  temp <- solve_N3_N4(Freund_2019_CoralSrCorald18O_NinoWarmPool_HadISST_tropicalIndoPacific_Nested_PCR_seasonal_Temp_detrended$val,
                      Freund_2019_CoralSrCorald18O_NinoColdTongue_HadISST_tropicalIndoPacific_Nested_PCR_seasonal_Temp_detrended$val)
  
  temp$year <- Freund_2019_CoralSrCorald18O_NinoWarmPool_HadISST_tropicalIndoPacific_Nested_PCR_seasonal_Temp_detrended$year
  
  Freund_2019_CoralSrCorald18O_Nino3_HadISST_tropicalIndoPacific_Nested_PCR_seasonal_Temp_detrended <- data.frame(year = temp$year, val = temp$N3)
  
  Freund_2019_CoralSrCorald18O_Nino4_HadISST_tropicalIndoPacific_Nested_PCR_seasonal_Temp_detrended <- data.frame(year = temp$year, val = temp$N4)
  
  rm(temp)
}

# https://www.ncei.noaa.gov/access/paleo-search/study/11194
Li_2011_TreeWidth_ENSO_NA_NorthAmerica_Nested_PCA_JunJulAug_Moist_PDSIReconHighPassFilter <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/enso-li2011-noaa.txt"),header = TRUE) %>%
  select(year = age, val = ensoi)

# https://www.ncei.noaa.gov/access/paleo-search/study/14632
Li_2013_TreeWidth_Nino34_Kaplan_CircumPacific_Nested_PCA_NovDecJan_Moist_DroughtAtlases <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/enso-li2013.txt"))[-1, ] %>%
  rename(year = V1, val = V2) %>%
  mutate_all(as.numeric)

# https://www.ncdc.noaa.gov/paleo/study/28417
Liu_2017_Treed18O_Nino4_Kaplan_Taiwan_Nested_LinReg_MarAprMay_Moist_raw <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/isotope/asia/china/taiwan2017d18o.txt"),header = TRUE) %>%
  select(year = age_CE, val = SSTA)

# https://www.ncei.noaa.gov/access/paleo-search/study/16459
Stahle_1993_TreeWidth_SOI_NA_SouthernNorthAmerica_Nested_LinReg_DecJanFeb_Moist_InstrumentalPersistence <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/southusa/soi.recon"),skip = 76) %>% 
  select(year = V2, val = V3)

# https://www.ncei.noaa.gov/access/paleo-search/study/6238
Stahle_1998_TreeWidth_SOI_NA_CircumPacific_Nested_LinReg_DecJanFeb_Moist_InstrumentalPersistence <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/soi_recon-noaa.txt"),header = TRUE) %>%
  select(year = age, val = soi)

# -----------------------------------------------------------------------------
# IMPORTANT NOTE: data in these recons is instrumental from 1979. 
# method and season details assumed similar to D'Arrigo but not stated on NOAA repo
# see also p67 of Wilson 2010
# -----------------------------------------------------------------------------

# https://www.ncei.noaa.gov/access/paleo-search/study/8704
Cook_2008_TreeWidth_Nino12_NA_SWUSMexico_Nested_PCR_DecJanFeb_Moist_InstrumentalPersistence <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/nino-cook2008-noaa.txt"),header = TRUE) %>%
  select(year = age, val = nino1.2)

Cook_2008_TreeWidth_Nino3_NA_SWUSMexico_Nested_PCR_DecJanFeb_Moist_InstrumentalPersistence <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/nino-cook2008-noaa.txt"),header = TRUE) %>%
  select(year = age, val = nino3)

Cook_2008_TreeWidth_Nino34_NA_SWUSMexico_Nested_PCR_DecJanFeb_Moist_InstrumentalPersistence <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/nino-cook2008-noaa.txt"),header = TRUE) %>%
  select(year = age, val = nino3.4)

Cook_2008_TreeWidth_Nino4_NA_SWUSMexico_Nested_PCR_DecJanFeb_Moist_InstrumentalPersistence <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/nino-cook2008-noaa.txt"),header = TRUE) %>%
  select(year = age, val = nino4)

# https://www.ncei.noaa.gov/access/paleo-search/study/8732
McGregor_2010_MP_ENSO_NA_Global_Nested_PCA_calyr_MoistTemp_various <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/contributions_by_author/mcgregor2010/mcgregor2010uep-noaa.txt"),header = TRUE) %>%
  select(year = age, val = uep)

# https://www.ncei.noaa.gov/access/paleo-search/study/11749
Wilson_2010_IceCoralMP_Nino34_HadISST_GlobalTropics_Nested_CPR_DecToNov_MoistTemp_various <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/enso-wilson2010-tel-noaa.txt"),header = TRUE) %>% 
  select(year = age, val = cprtel)

Wilson_2010_IceCoralMP_Nino34_HadISST_GlobalTropics_Nested_PCR_DecToNov_MoistTemp_various <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/enso-wilson2010-tel-noaa.txt"),header = TRUE) %>% 
  select(year = age, val = pcrtel)

Wilson_2010_IceCoralMP_Nino34_HadISST_GlobalTropics_Nested_RegEM_DecToNov_MoistTemp_various <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/enso-wilson2010-tel-noaa.txt"),header = TRUE) %>% 
  select(year = age, val = regemtel)

Wilson_2010_Corald18O_Nino34_HadISST_circumPacific_Nested_CPR_DecToNov_Temp_detrended <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/enso-wilson2010-coa-noaa.txt"),header = TRUE) %>%
  select(year = age, val = cprcoa)

Wilson_2010_Corald18O_Nino34_HadISST_circumPacific_Nested_PCR_DecToNov_Temp_detrended <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/enso-wilson2010-coa-noaa.txt"),header = TRUE) %>%
  select(year = age, val = pcrcoa)

Wilson_2010_Corald18O_Nino34_HadISST_circumPacific_Nested_RegEM_DecToNov_Temp_detrended <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/enso-wilson2010-coa-noaa.txt"),header = TRUE) %>%
  select(year = age, val = regemcoa)

# https://www.ncei.noaa.gov/access/paleo-search/study/17955
Tierney_2015_CoralMP_easternPacific_HadISST_GlobalTropics_Nested_CPS_AprToMar_Temp_raw <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/pages2k/tierney2015sst-epacific-noaa.txt"),header = TRUE) %>%
  select(year = age_CE, val = sst.apr.mar.anom)

# https://www.ncdc.noaa.gov/paleo/study/26710
Torbenson_2019_TreeWidth_MEI_HadCRUT_NorthAmerica_Nested_PCA_NovToFeb_MoistTemp_OnlyStableCors <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/northamerica/texas-mexico2019mei.txt"),header = TRUE) %>% 
  select(year = age_calCE, val = eMEIstable)

Torbenson_2019_TreeWidth_MEI_HadCRUT_NorthAmerica_Nested_PCA_NovToFeb_MoistTemp_AllRecs <- read.table(url("https://www.ncei.noaa.gov/pub/data/paleo/treering/reconstructions/northamerica/texas-mexico2019mei.txt"),header = TRUE) %>% 
  select(year = age_calCE, val = eMEIall)

# =============================================================================
# Liu et al 2024
# =============================================================================

# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023JD040491

# -----------------------------------------------------------------------------
# This is the sole reconstruction that is not publicly available.
# to run the code from here you will need to contact the authors 
# to request the data. Or see at the base of this script for a version
# that does not include the two reconstructions from this paper.
# -----------------------------------------------------------------------------

Liu_2024_Treed18OSpeld18O_Nino34_HadISST_Global_Nested_PCR_calyr_MoistTemp_nineYrHighPass <- read.csv(paste0(fpath, "Liu2024JGRA-ENSO-reconstruction.csv"), skip = 5) %>%
  select(year, val = PCR.reconstruction)

Liu_2024_Treed18OSpeld18O_Nino34_HadISST_Global_Nested_PCA_calyr_MoistTemp_nineYrHighPass <- read.csv(paste0(fpath, "Liu2024JGRA-ENSO-reconstruction.csv"), skip = 5) %>%
  select(year, val = DCC.reconstruction)

# =============================================================================
# get some metadata
# =============================================================================

recon_names <- c("Braganza_2009_TreeWidthIce_ENSO_NA_PacificBasin_NonNested_PCA_calyr_MoistTemp_threeYrLowPass",
                 "Braganza_2009_TreeWidthIceCoralMP_ENSO_NA_PacificBasin_NonNested_PCA_calyr_MoistTemp_threeYrLowPass",
                 "Datwyler_2019_MP_Nino34_ERSST_Global_Nested_PCA_calyr_MoistTemp_raw",
                 "DArrigo_2005_TreeWidth_Nino34_Kaplan_SWUSMexico_Nested_PCR_DecJanFeb_Moist_InstrumentalPersistence",
                 "EmileGeay_2013_MP_Nino34_ERSST_GlobalTropics_Nested_RegEM_calyr_MoistTemp_raw",
                 "EmileGeay_2013_MP_Nino34_HadISST_GlobalTropics_Nested_RegEM_calyr_MoistTemp_raw",
                 "EmileGeay_2013_MP_Nino34_Kaplan_GlobalTropics_Nested_RegEM_calyr_MoistTemp_raw",
                 "Falster_2023_MP_dSLP_Multi_Global_Nested_Multi_calyr_MoistTemp_raw",
                 "Freund_2019_CoralSrCorald18O_NinoWarmPool_HadISST_tropicalIndoPacific_Nested_PCR_seasonal_Temp_detrended",
                 "Freund_2019_CoralSrCorald18O_NinoColdTongue_HadISST_tropicalIndoPacific_Nested_PCR_seasonal_Temp_detrended",
                 "Freund_2019_CoralSrCorald18O_Nino3_HadISST_tropicalIndoPacific_Nested_PCR_seasonal_Temp_detrended",
                 "Freund_2019_CoralSrCorald18O_Nino4_HadISST_tropicalIndoPacific_Nested_PCR_seasonal_Temp_detrended",
                 "Li_2011_TreeWidth_ENSO_NA_NorthAmerica_Nested_PCA_JunJulAug_Moist_PDSIReconHighPassFilter",
                 "Li_2013_TreeWidth_Nino34_Kaplan_CircumPacific_Nested_PCA_NovDecJan_Moist_DroughtAtlases",
                 "Liu_2024_Treed18OSpeld18O_Nino34_HadISST_Global_Nested_PCR_calyr_MoistTemp_nineYrHighPass",
                 "Liu_2024_Treed18OSpeld18O_Nino34_HadISST_Global_Nested_PCA_calyr_MoistTemp_nineYrHighPass",
                 "Stahle_1993_TreeWidth_SOI_NA_SouthernNorthAmerica_Nested_LinReg_DecJanFeb_Moist_InstrumentalPersistence",
                 "Stahle_1998_TreeWidth_SOI_NA_CircumPacific_Nested_LinReg_DecJanFeb_Moist_InstrumentalPersistence",
                 "Cook_2008_TreeWidth_Nino12_NA_SWUSMexico_Nested_PCR_DecJanFeb_Moist_InstrumentalPersistence",
                 "Cook_2008_TreeWidth_Nino3_NA_SWUSMexico_Nested_PCR_DecJanFeb_Moist_InstrumentalPersistence",
                 "Cook_2008_TreeWidth_Nino34_NA_SWUSMexico_Nested_PCR_DecJanFeb_Moist_InstrumentalPersistence",
                 "Cook_2008_TreeWidth_Nino4_NA_SWUSMexico_Nested_PCR_DecJanFeb_Moist_InstrumentalPersistence",
                 "McGregor_2010_MP_ENSO_NA_Global_Nested_PCA_calyr_MoistTemp_various",
                 "Wilson_2010_IceCoralMP_Nino34_HadISST_GlobalTropics_Nested_CPR_DecToNov_MoistTemp_various",
                 "Wilson_2010_IceCoralMP_Nino34_HadISST_GlobalTropics_Nested_PCR_DecToNov_MoistTemp_various",
                 "Wilson_2010_IceCoralMP_Nino34_HadISST_GlobalTropics_Nested_RegEM_DecToNov_MoistTemp_various",
                 "Wilson_2010_Corald18O_Nino34_HadISST_circumPacific_Nested_CPR_DecToNov_Temp_detrended",
                 "Wilson_2010_Corald18O_Nino34_HadISST_circumPacific_Nested_PCR_DecToNov_Temp_detrended",
                 "Wilson_2010_Corald18O_Nino34_HadISST_circumPacific_Nested_RegEM_DecToNov_Temp_detrended",
                 "Tierney_2015_CoralMP_easternPacific_HadISST_GlobalTropics_Nested_CPS_AprToMar_Temp_raw",
                 "Torbenson_2019_TreeWidth_MEI_HadCRUT_NorthAmerica_Nested_PCA_NovToFeb_MoistTemp_OnlyStableCors",
                 "Torbenson_2019_TreeWidth_MEI_HadCRUT_NorthAmerica_Nested_PCA_NovToFeb_MoistTemp_AllRecs",
                 "Zhu_2022_CoralMP_Nino34_HadCRUT_GlobalTropics_Nested_DA_DecJanFeb_Temp_various",
                 "Zhu_2022_TreeWidthCoralMP_Nino34_HadCRUT_Global_Nested_DA_DecJanFeb_MoistTemp_various",
                 "Zhu_2022_TreeWidth_Nino34_HadCRUT_NH_Nested_DA_DecJanFeb_MoistTemp_various",
                 "Zhu_2022_TreeWidthCoralMP_Nino34_HadCRUT_Global_Nested_DA_DecJanFeb_MoistTemp_variousOcn2k"
)

recon_df <- matrix(NA, length(recon_names), 15) %>%
  set_colnames(c("Author","Year","ProxyType","Target","TargetSource","ProxyDist",
                 "Nested","Method","TemporalAveraging","ProxySensitivity",
                 "ProxyPretreatment","start", "end", "length", "res"))

for(i in 1:(length(recon_names))) {
  
  this_recon <- get(recon_names[i])
  
  recon_df[i, 1] <- strsplit(recon_names[i], "_")[[1]][1]
  recon_df[i, 2] <- strsplit(recon_names[i], "_")[[1]][2]
  recon_df[i, 3] <- strsplit(recon_names[i], "_")[[1]][3]
  recon_df[i, 4] <- strsplit(recon_names[i], "_")[[1]][4]
  recon_df[i, 5] <- strsplit(recon_names[i], "_")[[1]][5]
  recon_df[i, 6] <- strsplit(recon_names[i], "_")[[1]][6]
  recon_df[i, 7] <- strsplit(recon_names[i], "_")[[1]][7]
  recon_df[i, 8] <- strsplit(recon_names[i], "_")[[1]][8]
  recon_df[i, 9] <- strsplit(recon_names[i], "_")[[1]][9]
  recon_df[i, 10] <- strsplit(recon_names[i], "_")[[1]][10]
  recon_df[i, 11] <- strsplit(recon_names[i], "_")[[1]][11]
  recon_df[i, 12] <- min(this_recon$year)
  recon_df[i, 13] <- max(this_recon$year)
  recon_df[i, 14] <- length(na.omit((this_recon$val)))
  recon_df[i, 15] <- (max(this_recon$year)-min(this_recon$year))/length(na.omit((this_recon$val)))
  
  rm(this_recon)
  
}

recon_df <- as.data.frame(recon_df) %>%
  mutate_at(c("start", "end", "length", "res"), as.numeric) %>%
  mutate(recon_num = 1:n()) %>%
  mutate(recon_name = paste(recon_num, Author, Target, Method, ProxyType, sep = "_")) %>%
  # mutate(recon_name = paste(recon_num, ProxyType, Target, Method, ProxyDist, sep = "_")) %>%
  as_tibble()

# =============================================================================
# match analysis period of Freund
# =============================================================================

max_year <- 1900
min_year <- 1727

# ---------------------
# # indices for the PWC (dSLP, SOI) are inversely correlated with Nino 3.4 simply by construction.
# ---------------------

ref_ts <- filter(EmileGeay_2013_MP_Nino34_ERSST_GlobalTropics_Nested_RegEM_calyr_MoistTemp_raw,
                 year >= 1900 & year <= 1920)$val

all_recons <- matrix(NA, length(min_year:max_year), length(recon_names)) %>%
  set_colnames(recon_df$recon_name)

for(i in 1:(length(recon_names))) {
  
  this_recon <- get(recon_names[i])
  
  # first, align if necessary
  check_ts <- filter(this_recon,
                     year >= 1900 & year <= 1920)$val
  
  if(cor(ref_ts, check_ts) < 0) {
    
    this_recon$val <- this_recon$val*-1
    print(paste0("flipped ", recon_names[i]))
  }
  
  this_recon <- filter(this_recon, year <= max_year & year >= min_year)
  
  all_recons[ ,i] <- this_recon$val
  
  rm(this_recon)
}

# =============================================================================
# make a df
# =============================================================================

recs_df <- as.data.frame(all_recons) %>% mutate(year = min_year:max_year)

recs_df_long <- recs_df %>% pivot_longer(-year, names_to = "recon", values_to = "val")

# =============================================================================
# now match the 'hierarchical agglomerative clustering with the average linkage method'
# =============================================================================

# ---------------------
# preprocessing
# ---------------------

dat_scaled <- scale(all_recons)

# ---------------------
# clustering
# ---------------------

{
  dat_scaled_trans <- t(dat_scaled)
  d <- dist(dat_scaled_trans, method = "euclidean")
  hc1 <- hclust(d, method = "average" )
  hc1[["labels"]] <- gsub("^[^_]*_", "", hc1[["labels"]])
}

# ---------------------
# convert to dendrograms
# ---------------------

n_branch = 7

{
  dg_aligned <- as.dendrogram(hc1)
  
  dg_aligned <- color_branches(dg_aligned, k = n_branch)
  
  ggd_aligned <- as.ggdend(dg_aligned)
  
}

# =============================================================================
# plot
# =============================================================================

label_offset <- 0.02 * max(ggd_aligned$segments$y)

ggplot() +
  geom_segment(data = mutate(ggd_aligned$segments, col = ifelse(is.na(col), "black", col)),
    aes(x = x, y = y, xend = xend, yend = yend, colour = col),linewidth = 1) +
  geom_text(data = ggd_aligned$labels,aes(x = x, y = y-label_offset, label = label),
    size = 3.5,angle = 0,hjust = 0,vjust = 0.5) +
  labs(title = "Dendrogram summarising hierarchical relationships between different reconstructions of the\nEl Niño Southern Oscillation or Pacific Walker Circulation",
       subtitle = paste0("Calculated across ", min_year, "–", max_year, " CE")) +
  coord_flip(clip = "off") +
  scale_y_reverse() +
  scale_colour_identity() +
  # add author details
  annotate("text", x = -Inf, y = -Inf,
           label = paste0("Created by: Georgy Falster, ", format(Sys.Date(), "%d %b %Y")),
           hjust = -0.5, vjust = 01, size = 3, colour = "grey60") +
  theme_void() +
  theme(plot.margin = margin(t = 5, r = 250, b = 10, l = 5),
        axis.text = element_blank(),
        axis.title = element_blank())

ggsave(paste0("tropicalPacific_reconstructions_", min_year, "–", max_year, "hclust_dendrogram_includingLiu2024.png"),
       bg = "white", height = 17, width = 25, units = "cm")

# Inclusion of the 'transformed' Freund recons does not affect the dendrogram.
# They just tack on that end branch with the NCT and NWP recons

# =============================================================================
# Now a version only including publicaly available datasets 
# =============================================================================

# assumes you have run the line above which creates the 'recon_names' object

recons_pub <- recon_names[-which(grepl("Liu_2024", recon_names))]

# ---------------------
# get the details
# ---------------------

reconPub_df <- matrix(NA, length(recons_pub), 15) %>%
  set_colnames(c("Author","Year","ProxyType","Target","TargetSource","ProxyDist",
                 "Nested","Method","TemporalAveraging","ProxySensitivity",
                 "ProxyPretreatment","start", "end", "length", "res"))

for(i in 1:(length(recons_pub))) {
  
  this_recon <- get(recons_pub[i])
  
  reconPub_df[i, 1] <- strsplit(recons_pub[i], "_")[[1]][1]
  reconPub_df[i, 2] <- strsplit(recons_pub[i], "_")[[1]][2]
  reconPub_df[i, 3] <- strsplit(recons_pub[i], "_")[[1]][3]
  reconPub_df[i, 4] <- strsplit(recons_pub[i], "_")[[1]][4]
  reconPub_df[i, 5] <- strsplit(recons_pub[i], "_")[[1]][5]
  reconPub_df[i, 6] <- strsplit(recons_pub[i], "_")[[1]][6]
  reconPub_df[i, 7] <- strsplit(recons_pub[i], "_")[[1]][7]
  reconPub_df[i, 8] <- strsplit(recons_pub[i], "_")[[1]][8]
  reconPub_df[i, 9] <- strsplit(recons_pub[i], "_")[[1]][9]
  reconPub_df[i, 10] <- strsplit(recons_pub[i], "_")[[1]][10]
  reconPub_df[i, 11] <- strsplit(recons_pub[i], "_")[[1]][11]
  reconPub_df[i, 12] <- min(this_recon$year)
  reconPub_df[i, 13] <- max(this_recon$year)
  reconPub_df[i, 14] <- length(na.omit((this_recon$val)))
  reconPub_df[i, 15] <- (max(this_recon$year)-min(this_recon$year))/length(na.omit((this_recon$val)))
  
  rm(this_recon)
  
}

reconPub_df <- as.data.frame(reconPub_df) %>%
  mutate_at(c("start", "end", "length", "res"), as.numeric) %>%
  mutate(recon_num = 1:n()) %>%
  mutate(recon_name = paste(recon_num, Author, Target, Method, ProxyType, sep = "_")) %>%
  # mutate(recon_name = paste(recon_num, ProxyType, Target, Method, ProxyDist, sep = "_")) %>%
  as_tibble()

# -----------------------------------------------------------------------------
# all the other steps
# -----------------------------------------------------------------------------

ref_ts <- filter(EmileGeay_2013_MP_Nino34_ERSST_GlobalTropics_Nested_RegEM_calyr_MoistTemp_raw,
                 year >= 1900 & year <= 1920)$val

all_reconsPub <- matrix(NA, length(min_year:max_year), length(recons_pub)) %>%
  set_colnames(reconPub_df$recon_name)

for(i in 1:(length(recons_pub))) {
  
  this_recon <- get(recons_pub[i])
  
  # first, align if necessary
  check_ts <- filter(this_recon,
                     year >= 1900 & year <= 1920)$val
  
  if(cor(ref_ts, check_ts) < 0) {
    
    this_recon$val <- this_recon$val*-1
    print(paste0("flipped ", recons_pub[i]))
  }
  
  this_recon <- filter(this_recon, year <= max_year & year >= min_year)
  
  all_reconsPub[ ,i] <- this_recon$val
  
  rm(this_recon)
}


recsPub_df <- as.data.frame(all_reconsPub) %>% mutate(year = min_year:max_year)

recsPub_df_long <- recsPub_df %>% pivot_longer(-year, names_to = "recon", values_to = "val")

# ---------------------
# preprocessing
# ---------------------

datPub_scaled <- scale(all_reconsPub)

# ---------------------
# clustering
# ---------------------

{
  dat_scaled_trans <- t(datPub_scaled)
  d <- dist(dat_scaled_trans, method = "euclidean")
  hc1p <- hclust(d, method = "average" )
  hc1p[["labels"]] <- gsub("^[^_]*_", "", hc1p[["labels"]])
}

# ---------------------
# convert to dendrograms
# ---------------------

n_branch = 7

{
  dgp_aligned <- as.dendrogram(hc1p)
  
  dgp_aligned <- color_branches(dgp_aligned, k = n_branch)
  
  ggdp_aligned <- as.ggdend(dgp_aligned)
  
}

# ---------------------
# the plot
# ---------------------

ggplot() +
  geom_segment(data = mutate(ggdp_aligned$segments, col = ifelse(is.na(col), "black", col)),
               aes(x = x, y = y, xend = xend, yend = yend, colour = col),linewidth = 1) +
  geom_text(data = ggdp_aligned$labels,aes(x = x, y = y-label_offset, label = label),
            size = 3.5,angle = 0,hjust = 0,vjust = 0.5) +
  labs(title = "Dendrogram summarising hierarchical relationships between different reconstructions of the\nEl Niño Southern Oscillation or Pacific Walker Circulation",
       subtitle = paste0("Calculated across ", min_year, "–", max_year, " CE, using only publicly available reconstructions")) +
  coord_flip(clip = "off") +
  scale_y_reverse() +
  scale_colour_identity() +
  # add author details
  annotate("text", x = -Inf, y = -Inf,
           label = paste0("Created by: Georgy Falster, ", format(Sys.Date(), "%d %b %Y")),
           hjust = -0.5, vjust = 01, size = 3, colour = "grey60") +
  theme_void() +
  theme(plot.margin = margin(t = 5, r = 250, b = 10, l = 5),
        axis.text = element_blank(),
        axis.title = element_blank())

ggsave(paste0("tropicalPacific_reconstructions_", min_year, "–", max_year, "hclust_dendrogram_excludingLiu2024.png"),
       bg = "white", height = 17, width = 25, units = "cm")

# the results are apparently quite stable! 
# removing two unpublished timeseries otherwise does not alter the dendrogram