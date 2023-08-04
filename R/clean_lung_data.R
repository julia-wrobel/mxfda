## this file cleans the lung cancer dataset
library(tidyverse)
# library(VectraPolarisData) # Bioconductor data package
#
# # load lung data
# lung <- HumanLungCancerV3()
#
# ## Assays slots
# assays_slot <- assays(lung)
# intensities_df <- assays_slot$intensities
# nucleus_intensities_df<- assays_slot$nucleus_intensities
# rownames(nucleus_intensities_df) <- paste0("nucleus_", rownames(nucleus_intensities_df))
# membrane_intensities_df<- assays_slot$membrane_intensities
# rownames(membrane_intensities_df) <- paste0("membrane_", rownames(membrane_intensities_df))
#
# # colData and spatialData
# colData_df <- colData(lung)
# spatialCoords_df <- spatialCoords(lung)
#
# # clinical data
# patient_level_lung <- metadata(lung)$clinical_data
#
# cell_level_lung <- as_tibble(cbind(colData_df,
#                                    spatialCoords_df,
#                                    t(intensities_df),
#                                    t(nucleus_intensities_df),
#                                    t(membrane_intensities_df))
# )   %>%
#   dplyr::rename(cd19 = cd19_opal_650,
#                 cd3 = cd3_opal_520,
#                 cd14 = cd14_opal_540,
#                 cd8 = cd8_opal_620,
#                 hladr = hladr_opal_690,
#                 ck = ck_opal_570) %>%
#   dplyr::select(cell_id:slide_id, sample_id:dapi,
#                 entire_cell_axis_ratio:entire_cell_area_square_microns, contains("phenotype"))
#
# # data frame with clinical characteristics where each row is a different cell
# lung_df <- full_join(patient_level_lung, cell_level_lung, by = "slide_id") %>%
#   #mutate(slide_id = as.numeric(as.factor(slide_id))) %>%
#   dplyr::select(image_id = sample_id, patient_id = slide_id,
#                 cell_id, x = cell_x_position, y = cell_y_position,
#                 everything())
#
# rm(lung, assays_slot, intensities_df, nucleus_intensities_df, membrane_intensities_df, colData_df, spatialCoords_df, patient_level_lung,
#    cell_level_lung)
#

load(url("https://github.com/julia-wrobel/MI_tutorial/raw/main/Data/lung.RDA"))


# clean data
lung_df = lung_df %>%
  # subset to only analyze tumor areas
  filter(tissue_category == "Tumor") %>%
  # provide more intuitive patient and image IDs
  mutate(patient_id = as.numeric(factor(patient_id))) %>%
  group_by(patient_id) %>%
  mutate(image_id = as.numeric(factor(image_id))) %>%
  ungroup() %>%
  # create patient-image id
  mutate(patientImage_id = str_c(patient_id, "_", image_id)) %>%
  # define cell type 'immune', which groups all immune cells
  mutate(immune = ifelse(phenotype_cd19 == "CD19+" | phenotype_cd8 == "CD8+" |
                           phenotype_cd4 == "CD4+" | phenotype_cd14 == "CD14+", "immune", "other")) %>%
  select(contains("id"), x, y, gender, age = age_at_diagnosis, immune, survival_days, survival_status,
         contains("phenotype"), tissue_category, stage = stage_numeric)





# add dataset to data folder
usethis::use_data(lung_df, overwrite = TRUE)
