# ## this file cleans the ovarian cancer dataset
# # library(tidyverse)
# # library(VectraPolarisData) # Bioconductor data package
# #
# # # load spatial experiment object
# # oc <- HumanOvarianCancerVP()
# # ## Assays slots
# # assays_slot <- assays(oc)
# # intensities_df <- assays_slot$intensities
# # nucleus_intensities_df<- assays_slot$nucleus_intensities
# # rownames(nucleus_intensities_df) <- paste0("nucleus_", rownames(nucleus_intensities_df))
# # membrane_intensities_df<- assays_slot$membrane_intensities
# # rownames(membrane_intensities_df) <- paste0("membrane_", rownames(membrane_intensities_df))
# #
# # # colData and spatialData
# # colData_df <- colData(oc)
# # spatialCoords_df <- spatialCoords(oc)
# #
# # # clinical data
# # patient_level_ovarian <- metadata(oc)$clinical_data %>%
# #   # create binary stage variable
# #   dplyr::mutate(stage_bin = ifelse(stage %in% c("1", "2"), 0, 1))
# #
# # # only include samples for whom we have an id in the patient dataset, which is 128- the other 4 are controls
# # # give Markers shorter names
# # cell_level_ovarian <- as.data.frame(cbind(colData_df,
# #                                           spatialCoords_df,
# #                                           t(intensities_df),
# #                                           t(nucleus_intensities_df),
# #                                           t(membrane_intensities_df))
# # ) %>%
# #   dplyr::rename(cd19 = cd19_opal_480,
# #                 cd68 = cd68_opal_520,
# #                 cd3 = cd3_opal_540,
# #                 cd8 = cd8_opal_650,
# #                 ier3 = ier3_opal_620,
# #                 pstat3 = p_stat3_opal_570,
# #                 ck = ck_opal_780,
# #                 ki67 = ki67_opal_690) %>%
# #   # define cell type 'immune'
# #   mutate(immune = ifelse(phenotype_cd19 == "CD19+" | phenotype_cd8 == "CD8+" |
# #                            phenotype_cd3 == "CD3+" | phenotype_cd68 == "CD68+", "immune", "other")) %>%
# #   dplyr::select(contains("id"), tissue_category, contains("phenotype"),
# #                 contains("position"), ck:dapi, immune) %>%
# #   # only retain 128 subjects who have clinical data (other 4 are controls)
# #   dplyr::filter(sample_id %in% patient_level_ovarian$sample_id)
# #
# #
# # # data frame with clinical characteristics where each row is a different cell
# # ovarian_df <- full_join(patient_level_ovarian, cell_level_ovarian, by = "sample_id") %>%
# #   mutate(sample_id = as.numeric(as.factor(sample_id))) %>%
# #   filter(tissue_category != "Glass") %>%
# #   dplyr::select(sample_id, cell_id, tissue_category, x = cell_x_position, y = cell_y_position,
# #                 everything()) %>%
# #   select(-tma, -diagnosis, -grade)
# #
# # rm(oc, assays_slot, intensities_df, nucleus_intensities_df, membrane_intensities_df, colData_df, spatialCoords_df, patient_level_ovarian, cell_level_ovarian)
# #
#
#
# ##load processed ovarian cancer data
# load(url("https://github.com/julia-wrobel/MI_tutorial/raw/main/Data/ovarian.RDA"))
#
# # clean data
# ovarian_df_full = ovarian_df %>%
#   # subset to only analyze tumor areas
#   # provide more intuitive patient and image IDs
#   mutate(patient_id = as.numeric(factor(sample_id))) %>%
#   # define cell type 'immune', which groups all immune cells
#   mutate(immune = ifelse(phenotype_cd19 == "CD19+" | phenotype_cd8 == "CD8+" |
#                            phenotype_cd3 == "CD3+" | phenotype_cd68 == "CD68+", "immune", "other"),
#          phenotype = case_when(phenotype_cd19 == "CD19+" ~ "B-cell",
#                                phenotype_cd3 == "CD3+" ~ "T-cell",
#                                phenotype_cd68 == "CD68+" ~ "macrophage",
#                                TRUE ~ "other"),
#          phenotype = factor(phenotype),
#          x = x/1000,
#          y = y/1000) %>%
#   filter(tissue_category == "Tumor", immune == "immune") %>%
#   select(patient_id,x, y, age = age_at_diagnosis, immune, survival_time,
#          event = death, immune, stage = stage_bin)
#
# #
# # #
# # # enhance signal of original data for illustrating use of the package
# add_clustering <- function(mximg){
#   # create a convex hull as the observation window
#   w = spatstat.geom::convexhull.xy(mximg[["x"]], mximg[["y"]])
#
#   # create ppp object
#   pp_obj = spatstat.geom::ppp(mximg[["x"]], mximg[["y"]], window = w, checkdup = FALSE)
#
#   if(mximg[["cluster"]][1] == TRUE){
#     pp_obj_clust = spatstat.random::rMatClust(7, 1, 20, w)
#     pp_obj_clust2 = spatstat.random::rMatClust(4, 5, 15, w)
#     pp_obj_clust = spatstat.geom::superimpose(pp_obj_clust2,pp_obj_clust)
#   }else{
#     pp_obj_clust = spatstat.random::rpoispp(50, win = w)
#   }
#
#   pp_obj = spatstat.geom::superimpose(pp_obj,pp_obj_clust)
#
#   as_tibble(pp_obj) %>% mutate(immune = "immune", x = x * 1000, y = y * 1000)
# }
#
# df_nest = ovarian_df_full %>%
#   mutate(cluster = ifelse(survival_time > median(survival_time), TRUE, FALSE)) %>%
#   nest(data = c(x, y, immune, cluster))
#
# set.seed(1323)
# ovarian_df = df_nest %>% mutate(new_pp = map(df_nest$data, add_clustering)) %>%
#   select(-data) %>%
#   unnest(new_pp) %>%
#   distinct()
#
# ### Make mxFDA object
# clinical = ovarian_df %>%
#   select(patient_id, age, survival_time, event, stage) %>% distinct() %>%
#   mutate(sample_id = patient_id)
# spatial = ovarian_df %>%
#   select(-survival_time, -event, -age, -stage) %>%
#   rename("sample_id" = patient_id)
#
# ovarian_FDA = make_mxfda(clinical,
#                          spatial,
#                          subject_key = "patient_id",
#                          sample_key = "sample_id")
#
#
# # extract gfunctions from the ovarian data
# ovarian_FDA = extract_summary_functions(ovarian_FDA,
#                                         extract_func = extract_univariate,
#                                         summary_func = Gest,
#                                         r_vec = seq(0, 50, by = 1),
#                                         edge_correction = "rs",
#                                         markvar = "immune",
#                                         mark1 = "immune")
#
#
# ## write.csv(ovarian_gfun, file=gzfile("data-raw/ovarian_FDA.csv.gz", compression = 9))
# # #
# # # ovarian_gfun = read.csv(gzfile("data-raw/ovarian_FDA.csv.gz"))
#
# # add dataset to data folder
# usethis::use_data(ovarian_FDA, overwrite = TRUE)

