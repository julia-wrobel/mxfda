#creating the mxFDA object for the VectraPolarisData
#load processed ovarian cancer data
load(url("https://github.com/julia-wrobel/MI_tutorial/raw/main/Data/ovarian.RDA"))

#clean data
ovarian_df_full = ovarian_df %>%
  #subset to only analyze tumor areas
  #provide more intuitive patient and image IDs
  mutate(patient_id = as.numeric(factor(sample_id))) %>%
  #define cell type 'immune', which groups all immune cells
  mutate(immune = ifelse(phenotype_cd19 == "CD19+" | phenotype_cd8 == "CD8+" |
                           phenotype_cd3 == "CD3+" | phenotype_cd68 == "CD68+", "immune", "other"),
         phenotype = case_when(phenotype_cd19 == "CD19+" ~ "B-cell",
                               phenotype_cd3 == "CD3+" ~ "T-cell",
                               phenotype_cd68 == "CD68+" ~ "macrophage",
                               TRUE ~ "other"),
         phenotype = factor(phenotype),
         x = x/1000,
         y = y/1000) %>%
  filter(tissue_category == "Tumor", immune == "immune") %>%
  select(patient_id,x, y, age = age_at_diagnosis, immune, survival_time,
         event = death, immune, stage = stage_bin)


#enhance signal of original data for illustrating use of the package
add_clustering <- function(mximg){
  #create a convex hull as the observation window
  w = spatstat.geom::convexhull.xy(mximg[["x"]], mximg[["y"]])
  #create ppp object
  pp_obj = spatstat.geom::ppp(mximg[["x"]], mximg[["y"]], window = w, checkdup = FALSE)

  if(mximg[["cluster"]][1] == TRUE){
    pp_obj_clust = spatstat.random::rMatClust(7, 1, 20, w)
    pp_obj_clust2 = spatstat.random::rMatClust(4, 5, 15, w)
    pp_obj_clust = spatstat.geom::superimpose(pp_obj_clust2,pp_obj_clust)
  }else{
    pp_obj_clust = spatstat.random::rpoispp(50, win = w)
  }

  pp_obj = spatstat.geom::superimpose(pp_obj,pp_obj_clust)

  as_tibble(pp_obj) %>% mutate(immune = "immune", x = x * 1000, y = y * 1000)
}

df_nest = ovarian_df_full %>%
  mutate(cluster = ifelse(survival_time > median(survival_time), TRUE, FALSE)) %>%
  nest(data = c(x, y, immune, cluster))

set.seed(1323)
ovarian_df = df_nest %>% mutate(new_pp = map(df_nest$data, add_clustering)) %>%
  select(-data) %>%
  unnest(new_pp) %>%
  distinct()

#Make mxFDA object
clinical = ovarian_df %>%
  select(patient_id, age, survival_time, event, stage) %>% distinct() %>%
  mutate(sample_id = patient_id)
spatial = ovarian_df %>%
  select(-survival_time, -event, -age, -stage) %>%
  rename("sample_id" = patient_id)

ovarian_FDA = make_mxfda(clinical,
                         spatial,
                         subject_key = "patient_id",
                         sample_key = "sample_id")


#extract gfunctions from the ovarian data
ovarian_FDA = extract_summary_functions(ovarian_FDA,
                                        extract_func = univariate,
                                        summary_func = Gest,
                                        r_vec = seq(0, 50, by = 1),
                                        edge_correction = "rs",
                                        markvar = "immune",
                                        mark1 = "immune")

#if wanting to export for another reason
# write.csv(ovarian_gfun, file=gzfile("data-raw/ovarian_FDA.csv.gz", compression = 9))
# ovarian_gfun = read.csv(gzfile("data-raw/ovarian_FDA.csv.gz"))

#add dataset to data folder
usethis::use_data(ovarian_FDA, overwrite = TRUE)

