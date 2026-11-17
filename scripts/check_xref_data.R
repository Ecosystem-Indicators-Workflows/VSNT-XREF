library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ecoSynthesis)
library(sf)
library(fasterize)
library(raster)
library(terra)
# library(pointblank)

# ---- Create directories if not present -------------------------------------------------

### Paths
xrefdata <- here::here("data/xrefdata")
gisdata <- here::here("data/gisdata")
out_dir <- here::here("data/gisout")

for (check_dir in c(xrefdata, gisdata, out_dir)) {
  if (!dir.exists(check_dir)) {
    dir.create(check_dir)
  }
}

### Predefined parameters

the_round <- 3

### Custom functions
source(here::here("scripts", "functions.R"))

# ---- Download and unzip the NT NVIS & vegetation data ------------------------------------------------

list(
  nvis_data_url = "https://ftp-dlrm.nt.gov.au/main.html?download&weblink=0174c920a0189d8c13f18dd3121c9675&realfilename=VegResources_NVIS_1M.zip",
  veg_data_url = "https://ftp-dlrm.nt.gov.au/main.html?download&weblink=e78bd012813d86c84c37681c72021704&realfilename=VegResources_ntveg_1M.zip"
) |>
  purrr::walk(
    ~ get_data(
      data_url = .x,
      output_file = {
        here::here(gisdata, urltools::param_get(.x, "realfilename"))
      },
      gisdata_path = gisdata,
      overwrite = FALSE
    )
  )

NT_ntveg <- read_sf(here::here(gisdata, "ntveg_1M.shp"))

# ---- Create a lookup tables ------------------------------------------------

## Lookup table source classification - NT vegetation data

ntveg_LUT <- NT_ntveg |>
  st_drop_geometry() |>
  filter(!is.na(MAPUNIT)) |>
  dplyr::select(
    MAPUNIT,
    ENV_DESC,
    L5_ASSOC_1,
    L5_ASSOC_2,
    L5_ASSOC_3,
    L5_ASSOC_4,
    L5_ASSOC_5,
    L5_ASSOC_6
  ) |>
  unique() |>
  lookup_tbl("MAPUNIT", "ENV_DESC") |>
  mutate(MAPUNIT = as.numeric(MAPUNIT))

## Lookup table target classification

IUCN_GET <- lookup_tbl(IUCN_GET_LUT, "efg_id", "efg_name")

# ---- Load the crosswalk data --------------------------------------------

# read data from the xrefdata directory
read_xref_data <- function(file_name) {
  file_path <- here::here(xrefdata, file_name)
  if (!file.exists(file_path)) {
    cli::cli_abort(c(
      "File does not exist in xrefdata directory",
      file_name
    ))
  }
  readr::read_csv(file_path, show_col_types = FALSE, skip = 12)
}

ntveg_memb <-
  list(
    expert = c("MGT", "MGT", "MGT", "DAK"),
    round_nr = c(1, 2, 3, 3),
    file_name = c(
      "MGT_V1_ntveg_1M_290925.csv",
      "MGT_V2_ntveg_1M_151026.csv",
      "MGT_V3_ntveg_1M_091126.csv",
      "DAK_V1_ntveg_1M_271026.csv"
    )
  ) |>
  purrr::pmap_dfr(
    ~ read_xref_data(file_name = ..3) |>
      pivot_wide_membership_table(
        original_names = "MAPUNIT",
        target_names = sprintf("EFG_%s", 1:3),
        membership_values = sprintf("memb_%s", 1:3),
        expert_name = ..1,
        round_nr = ..2
      ) |>
      membership_tbl(
        source_id = "MAPUNIT",
        target_id = "EFG",
        membership = "memb",
        interpreter = "expert",
        round = "round"
      )
  ) |>
  arrange(MAPUNIT) |>
  filter(!is.na(EFG))

## Explore and validate cross-reference data

# checked <- check_membership_pb(
#   x = ntveg_memb,
#   source = ntveg_LUT,
#   target = IUCN_GET
# )

# checked |>
#   export_report(filename = here::here(out_dir, "xref_data_check_report.html"))

# checked |> get_sundered_data(type = "pass")
# checked |> get_sundered_data(type = "fail")

ntveg_memb |>
  filter(!is.na(memb), memb > 0) |>
  group_by(round, expert) |>
  summarise(
    `source ids` = n_distinct(MAPUNIT),
    `target ids` = n_distinct(EFG),
    .groups = "drop"
  )

# debugonce(check_membership)

check_membership(ntveg_memb, source = ntveg_LUT, target = IUCN_GET)

# Every cross-ref row can be linked to spatial feature?
# should be 0, if not then there are some MAPUNITs in the membership table that are not in the source lookup table
ntveg_memb |>
  anti_join(ntveg_LUT) |>
  nrow()

# Inversely, every spatial feature can be linked to cross-ref row?
# should be 0, if not then there are some MAPUNITs in the source lookup table that are not in the membership table
ntveg_LUT |>
  anti_join(ntveg_memb) |>
  nrow()

# Every cross-ref row can be linked to target lookup table?
# should be 0, if not then there are some EFGs in the membership table that are not in the target lookup table

ntveg_memb_missing_efgs <- ntveg_memb |>
  anti_join(IUCN_GET, by = c("EFG" = "efg_id"))

ntveg_memb_missing_efgs |>
  summarise(n_distinct(MAPUNIT), .by = c(expert, round, EFG))

ntveg_memb_missing_efgs |>
  nrow()

distinct(ntveg_memb_missing_efgs, EFG) |> unlist() |> gluedown::md_bullet()
distinct(ntveg_memb_missing_efgs, MAPUNIT) |> unlist() |> gluedown::md_bullet()

## Clean membership table by filtering out records with invalid membership values, source IDs, or target IDs

ntveg_memb_check <- clean_membership(
  x = ntveg_memb,
  source = ntveg_LUT,
  target = IUCN_GET,
  action = "check"
) # this matches what we expect from the previous checks, so we can proceed to clean the membership table

ntveg_memb_check |>
  drop_na(flag) |>
  summarise(
    source_units = n_distinct(MAPUNIT),
    .by = c(expert, round, flag)
  )

ntveg_memb_clean <- clean_membership(
  x = ntveg_memb,
  source = ntveg_LUT,
  target = IUCN_GET,
  action = "filter"
)

# ---- Synthesis -----

ntveg_best_summary <- summarise_membership(
  ntveg_memb_clean |> mutate(MAPUNIT = as.character(MAPUNIT)),
  keep_all = FALSE,
  which_round = the_round
)

ntveg_best_summary
