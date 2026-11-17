library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ecoSynthesis)
library(sf)
library(fasterize)
library(raster)
library(terra)

# ---- Create directories if they don't exist -------------------------------------------------

### Paths
xrefdata <- here::here("data/xrefdata")
gisdata <- here::here("data/gisdata")
out_dir <- here::here("data/gisout")

for (check_dir in c(xrefdata, gisdata, out_dir)) {
  if (!dir.exists(check_dir)) {
    dir.create(check_dir)
  }
}

### Custom functions
source(here::here("scripts", "functions.R"))

# function to create a raster layer for a given target_id (EFG) and save it to the output directory

make_target_raster_layer <- function(
  target_id,
  survey_raster,
  survey_xy,
  survey_slc,
  rez = 100,
  dir = NULL,
  print = FALSE,
  plot = FALSE,
  write = TRUE,
  overwrite = FALSE
) {
  if (sum(print, plot, write) == 0) {
    cli::cli_abort("At least one of print, plot, or write must be TRUE")
  }

  if (!target_id %in% survey_slc$EFG) {
    cli::cli_abort(sprintf("%s not found in survey_slc$EFG", target_id))
  }

  if (write) {
    if (!fs::dir_exists(dir)) {
      cli::cli_abort(sprintf("Directory %s does not exist", dir))
    }
    path <- here::here(out_dir, sprintf("%s-res%s.tiff", target_id, rez))
    cli::cli_alert(
      c("write = TRUE", "writing {path}")
    )
  }

  sum_membs <- survey_slc |>
    group_by(MAPUNIT) |>
    summarise(
      k = sum(if_else(EFG == target_id, memb, 0)),
      total = sum(memb),
      prop = k / total
    )

  ntveg_k <- survey_xy |> left_join(sum_membs, by = "MAPUNIT")
  k_raster <- fasterize(ntveg_k, survey_raster, field = "prop", fun = "sum")

  if (print) {
    print(sum(k_raster@data@values, na.rm = T))
  }

  if (plot) {
    plot(k_raster)
  }

  if (write) {
    if (file.exists(path) && !overwrite) {
      cli::cli_abort("{path} already exists, skipping write")
    } else if (file.exists(path) && overwrite) {
      cli::cli_inform("{path} already exists, overwriting")
    }
    writeRaster(k_raster, path, overwrite = overwrite)
  }
}

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

ntveg_memb |>
  filter(!is.na(memb), memb > 0) |>
  group_by(round, expert) |>
  summarise(
    `source ids` = n_distinct(MAPUNIT),
    `target ids` = n_distinct(EFG),
    .groups = "drop"
  )


debugonce(check_membership)

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
  nrow()

distinct(ntveg_memb_missing_efgs, EFG) |> unlist() |> gluedown::md_bullet()
distinct(ntveg_memb_missing_efgs, MAPUNIT) |> unlist() |> gluedown::md_bullet()
