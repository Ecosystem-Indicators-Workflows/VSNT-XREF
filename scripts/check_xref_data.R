library(here)
library(readr)
library(dplyr)
library(tidyr)
library(ecoSynthesis)
library(sf)
library(fasterize)
library(raster)
library(terra)
library(fs)
library(glue)
library(geco)
library(stringr)
library(pointblank)

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

IUCN_GET <- tibble(
  efg_id = "T5.6",
  efg_name = "Semi-arid woodlands",
  biome_id = "T5",
  biome_name = "Deserts and semi-deserts biome",
  realm_id = "T",
  realm_name = "Terrestrial"
) |>
  bind_rows(IUCN_GET_LUT) |>
  lookup_tbl("efg_id", "efg_name")

# ---- Load the crosswalk data --------------------------------------------

# read data from the xrefdata directory
read_xref_data <- function(file_path) {
  if (!file.exists(file_path)) {
    cli::cli_abort(c(
      "File does not exist in xrefdata directory",
      file_path
    ))
  }
  readr::read_csv(file_path, show_col_types = FALSE, skip = 12)
}

## shift each interpreter's version up by the difference between the global max version and the max version for that interpreter, so that the latest version for each interpreter is aligned to the same round number
## pointblank validates that the version number for each interpreter is less than or equal to the global max version number without breaking the pipe
ntveg_memb <-
  tibble(file_path = fs::dir_ls(xrefdata)) |>
  mutate(
    expert = str_extract(fs::path_file(file_path), "^[^_]+"),
    version_nr = str_extract(fs::path_file(file_path), "(?<=_V)[0-9]+") |>
      as.numeric(),
    global_max_version_nr = max(version_nr)
  ) |>
  mutate(
    round_nr = version_nr + (global_max_version_nr - max(version_nr)),
    interp_min_round_nr = min(round_nr),
    .after = version_nr,
    .by = expert
    # .keep = "unused"
  ) |>
  pointblank::col_vals_lte(
    columns = round_nr,
    value = vars(global_max_version_nr),
    segment = vars(expert),
  ) |>
  pointblank::col_vals_gte(
    columns = version_nr,
    value = vars(round_nr),
    segment = vars(expert),
    preconditions = \(x) {
      filter(x, version_nr == last(version_nr), .by = expert)
    },
    actions = pointblank::action_levels(
      warn_at = 1,
      stop_at = NULL,
      fns = list(
        warn = ~ cli::cli_warn(
          "Version number for is less than the round number for the latest version for that interpreter. This may indicate that the crosswalk data is not up to date."
        )
      )
    )
  ) |>
  purrr::pmap_dfr(
    ~ read_xref_data(file_path = ..1) |>
      pivot_wide_membership_table(
        original_names = "MAPUNIT",
        target_names = sprintf("EFG_%s", 1:3),
        membership_values = sprintf("memb_%s", 1:3),
        expert_name = ..2,
        round_nr = ..3
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

# ----- Step 4 - read additional data for outputs -----

## 4a. Read GIS data
## Convert vector to raster using 'fasterize'

## template raster file
ntveg_xy <- sf::st_transform(NT_ntveg, "EPSG:3112") |>
  dplyr::filter(!is.na(MAPUNIT)) |>
  dplyr::mutate(MAPUNIT = as.character(MAPUNIT))

## 4b. Standard pixel values for Global Ecosystem Atlas + new T5.6

t5_6_info <- tibble(
  EFG = "T5.6",
  efg_name = "Semi-arid woodlands",
  efg_pixel_value = 35,
  efg_hex_code = "#EEB4B4",
  efg_rgb_code = "238,180,180"
) #rosybronw2

efg_info <- geco::geco_efg_data() |>
  dplyr::select(
    efg_code,
    efg_name,
    efg_pixel_value,
    efg_hex_code,
    efg_rgb_code
  ) |>
  dplyr::rename("EFG" = "efg_code") |>
  dplyr::bind_rows(t5_6_info) |>
  dplyr::arrange(efg_pixel_value)

efg_pixel_values <- efg_info |>
  dplyr::select(EFG, efg_pixel_value)

## Create a raster stack of all the EFG rasters

# ----- Step 5 - export outputs -----

xref_best_file <- "ntveg-xref-best-candidate.csv"
efg_lookup_file <- "IUCN-GET-pixel-lookup.csv"
xref_clr_file <- "ntveg-xref-target.clr"
xref_aux_file <- "ntveg-xref-target.tiff.aux.xml"
xref_qml_file <- "ntveg-xref-target.qml"
xref_tiff_file <- "ntveg-xref-target.tiff"
memb_tiff_file <- "ntveg-xref-membership.tiff"
supp_tiff_file <- "ntveg-xref-supporting.tiff"

## 5a. Export tabular synthesis tables
xref_best_tbl <- ntveg_LUT |>
  mutate(MAPUNIT = as.character(MAPUNIT)) |>
  left_join(ntveg_best_summary, by = c("MAPUNIT")) |>
  left_join(efg_info, by = "EFG") |>
  dplyr::select(
    MAPUNIT,
    ENV_DESC,
    eligible,
    EFG,
    efg_name,
    efg_pixel_value,
    efg_hex_code,
    efg_rgb_code,
    round,
    prop_memb,
    total_memb,
    supporting,
    total_interp
  )

write_csv(xref_best_tbl, here::here(out_dir, xref_best_file))

## 5b. Export tabular lookup  tables

efg_lookup_tbl <- xref_best_tbl |>
  dplyr::filter(!is.na(EFG)) |>
  summarise(
    n_ntveg_units = n_distinct(`MAPUNIT`),
    .by = c(efg_pixel_value, EFG, efg_name, efg_hex_code, efg_rgb_code)
  ) |>
  left_join(
    as_tibble(IUCN_GET) |>
      dplyr::select(efg_id, realm_id, realm_name, biome_id, biome_name),
    by = c("EFG" = "efg_id")
  ) |>
  dplyr::select(
    efg_pixel_value,
    EFG,
    efg_name,
    realm_id,
    realm_name,
    biome_id,
    biome_name,
    efg_hex_code,
    efg_rgb_code,
    n_ntveg_units
  ) |>
  arrange(efg_pixel_value)

write_csv(
  efg_lookup_tbl,
  here::here(out_dir, efg_lookup_file)
)

## 5c. Export colour file

clr_lines <- efg_lookup_tbl |>
  tidyr::separate(efg_rgb_code, into = c("R", "G", "B"), sep = ",") |>
  transmute(
    line = sprintf(
      "%d %s %s %s 255 %s %s",
      efg_pixel_value,
      R,
      G,
      B,
      EFG,
      efg_name
    )
  ) |>
  pull(line)

writeLines(
  clr_lines,
  here::here(out_dir, xref_clr_file)
)

## 5d. Export style file

qml_entries <- efg_lookup_tbl |>
  transmute(
    entry = sprintf(
      '        <paletteEntry value="%d" color="%s" alpha="255" label="%s %s"/>',
      efg_pixel_value,
      efg_hex_code,
      EFG,
      efg_name
    )
  ) |>
  pull(entry)

qml_lines <- c(
  "<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>",
  '<qgis version="3.34.0" minScale="1e+08" maxScale="0" hasScaleBasedVisibilityFlag="0" styleCategories="AllStyleCategories">',
  "  <pipe>",
  "    <provider>",
  '      <resampling enabled="false" maxOversampling="2" zoomedInResamplingMethod="nearestNeighbour" zoomedOutResamplingMethod="nearestNeighbour"/>',
  "    </provider>",
  '    <rasterrenderer band="1" type="paletted" opacity="1" alphaBand="-1" nodataColor="">',
  "      <rasterTransparency/>",
  "      <minMaxOrigin>",
  "        <limits>None</limits>",
  "        <extent>WholeRaster</extent>",
  "        <statAccuracy>Estimated</statAccuracy>",
  "        <cumulativeCutLower>0.02</cumulativeCutLower>",
  "        <cumulativeCutUpper>0.98</cumulativeCutUpper>",
  "        <stdDevFactor>2</stdDevFactor>",
  "      </minMaxOrigin>",
  "      <colorPalette>",
  qml_entries,
  "      </colorPalette>",
  "    </rasterrenderer>",
  '    <brightnesscontrast brightness="0" contrast="0" gamma="1"/>',
  '    <huesaturation colorizeGreen="128" colorizeOn="0" colorizeRed="255" colorizeBlue="128" grayscaleMode="0" saturation="0" colorizeStrength="100"/>',
  '    <rasterresampler maxOversampling="2"/>',
  "  </pipe>",
  "  <blendMode>0</blendMode>",
  "</qgis>"
)
writeLines(qml_lines, here::here(out_dir, xref_qml_file))

## 5f. Export raster files

## Check for non-matching MAPUNITs between the spatial data match and those in the membership table
ntveg_best_summary |>
  anti_join(ntveg_xy)

ntveg_best_summary |>
  anti_join(ntveg_xy)

#  Combined attribute table (vector data)
ntveg_xy_best <- ntveg_xy |>
  left_join(
    ntveg_best_summary,
    by = "MAPUNIT",
    relationship = "many-to-one"
  ) |>
  left_join(efg_pixel_values, by = "EFG")

# template raster file
### template raster file
rez <- 100
ntveg_raster <- raster(ntveg_xy, res = rez)

#### Best candidate

target_raster <- fasterize(
  ntveg_xy_best,
  ntveg_raster,
  field = "efg_pixel_value",
  fun = "first"
) |>
  rast()
## Supporting interpreters

supp_raster <- fasterize(
  ntveg_xy_best,
  ntveg_raster,
  field = "supporting",
  fun = "max"
) |>
  rast()

rast_metadata["TITLE"] <- "Nr. of interpreters supporting assignment"

metags(supp_raster) <- rast_metadata

writeRaster(supp_raster, here::here(out_dir, supp_tiff_file), overwrite = TRUE)
