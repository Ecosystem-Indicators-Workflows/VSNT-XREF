# ---- Setup ----

# Libraries

## tidyverse libraries
library(dplyr)
library(readr)

## Read and write other formats
library(jsonlite)

## Spatial data
library(terra)

## paths
library(here)

out_dir <- here::here("results", "output1-best-candidate-raster")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ----- Metadata -----

# Extract metadata from codemeta
codemeta <- read_json(here::here("codemeta.json"))

# Update these:
metadata <- read_json(here::here("results/template-zenodo.json"))
datapkg <- read_json(here::here("results/template-datapackage.json"))

# Write here:

metadata_out <- here::here(out_dir, ".zenodo.json")
datapkg_out <- here::here(out_dir, "datapackage.json")

# Input files

input_filenames <- c(
  xref_best_file = "ntveg-xref-best-candidate.csv",
  efg_lookup_file = "IUCN-GET-pixel-lookup.csv",
  xref_tiff_file = "ntveg-xref-target.tiff",
  memb_tiff_file = "ntveg-xref-membership.tiff",
  supp_tiff_file = "ntveg-xref-supporting.tiff",
  xref_clr_file = "ntveg-xref-target.clr",
  xref_aux_file = "ntveg-xref-target.tiff.aux.xml",
  xref_qml_file = "ntveg-xref-target.qml"
)
input_files <- here::here(out_dir, input_filenames)
names(input_files) <- names(input_filenames)

input_files_info <- file_checks(
  input_files,
  here::here("data/metadata/data-checksums.txt"),
  "input"
)

xref_tiff_filename <- input_filenames["xref_tiff_file"]
efg_lookup_filename <- input_filenames["efg_lookup_file"]

xref_best_tbl <- read_csv(input_files["xref_best_file"])
efg_lookup_tbl <- read_csv(input_files["efg_lookup_file"])
r_target <- rast(input_files["xref_tiff_file"])


## ----- resources for datapackage -----

rs <- list()

the_round <- unique(xref_best_tbl$round) |> purrr::discard(is.na)
total_experts <- max(xref_best_tbl$total_interp, na.rm = TRUE)
min_memb <- min(xref_best_tbl$total_memb, na.rm = TRUE)
max_memb <- max(xref_best_tbl$total_memb, na.rm = TRUE)
min_supporting <- min(xref_best_tbl$supporting, na.rm = TRUE)
max_supporting <- max(xref_best_tbl$supporting, na.rm = TRUE)
min_prop <- min(xref_best_tbl$prop_memb, na.rm = TRUE)
max_prop <- max(xref_best_tbl$prop_memb, na.rm = TRUE)

how_many_units <- n_distinct(xref_best_tbl$MAPUNIT)
nr_assigned <- xref_best_tbl |>
  dplyr::filter(!is.na(EFG)) |>
  n_distinct("MAPUNIT")
total_rounds <- max(the_round, na.rm = TRUE)
nr_eligible <- xref_best_tbl |>
  dplyr::filter(eligible) |>
  n_distinct("MAPUNIT")
nr_targets <- n_distinct(xref_best_tbl$EFG, na.rm = TRUE)

# Spatial grid statistics
crs_info <- crs(r_target, describe = TRUE)
crs_label <- sprintf(
  "%s (%s:%s)",
  crs_info$name[1],
  crs_info$authority[1],
  crs_info$code[1]
)

res_xy <- res(r_target)

n_rows <- nrow(r_target)
n_cols <- ncol(r_target)
n_bands <- nlyr(r_target)
total_cells <- ncell(r_target) #TODO: check downstream usage, does not appear to be used
ext_bbox <- as.vector(ext(r_target)) # xmin, xmax, ymin, ymax

dimensions <- list(
  rows = n_rows,
  columns = n_cols,
  bands = n_bands
)

datapkg$spatial <-
  c(
    crs = paste(crs_info[2:3], collapse = ":"),
    crs_name = crs_info$name,
    bbox = list(ext_bbox),
    resolution_meters = list(res_xy),
    dimensions = dimensions
  )

### ----- 'Cross-reference summary' (`xref_best_tbl`) datapackage resource -----

descriptions <- c(
  MAPUNIT = "NT Vegetation survey unit identifier for the source units.",
  ENV_DESC = "Description of the source unit from the NT Vegetation survey.",
  eligible = "Whether the source unit is eligible for cross-referencing to an ecosystem concept in the IUCN GET.",
  EFG = "Best-candidate IUCN Global Ecosystem Typology Level 3 Ecosystem Functional Group code (e.g., T2.5). Empty when no interpreter cross-referenced the MAPUNITS.",
  efg_name = "Full name of the best-candidate IUCN GET Ecosystem Functional Group.",
  efg_pixel_value = glue::glue(
    "Standard integer pixel value used in {xref_tiff_filename} and the Global Ecosystems Atlas (geco package)."
  ),
  efg_hex_code = "Standard hexadecimal color code (#RRGGBB) for the EFG from the geco package.",
  efg_rgb_code = "Comma-separated RGB color values (R,G,B) for the EFG.",
  round = glue::glue(
    "Expert elicitation round number used for synthesis ({the_round} = revised assessment after group discussion)."
  ),
  prop_memb = glue::glue(
    "Proportion of total Round {the_round} membership values across all interpreters assigned to the best-candidate EFG (0 to 1)."
  ),
  total_memb = cli::pluralize(
    "Sum of membership values across all interpreters for the source unit in Round {the_round} (data range: {min_memb} to {max_memb} indicates errors, should be equal to {total_experts * 1.0} when {total_experts} interpreter{?s} assessed the unit with weight{?s} summing to 1)."
  ),
  supporting = cli::pluralize(
    "Number of interpreters in Round {the_round} who assigned a non-zero membership (> 0) to the best-candidate EFG (range: {min_supporting} to {max_supporting}).",
    "{qty(as.integer(min_supporting == max_supporting))} {?Not every} interpreter assessed every MAPUNITS."
  ),
  total_interp = pluralize(
    "Number of interpreters who assessed the source unit (1 to {total_experts}). {?Not every} interpreter assessed every MAPUNITS."
  )
)

finfo <- input_files_info |> slice(match("xref_best_file", shortname))

rs[[1]] <- list(
  name = "Crossreference summary",
  path = finfo$file_name,
  profile = "tabular-data-resource",
  format = "csv",
  mediatype = "text/csv",
  hash = unname(finfo$md5sum),
  encoding = guess_encoding(names(finfo$md5sum))$encoding[1],
  description = glue::glue(
    "Synthesised best-candidate cross-reference table linking all {how_many_units} MAPUNITS of Vegetation survey of the Northern Territory to IUCN GET Level 3 Ecosystem Functional Groups from Round {total_rounds} expert elicitation."
  ),
  schema = list(
    primaryKey = "MAPUNIT",
    fields = list_schema_fields_tbl(xref_best_tbl, descriptions)
  )
)

### ----- 'EFG Lookup Table' (`efg_lookup_tbl`) datapackage resource -----

min_pixel <- min(efg_lookup_tbl$efg_pixel_value)
max_pixel <- max(efg_lookup_tbl$efg_pixel_value)

descriptions <- c(
  efg_pixel_value = glue::glue(
    "Integer pixel value encoded in {xref_tiff_filename} ({min_pixel} to max_pixel)."
  ),
  EFG = "IUCN GET Level 3 Ecosystem Functional Group code (e.g., T2.3, TF1.6).",
  efg_name = "Full name of the Ecosystem Functional Group.",
  realm_id = "IUCN GET Level 1 Realm code (T, TF, F, MT, MFT).",
  realm_name = "IUCN GET Level 1 Realm name.",
  biome_id = "IUCN GET Level 2 Functional Biome code.",
  biome_name = "IUCN GET Level 2 Functional Biome name.",
  efg_hex_code = "Standard hexadecimal color code (#RRGGBB).",
  efg_rgb_code = "Comma-separated RGB color components (R,G,B).",
  n_ntveg_units = "Number of vegetation units (`MAPUNITS`) on the map assigned to this EFG as their best candidate."
)


finfo <- input_files_info |> slice(match("efg_lookup_file", shortname))

rs[[2]] <- list(
  name = "Pixel Lookup table",
  path = finfo$file_name,
  profile = "tabular-data-resource",
  format = "csv",
  mediatype = "text/csv",
  hash = unname(finfo$md5sum),
  encoding = guess_encoding(names(finfo$md5sum))$encoding[1],
  description = glue::glue(
    "Lookup table for the {nr_targets} IUCN Global Ecosystem Typology Level 3 Ecosystem Functional Groups mapped as best candidates in Northern Territory, including integer raster pixel values, realm/biome hierarchy, and color codes."
  ),
  schema = list(
    primaryKey = "efg_pixel_value",
    fields = list_schema_fields_tbl(efg_lookup_tbl, descriptions)
  )
)

### ----- 'Best candidate EFG raster' datapackage resource -----

finfo <- input_files_info |> slice(match("xref_tiff_file", shortname))

rs[[3]] <- list(
  name = "Best candidate EFG raster",
  path = finfo$file_name,
  format = "GeoTIFF",
  mediatype = "image/tiff; application=geotiff",
  hash = unname(finfo$md5sum),
  description = glue::glue(
    "90 m raster of the best-candidate IUCN GET Ecosystem Functional Group (integer pixel values {min_pixel} to {max_pixel} matching {efg_lookup_filename}, with embedded color table) in {crs_label}."
  )
)

### ----- 'Best candidate EFG sidecar file' datapackage resource -----

finfo <- input_files_info |> slice(match("xref_aux_file", shortname)) # TODO: consider alternative approach because if line fails, previous value of finfo will be used, which will silently generate incorrect metadata for the sidecar file. Consider using a separate variable name for each file's finfo or calling without variable name for the 'hash' variable.

rs[[4]] <- list(
  name = "Best candidate EFG sidecar file",
  path = finfo$file_name,
  format = "XML",
  mediatype = "application/xml",
  hash = unname(finfo$md5sum),
  description = glue::glue(
    "Persistent Auxiliary Metadata (PAM) with categorical information for {xref_tiff_filename}."
  )
)

### ----- 'Best candidate EFG colour file' datapackage resource -----

finfo <- input_files_info |> slice(match("xref_clr_file", shortname))

rs[[5]] <- list(
  name = "Best candidate EFG colour file",
  path = finfo$file_name,
  format = "CLR",
  mediatype = "text/plain",
  hash = unname(finfo$md5sum),
  description = glue::glue(
    "GDAL / ArcGIS RGBA color map file for {xref_tiff_filename}."
  )
)

### ----- 'Best candidate membership raster' datapackage resource -----

finfo <- input_files_info |> slice(match("memb_tiff_file", shortname))

rs[[7]] <- list(
  name = "Best candidate membership raster",
  path = finfo$file_name,
  format = "GeoTIFF",
  mediatype = "image/tiff; application=geotiff",
  hash = unname(finfo$md5sum),
  description = glue::glue(
    "90 m floating-point raster of the proportion of total membership votes (prop_memb, range {min_prop} to {max_prop}) supporting the best-candidate EFG in {crs_label}."
  )
)

### ----- 'Best candidate membership raster' datapackage resource -----

finfo <- input_files_info |> slice(match("supp_tiff_file", shortname))

rs[[8]] <- list(
  name = "Best candidate support raster",
  path = finfo$file_name,
  format = "GeoTIFF",
  mediatype = "image/tiff; application=geotiff",
  hash = unname(finfo$md5sum),
  description = glue::glue(
    "90 m raster of the number of interpreters (supporting, range {min_supporting} to {max_supporting} out of {total_experts}) assigning non-zero membership to the best-candidate EFG in {crs_label}."
  )
)

### ----- Best candidate -----

metadata$title <- datapkg$title <- "Cross-referencing vegetation units of the Vegetation Survey of the Northern Territory to the IUCN Global Ecosystem Typology: 90m Best-Candidate Spatial Synthesis"

datapkg$name <- "ntveg-xref-iucn-get-90m-best-candidate"

datapkg$description <- cli::pluralize(
  "Spatial synthesis (90 m GeoTIFF rasters) and tabular lookup tables cross-referencing the vegetation units of the Vegetation Survey of the Northern Territory to Level 3 Ecosystem Functional Groups (EFGs) of the IUCN Global Ecosystem Typology; based on Round {the_round} structured expert elicitation across {xfun::numbers_to_words(total_experts)}{qty(total_experts)} interpreter{?s}."
)
datapkg$contributors <- metadata$creators
datapkg$created <- metadata$publication_date <- format(
  Sys.Date(),
  format = "%Y-%m-%d"
)

metadata$description <- glue::glue(
  "<p>This dataset provides a 90 m spatial synthesis and tabular lookup cross-referencing the vegetation units (`MAPUNITS`) of the Vegetation Survey of the Northern Territory to Level 3 Ecosystem Functional Groups (EFGs) of the <strong>IUCN Global Ecosystem Typology (IUCN GET)</strong>.</p>

<p>Cross-referencing relationships were established by {total_experts} domain experts using fuzzy membership values; the latest assessment of each interpreter is used. This package provides the synthesised <strong>best-candidate ('winner-takes-all')</strong> spatial layers ({nr_assigned} MAPUNITS assigned to {nr_targets} IUCN GET Ecosystem Functional Groups, of which {nr_assigned_map} MAPUNITS are on the map), along with pixel-level membership proportions, number of supporting interpreters, tabular cross-reference tables, and GIS styling files. MAPUNITS that no interpreter cross-referenced are not mapped.</p>

<p>The EFG <strong>T5.6 Semi-arid woodlands</strong> is used by the interpreters but is not in the IUCN GET Level 3 list; it is mapped with pixel value 35.</p>

<p><strong>Attribution requirement:</strong> Spatial layers are derived from the Department of Environment and Natural Resources &copy; Northern Territory Government, available via the <a href=\"https://datasets.seed.nsw.gov.au/dataset/nsw-state-vegetation-type-map\">Northern Territory Land Information System</a> (NTLIS).</p>"
) # TODO replace link with NT dataset

metadata$related_identifiers <- c(
  metadata$related_identifiers,
  list(
    list(
      identifier = codemeta$codeRepository,
      relation = "isDerivedFrom",
      resource_type = "software",
      scheme = "url"
    ),
    list(
      identifier = "https://ftp-dlrm.nt.gov.au/main.html?download&weblink=e78bd012813d86c84c37681c72021704&realfilename=VegResources_ntveg_1M.zip",
      relation = "isDerivedFrom",
      resource_type = "dataset",
      scheme = "url"
    ),
    list(
      identifier = "urn:anzlic:dataset:ANZNT0001000081",
      relation = "isDerivedFrom",
      resource_type = "dataset",
      scheme = "urn"
    )
  )
)

source_note <- glue::glue(
  "Vegetation Survey of the Northern Territory ",
  "\u00A9 Department of Environment and Natural Resources; ",
)
metadata$notes <- sprintf(metadata$notes, source_note)
metadata$version <- datapkg$version <- codemeta$version
datapkg$homepage <- codemeta$codeRepository

datapkg$sources <- c(
  datapkg$sources,
  list(list(
    title = "Vegetation Survey of the Northern Territory (urn:anzlic:dataset:ANZNT0001000081)",
    path = "https://ftp-dlrm.nt.gov.au/main.html?download&weblink=e78bd012813d86c84c37681c72021704&realfilename=VegResources_ntveg_1M.zip",
    email = "Datarequests.DLPE@nt.gov.au"
  ))
)

datapkg$licenses <- c(
  datapkg$licenses,
  list(list(
    name = "CC-BY-4.0",
    title = "Creative Commons Attribution 4.0 International (Supplied by the Department of Environment and Natural Resources. Copyright Northern Territory Government)",
    path = "https://creativecommons.org/licenses/by/4.0/legalcode"
  ))
)

metadata$keywords <- c(unlist(metadata$keywords), "Northern Territory")

datapkg$resources <- rs

write_json(
  metadata,
  metadata_out,
  auto_unbox = TRUE,
  pretty = TRUE
)

write_json(
  datapkg,
  datapkg_out,
  auto_unbox = TRUE,
  pretty = TRUE
)

if (fs::file_exists(datapkg_out) & fs::file_exists(metadata_out)) {
  cli::cli_alert_success(
    "Output metadata written: {basename(metadata_out)}, {basename(datapkg_out)}"
  )
} else {
  cli::cli_alert_danger(
    "Output metadata not successfully written: {basename(metadata_out)}, {basename(datapkg_out)}"
  )
}
