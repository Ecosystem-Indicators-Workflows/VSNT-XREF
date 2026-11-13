library(ecoSynthesis)
library(rlang)
library(sf)
library(readr)
library(dplyr)
library(fasterize)
library(raster)
library(terra)

# ---- Functions ----------------------------------------------------------------

# Create a function to download and unzip the NT vegetation data

get_data <- function(data_url, output_file, gisdata, ..., overwrite = FALSE) {
  rlang::check_required(data_url, output_file, gisdata)

  defaults <- list(
    mode = "wb",
      method = "wget",
      extra = c("-c")
  )
  
  user_args <- dots_list(..., .homonyms = "error")
  
  final_args <- modifyList(defaults, user_args, keep.null = TRUE)

  if (file.exists(output_file) && !overwrite) {
    cli::cli_abort(sprintf("%s already exists, skipping download\n", output_file))
  } else if (file.exists(output_file) && overwrite) {
    cli::cli_inform(sprintf("%s already exists, overwriting\n", output_file))
  } 
  cli::cli_alert("Downloading {data_url} to {output_file}") 
    rlang::inject(
      download.file(
      data_url,
      destfile = output_file,
      !!!final_args
    )
  )
  cli::cli_alert("Unzipping {output_file} to {gisdata}")
    rlang::inject(
      unzip(output_file, exdir = gisdata)
    )
}

# Create a function to pivot the membership tables from wide to long format 

pivot_wide_membership_table <- function(
  xref,
  original_names,
  target_names,
  membership_values,
  expert_name = NULL,
  round_nr = NULL
) {
  memb_table <- xref |>
    dplyr::select(
      dplyr::any_of(original_names),
      dplyr::any_of(target_names),
      dplyr::any_of(membership_values)
    ) |>
    tidyr::pivot_longer(
      dplyr::any_of(c(target_names, membership_values)),
      cols_vary = "slowest",
      names_to = c(".value"),
      names_pattern = "([[:alpha:]]+)"
    )

  if (!is.null(expert_name)) {
    memb_table <- mutate(memb_table, expert = expert_name)
  }
  if (!is.null(round_nr)) {
    memb_table <- mutate(memb_table, round = round_nr)
  }
  return(memb_table)
}

ntveg_round1 <- pivot_wide_membership_table(
  xref1,
  "MAPUNIT",
  sprintf("EFG_%s", 1:3),
  sprintf("memb_%s", 1:3),
  expert_name = "MGT",
  round_nr = 1
) |>
  membership_tbl("MAPUNIT", "EFG", "memb", "expert", "round")


ntveg_round2 <- pivot_wide_membership_table(
  xref2,
  "MAPUNIT",
  sprintf("EFG_%s", 1:3),
  sprintf("memb_%s", 1:3),
  expert_name = "MGT",
  round_nr = 2
) |>
  membership_tbl("MAPUNIT", "EFG", "memb", "expert", "round")

ntveg_memb <-
  bind_rows(ntveg_round1, ntveg_round2) |>
  arrange(MAPUNIT) |>
  filter(!is.na(EFG))

IUCN_GET <- lookup_tbl(IUCN_GET_LUT, "efg_id", "efg_name")

check_membership(ntveg_memb, source = ntveg_LUT, target = IUCN_GET)

unique(ntveg_memb$EFG)


rez <- 100
ntveg_xy <- st_transform(NT_ntveg, "EPSG:3112") |>
  filter(!is.na(MAPUNIT)) |>
  mutate(MAPUNIT = as.numeric(MAPUNIT))
ntveg_raster <- raster(ntveg_xy, res = rez)

ntveg_slc <- filter(ntveg_memb, round == 2)

for (tid in unique(ntveg_slc$EFG)) {
  cat(sprintf("Starting %s: %s // ", tid, date()))
  out_file <- here::here(out_dir, sprintf("%s-res%s.tiff", tid, rez))
  if (file.exists(out_file)) {
    cat(sprintf("%s already exists\n", out_file))
  } else {
    if (tid %in% ntveg_slc$EFG) {
      sum_membs <- ntveg_slc |>
        group_by(MAPUNIT) |>
        summarise(
          k = sum(if_else(EFG == tid, memb, 0)),
          total = sum(memb),
          prop = k / total
        )

      ntveg_k <- ntveg_xy |> left_join(sum_membs, by = "MAPUNIT")
      k_raster <- fasterize(ntveg_k, ntveg_raster, field = "prop", fun = "sum")

      print(sum(k_raster@data@values, na.rm = T))
      plot(k_raster)

      writeRaster(k_raster, out_file)
    } else {
      cat(sprintf("%s not found\n", tid))
    }
  }
}


nt_all <- rast(dir(out_dir, pattern = "res100", full.names = TRUE))

nt_which <- app(nt_all, fun = "which.max")
nt_max <- app(nt_all, fun = "max")
plot(nt_which)
plot(nt_max)
