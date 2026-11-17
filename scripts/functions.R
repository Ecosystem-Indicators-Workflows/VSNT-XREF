# ---- Functions ----------------------------------------------------------------

# Create a function to download and unzip the NT vegetation data

get_data <- function(
  data_url,
  output_file,
  gisdata_path,
  ...,
  overwrite = FALSE
) {
  rlang::check_required(data_url, output_file, gisdata_path)

  defaults <- list(
    mode = "wb",
    method = "wget",
    extra = c("-c")
  )

  user_args <- rlang::dots_list(..., .homonyms = "error")

  final_args <- utils::modifyList(defaults, user_args, keep.null = TRUE)

  if (file.exists(output_file) && !overwrite) {
    cli::cli_abort(sprintf(
      "%s already exists, skipping download\n",
      output_file
    ))
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
  cli::cli_alert("Unzipping {output_file} to {gisdata_path}")
  rlang::inject(
    unzip(output_file, exdir = gisdata_path)
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
    memb_table <- dplyr::mutate(memb_table, expert = expert_name)
  }
  if (!is.null(round_nr)) {
    memb_table <- dplyr::mutate(memb_table, round = round_nr)
  }
  return(memb_table)
}
