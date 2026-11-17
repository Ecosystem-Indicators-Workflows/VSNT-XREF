# ---- Functions ----------------------------------------------------------------

#' Download and unzip the NT vegetation data
#'
#' @description This function downloads and unzips the NT vegetation data from a specified URL to a specified output file and directory.
#' @param data_url The URL of the data to download.
#' @param output_file The path to the output file where the downloaded data will be saved.
#' @param gisdata_path The path to the directory where the unzipped data will be saved.
#' @param ... Additional arguments to pass to the download.file function.
#' @param overwrite Logical. If TRUE, overwrite existing files. Default is FALSE.
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
#' @description This function pivots the membership tables from wide to long format.
#' @param xref The cross-reference data.
#' @param original_names The names of the original columns.
#' @param target_names The names of the target columns.
#' @param membership_values The names of the membership value columns.
#' @param expert_name The name of the expert/interpreter. Default is NULL.
#' @param round_nr The round number. Default is NULL.
#' @return A long-format membership table.
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

#' Clean membership table by filtering out records with invalid membership values, source IDs, or target IDs.
#' @author JRFP
#' @source Copied from https://github.com/Ecosystem-Indicators-Workflows/WAVSA-XREF/blob/1a9cb836d0120ac0f2276557804f07584377b6a9/scripts/functions.R#L33-L128
clean_membership <- function(x, source, target, action = "filter") {
  message("Using `clean_membership` from scripts/functions.R")
  stopifnot(inherits(x, "membership_tbl"))
  stopifnot(inherits(source, "lookup_tbl"))
  stopifnot(inherits(target, "lookup_tbl"))

  membership <- attr(x, "membership")
  source_id <- attr(x, "source_id")
  target_id <- attr(x, "target_id")
  round_nr <- attr(x, "round")
  interp_name <- attr(x, "interpreter")
  m <- x[[membership]]
  f1 <- is.na(m)
  f2 <- !f1 & (m < 0 | m > 1)

  sid <- x[[source_id]]
  vid <- rlang::sym(attr(source, "id"))
  elig <- rlang::sym(attr(source, "eligible"))
  eid <- dplyr::filter(source, !!elig) |> dplyr::pull(!!vid)
  f3 <- !(sid %in% eid)

  tid <- x[[target_id]]
  vid <- attr(target, "id")
  elig <- attr(target, "eligible")
  eid <- target[[vid]][target[[elig]]]
  f4 <- !(tid %in% eid)

  f5 <- duplicated(x)
  f6 <- duplicated(select(
    x,
    all_of(c(source_id, target_id, interp_name, round_nr))
  ))

  res <- x |> mutate(flag = "")

  if (any(f1)) {
    message(sprintf("Flagged %s records with `NA`", sum(f1)))
    res <- res |>
      mutate(flag = if_else(f1, "NA membership", flag))
  }
  if (any(f2)) {
    message(sprintf("Flagged %s records outside bounds [0,1]", sum(f2)))
    res <- res |>
      mutate(
        flag = if_else(
          f2,
          paste(flag, "membership out of bounds", sep = ";"),
          flag
        )
      )
  }
  if (any(f3)) {
    message(sprintf("Flagged %s records with not eligible source ids", sum(f3)))
    res <- res |>
      mutate(
        flag = if_else(f3, paste(flag, "not eligible source", sep = ";"), flag)
      )
  }
  if (any(f4)) {
    message(sprintf("Flagged %s records with invalid target ids", sum(f4)))
    res <- res |>
      mutate(flag = if_else(f4, paste(flag, "invalid target", sep = ";"), flag))
  }

  if (any(f5)) {
    message(sprintf("Flagged %s duplicated records ", sum(f5)))
    res <- res |>
      mutate(
        flag = if_else(f5, paste(flag, "duplicated record", sep = ";"), flag)
      )
  }
  if (any(f6)) {
    message(sprintf("Flagged %s duplicated assignments ", sum(f6)))
    res <- res |>
      mutate(
        flag = if_else(
          f6,
          paste(flag, "duplicated assignment", sep = ";"),
          flag
        )
      )
  }

  if (action == "filter") {
    res <- dplyr::filter(x, !f1 & !f2 & !f3 & !f4 & !f5 & !f6)
    if (any(f1 | f2 | f3 | f4 | f5 | f6)) {
      message(sprintf(
        "Removed %s records with one or more flags",
        sum(f1 | f2 | f3 | f4 | f5 | f6)
      ))
    } else {
      message("No records removed")
    }
  }
  res
}
