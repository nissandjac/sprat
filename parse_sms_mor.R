parse_sms_mor <- function(path, species_map = NULL) {
  raw_lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  # sanitize to UTF-8 text; drop invalid bytes
  lines <- iconv(raw_lines, from = "", to = "UTF-8", sub = "")
  
  cur_species_id   <- NA_integer_
  cur_species_name <- NA_character_
  cur_year         <- NA_integer_
  ages             <- NULL
  season_idx       <- 0L
  out <- list()
  
  split_nums <- function(x) {
    # take only numeric tokens; ignore stray text
    toks <- unlist(strsplit(trimws(x), "\\s+"))
    suppressWarnings(as.numeric(toks))
  }
  
  clean_fence <- function(s) {
    s <- gsub("#+$", "", s)         # strip trailing #####
    s <- gsub("\\s+#+$", "", s)     # strip '  #####'
    trimws(s)
  }
  
  for (ln in lines) {
    tl <- trimws(ln)
    if (tl == "") next
    
    if (startsWith(tl, "#")) {
      # header lines
      if (grepl("^#\\s*Species:", tl, ignore.case = TRUE)) {
        lab <- sub(".*Species:\\s*", "", tl, ignore.case = TRUE)
        lab <- clean_fence(lab)
        
        # Try to read a leading numeric id
        id_try <- suppressWarnings(as.integer(sub("^\\s*(\\d+).*$", "\\1", lab)))
        if (!is.na(id_try)) {
          cur_species_id <- id_try
          # Try to extract a name after the id (e.g., "16  Sandeel")
          name_part <- sub("^\\s*\\d+\\s*[-:]?\\s*", "", lab)
          name_part <- trimws(name_part)
          if (nzchar(name_part)) {
            cur_species_name <- name_part
          } else if (!is.null(species_map) && as.character(cur_species_id) %in% names(species_map)) {
            cur_species_name <- species_map[[as.character(cur_species_id)]]
          } else {
            cur_species_name <- as.character(cur_species_id)
          }
        } else {
          # No leading id; treat whole thing as name
          cur_species_id <- NA_integer_
          cur_species_name <- lab
        }
        
      } else if (grepl("^#\\s*Year:", tl, ignore.case = TRUE)) {
        cur_year   <- as.integer(sub(".*Year:\\s*", "", tl, ignore.case = TRUE))
        season_idx <- 0L
        
      } else if (grepl("^#\\s*age\\b", tl, ignore.case = TRUE)) {
        age_part <- sub("^#\\s*age\\s*", "", tl, ignore.case = TRUE)
        ages <- split_nums(age_part)
      }
      next
    }
    
    # data rows
    if (is.null(ages)) next  # no age header seen yet
    
    vals <- split_nums(tl)
    if (length(vals) == 0) next
    
    # pad/truncate to age length
    if (length(vals) < length(ages)) {
      vals <- c(vals, rep(NA_real_, length(ages) - length(vals)))
    } else if (length(vals) > length(ages)) {
      vals <- vals[seq_along(ages)]
    }
    
    season_idx <- season_idx + 1L
    
    out[[length(out) + 1L]] <- data.frame(
      species_id   = cur_species_id,
      species_name = cur_species_name,
      year         = cur_year,
      season       = season_idx,
      age          = ages,
      value        = vals,
      stringsAsFactors = FALSE
    )
  }
  
  if (!length(out)) {
    return(data.frame(
      species_id=integer(), species_name=character(),
      year=integer(), season=integer(), age=numeric(), value=numeric()
    ))
  }
  do.call(rbind, out)
}
