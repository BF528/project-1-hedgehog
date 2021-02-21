library(tibble)
library(magrittr)
library(purrr)
library(stringr)
library(tidyr)



read_gmt <- function(file, tidy = FALSE) {
  con <- file(file, "r")
  gmt_lines <- readLines(file, warn = FALSE)
  close(con)
  rlist <- purrr::map(gmt_lines, parse_gmt_lines)
  rlist_names <- purrr::map_chr(gmt_lines, get_gmt_names)
  names(rlist) <- rlist_names
  if (tidy) rlist <- tidy_gmt(rlist)
  return(rlist)
}


