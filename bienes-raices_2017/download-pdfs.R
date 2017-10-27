library(rvest)
library(stringr)
library(tidyverse)

nominasUrl <- "http://www.comunidadesagricolas.cl/wp-content/uploads/nominas/"
getTable <- function(url) {
  read_html(url) %>% 
    html_nodes("table") %>% 
    .[[1]] %>%
    html_table()
}
# first and last items are empty strings,
# second item is "Parent directory"
getItems <- function(url)
  getTable(url)$Name %>% .[3:(length(.)-1)]
getChildUrls <- function(url)
  map_chr(getItems(url), ~str_c(url, .))
getPDFUrls <- function(nominasUrl) {
  yurls <- getChildUrls(nominasUrl)
  curls <- unlist(map(yurls, getChildUrls))
  unlist(map(curls, getChildUrls))
}

typosFn <- "typos.csv"
typoMap <- function(filename) {
  tbl <- read_csv(filename, col_types = "cc")
  typos <- tbl$right
  names(typos) <- tbl$wrong
  typos
}

urlInfo <- function(url, typos) {
  info <- tibble(year=NA,comuna=NA,name=NA,url=url)
  info[c("year","comuna","name")] =
    str_split(url, "/")[[1]] %>% tail(3)
  nameParts = info$name %>% str_split("\\.") %>% .[[1]]
  if (length(nameParts) == 2 && nameParts[[2]] == "pdf") {
    fst = nameParts[[1]]
    if (fst %in% names(typos))
      info$name = typos[[fst]]
    else
      info$name = fst
  } else
    stop(str_interp(
      c("file \"${info$name}\" contains more than one dot ",
        "or its extension is different than \"pdf\"")))
  info
}

infoTable <- function(pdfs, typos) {
  pdfs %>% map_dfr(~urlInfo(., typos)) %>% as.tibble %>%
    group_by(name, comuna) %>%
    summarise(years=str_c(year, collapse=","),
              last_year=last(year),
              url=last(url))
}

download <- function(infoTbl) {
  infoTbl %>% by(seq_len(nrow(infoTbl)), function(i) {
    fname <- str_c("pdfs/", i$comuna, "-", i$last_year,
                   "-", i$name, ".pdf")
    download.file(i$url, fname, mode="wb")
  })
}

# TODO: change camelCase naming to snake_case

download_pdfs <- function() {
  pdfs <- getPDFUrls(nominasUrl)
  typos <- typoMap(typosFn)
  infoTbl <- infoTable(pdfs, typos)
  dir.create("pdfs", showWarnings=FALSE)
  download(infoTbl)
  infoTbl %>% select(-last_year) %>%
    write_csv("info.csv")
}

