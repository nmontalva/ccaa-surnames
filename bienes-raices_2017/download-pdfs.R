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

# basename <- function(url)
#   url %>% str_split("/") %>% .[[1]] %>% tail(1)
# 
# # pdfnames <- map(pdfs, basename)
# download <- function(url) {
#   bn <- basename(url)
#   # fn <- function(i) {
#   #   fname <- if (i == 1) bn else
#   #     bn %>% str_split("\\.") %>% .[[1]] %>%
#   #       c(str_c("-", i, ".")) %>% .[c(1,3,2)] %>%
#   #       str_c(collapse="")
#   #   if (file.exists(fname)) fn(i+1)
#   #   else fname
#   # }
#   # fname <- fn(1)
#   comuna <- url %>% str_split("/") %>% .[[1]] %>%
#     .[(length(.)-1)] # tail(2) %>% head(1)
#   fname <- str_c("pdfs/", comuna, "-", bn)
#   download.file(url, fname, mode="wb")
# }

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

main <- function() {
  pdfs <- getPDFUrls(nominasUrl)
  typos <- typoMap(typosFn)
  infoTbl <- infoTable(pdfs, typos)
  download(infoTbl)
  infoTbl %>% select(-last_year) %>%
    write_csv("info.csv")
}

