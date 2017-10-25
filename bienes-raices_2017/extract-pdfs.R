library(stringr)
library(stringi)
library(tidyverse)
library(tabulizer)

on <- function(f, x) { f(x); x }
view <- function(x) on(View, x)
# n <- 0
# printn <- function(x) {
#   n <- n + 1
#   cat(n, ": ")
#   print(x)
# }

# read.csv(..., fill=TRUE) doesn't work!
my_read_csv <- function(filename, ncol=5) {
  read_lines(filename) %>%
    map(function(l) {
      xs <- l %>%
        # FIXME: quotes should be handled properly
        str_replace_all("\"", "") %>%
        str_split(",") %>% .[[1]]
      c(xs, rep("", 5-length(xs)))
    }) %>%
    unlist %>%
    matrix(ncol=ncol, byrow=TRUE)
}

num_shares <- function(s, fname) {
  if (s == "UN DERECHO") 1
  else if (s == "UNO") 1
  else if (s == "1 DERECHO") 1
  else if (s == "UN DEREHO") 1
  else if (s == "UN DEDECHO") 1
  else if (s == "UN  DERECHO") 1
  else if (s == "DOS DERECHOS") 2
  else if (s == "DOS DERECHO") 2
  else if (s == "TRES DERECHOS") 3
  else if (s == "TRES DERECHO") 3
  else if (s == "CUATRO DERECHOS") 4
  else if (s == "CINCO DERECHOS") 5
  else if (s == "1/2 DERECHO") 1/2
  else if (s == "1/3 DERECHO") 1/3
  else if (s == "1/4 DERECHO") 1/4
  else if (s == "3/4 DERECHO") 3/4
  else if (s == "1/8 DERECHO") 1/8
  else if (s == "1 1/2 DERECHO") 1+(1/2)
  else if (s == "1 1/2 DERECHOS") 1+(1/2)
  else if (s == "1 1/2  DERECHO") 1+(1/2)
  else if (s == "1 1/2  DERECHOS") 1+(1/2)
  else if (s == "1/1/2 DERECHO") 1+(1/2)
  else if (s == "1/1/2 DERECHOS") 1+(1/2)
  else if (s == "1 /1/2 DERECHOS") 1+(1/2)
  else if (s == "1 1/4 DERECHOS") 1+(1/4)
  else if (s == "1/1/4 DERECHOS") 1+(1/4)
  else if (s == "2 1/2 DERECHOS") 2+(1/2)
  else if (s == "2/1/2 DERECHOS") 2+(1/2)
  else if (s == "") 0
  else {
    warning(str_interp(
      "\"${s}\" isn't recognised in file \"${fname}\""))
    -1
  }
}

# is.complete <- function(df, fname) {
#   check_ids <- function(col, i, j, maybe_bigger) {
#     if (i >= length(col)) TRUE
#     else {
#       jc <- as.character(j)
#       if (col[[i]] == jc)
#         check_ids(col, i+1, j+1, FALSE)
#       else if (startsWith(col[[i]], jc))
#         check_ids(col, i+1, j, TRUE)
#       else if (maybe_bigger) {
#         jc <- as.character(j+1)
#         if (col[[i]] == jc)
#           check_ids(col, i+1, j+2, FALSE)
#         else if (startsWith(col[[i]], jc))
#           check_ids(col, i+1, j+1, TRUE)
#         else {
#           warning(str_interp(
#             c("expecting \"${jc}\" at row ${i} ",
#               "but found \"${col[[i]]}\" ",
#               "in file \"${fname}\"")))
#           FALSE
#         }
#       } else {
#         warning(str_interp(
#           c("expecting \"${jc}\" at row ${i} ",
#             "but found \"${col[[i]]}\" ",
#             "in file \"${fname}\"")))
#         FALSE
#       }
#     }
#   }
#   check_ids(df$right_id, 1, 1, FALSE)
# }

is.complete <- function(df, fname, flexible=TRUE) {
  check_ids <- function(col, i, j) {
    if (i >= length(col)) TRUE
    else if (col[[i]] == as.character(j))
      check_ids(col, i+1, j)
    else if (col[[i]] == as.character(j+1))
      check_ids(col, i+1, j+1)
    else if (flexible && col[[i]] == as.character(j+2))
      check_ids(col, i+1, j+2)
    else FALSE
    #   warning(str_interp(
    #     c("expecting \"${j}\" or \"${j+1}\" ",
    #       "at row ${i} but found \"${col[[i]]}\" ",
    #       "in file \"${fname}\"")))
    #   FALSE
    # }
  }
  length(df) == 5 && !is_empty(df) &&
    if (str_detect(fname, "limari-2010-oruro"))
      check_ids(df$right_id, 1, 7)
    else if (str_detect(fname, "limari-2010-alcones"))
      check_ids(df$right_id, 1, 4)
    else
      check_ids(df$right_id, 1, 1)
}

communityTable <- function(m, fname) {
  if (is_empty(m) || ncol(m) != 5)
    return(NULL)
  tibble(comunidad=m[,1],
         right_id= # the identifier of the commoner right
           str_replace(m[,2], "^0*([0-9]+).*", "\\1"),
         firstname=m[,3],
         surname=m[,4],
         shares=
           map_dbl(m[,5], ~num_shares(., fname)))
}

tidy <- function(df) {
  if (!is.data.frame(df) || is_empty(df))
    return(NULL)
  df %>%
    filter(right_id != "") %>% # firstname != "ELIMINADO") %>%
    mutate(right_id=as.integer(right_id)) %>%
    arrange(right_id)
}

formatTable <- function(m, filename) {
  if (str_detect(filename, "totoralillo")) {
    m <- m[,-c(1,2)]
    m[,1] <- rep("TOTORALILLO", nrow(m))
  } else if (str_detect(filename, "el_potrero_alto"))
    m[,1] <- rep("EL POTRERO ALTO", nrow(m))
  else if (str_detect(filename, "los_clonquis"))
    m <- cbind(rep("LOS CLONQUIS", nrow(m)-1), m[-1,1:4])
  else if (str_detect(filename, "los_morales"))
    m <- cbind(rep("LOS MORALES", nrow(m)-1), m[-1,1:4])
  # else if (str_detect(filename, "peral_ojo_de_agua"))
  #   m <- cbind(rep("PERAL OJO DE AGUA", nrow(m)-1), m[-1,])
  else if (str_detect(filename, "diaz_y_ocaranza"))
    m[,1] <- rep("DIAZ Y OCARANZA", nrow(m))
    # m <- cbind(rep("DIAZ Y OCARANZA", nrow(m)), m)
  else if (str_detect(filename, "atelcura|carrizo"))
    m <- cbind(m[,1:4], rep("UN DERECHO", nrow(m)))
  for (i in c(1,2,3)) {
    if (nrow(m) >= i && endsWith(m[i,2], "NOMBRES")) {
      # x <- m[,2]
      # ids <- str_replace(x, "^0*([0-9]*).*$", "\\1")
      # nms <- str_replace(x, "^([0-9]*)[^ ]* ", "")
      # if (ncol(m) == 5) {
      #   m[,2] <- ids
      #   m[,3] <- nms
      # } else if (ncol(m) == 4) {
      #   m <- cbind(m[,1], ids, nms, m[,3:4])
      # } else
      #   warning("table has ", ncol(m), " columns")
      m <- m[-(1:i),]
      # return(m)
    }
  }
  for (i in c(1,2,3)) {
    if (nrow(m) >= i && startsWith(m[i,3], "NOMBRES")) {
      m <- m[-(1:i),]
      # return(m)
    }
  }
  if (ncol(m) == 5 &&
      any(str_detect(m[,2], "^[0-9].* [A-Z][A-Z]+"))) {
    x <- m[,2]
    w <- which(str_detect(x, "^[0-9].* [A-Z][A-Z]+"))
    ids <- str_replace(x[w], "^0*([0-9]*).*$", "\\1")
    nms <- str_replace(x[w], "^([0-9]*)[^ ]* ", "")
    for (i in 4:5) {
      if (all(m[w,i] == ""))
        for (j in i:4) {
          m[w,j] <- m[w,j-1]
        }
    }
    # m[w,5] <- m[w,4]
    # m[w,4] <- m[w,3]
    m[w,3] <- nms
    m[w,2] <- ids
  }
  if (ncol(m) == 4 &&
      all(str_detect(m[,2], "^[0-9].* [A-Z][A-Z]+"))) {
    x <- m[,2]
    ids <- str_replace(x, "^0*([0-9]*).*$", "\\1")
    nms <- str_replace(x, "^([0-9]*)[^ ]* ", "")
    m <- cbind(m[,1], ids, nms, m[,3:4])
  }
  if (any(m[,2] == "")) {
    x <- m[,3]
    ids <- str_replace(x, "^0*([0-9]*).*$", "\\1")
    nms <- str_replace(x, "^([0-9]*)[^ ]* ", "")
    w <- which(m[,2] == "")
    if (all(ids[w] != "")) {
      m[,2][w] <- ids[w]
      m[,3][w] <- nms[w]
    }
  }
  if (ncol(m) >= 4 && all(is.na(m[,4])))
    m[,4] <- rep("", nrow(m))
  if (ncol(m) >= 4 && any(endsWith(m[,4], "DERECHO"))) {
    x <- m[,4]
    t <- str_match(x, "^(.*) ([^ ]* DERECHO)$")
    w <- which(!is.na(t[,1]))
    m[,4][w] <- t[,2][w]
    if (ncol(m) == 4) {
      m <- cbind(m, t[,3])
      m[,5][which(is.na(m[,5]))] <- ""
    } else
      m[,5][w] <- t[,3][w]
  }
  # if (m[1,3] == "NOMBRES")
  #   m <- m[-1,]
  # else if (nrow(m) >= 2 && m[2,3] == "NOMBRES")
  #   m <- m[-c(1,2),]
  # else if (nrow(m) >= 3 && m[3,3] == "NOMBRES")
  #   m <- m[-c(1,2,3),]
  # if (m[1,5] == "DERECHOS")
  #   m <- m[-1,]
  # else if (nrow(m) >= 2 && m[2,5] == "DERECHOS")
  #   m <- m[-c(1,2),]
  # else if (nrow(m) >= 3 && m[3,5] == "DERECHOS")
  #   m <- m[-c(1,2,3),]
  m
}

read_pdf <- function(filename, method=c("tabulizer",
                                        "tabula",
                                        "pdftohtml",
                                        "pdftotext")) {
  if (length(method) == 0) {
    warning("no methods left to try on \"", filename, "\"")
    data.frame(stringsAsFactors=FALSE)
  } else {
    m <- method[[1]]
    df <- if (m == "tabula") {
      read_pdf_tabula(filename)
    } else if (m == "tabulizer") {
      read_pdf_tabulizer(filename)
    } else if (m == "pdftohtml") {
      read_pdf_pdftohtml(filename)
    } else if (m == "pdftotext") {
      read_pdf_pdftotext(filename)
    } else stop("method \"${method}\" not recognised")
    if (is.data.frame(df) && is.complete(df, filename)) df
    else read_pdf(filename, method[-1])
  }
}

read_pdf_tabulizer <- function(filename) {
  if (str_detect(filename, "potrerillo_alto"))
    return(NULL)
  tbls <- tryCatch(
    extract_tables(filename),
    error=function(e) {
      warning("error executing tabulizer::extract_tables(\"",
              filename, "\"): ", e)
    })
  if (is.character(tbls)) tbls # tbls contains the warning
  else if (length(tbls) > 0) {
    tbls %>% map_dfr(function(m) {
      m %>% formatTable(filename) %>%
        communityTable(filename)
    }) %>% tidy
  } else NULL
    # warning(str_interp("\"${filename}\" has a weird format"))
}

read_pdf_tabula <- function(filename) {
  basename <- str_split(filename, "\\.")[[1]][[1]]
  if (!file.exists(str_c(basename, ".csv")))
    system(str_interp(
      c("java -jar tabula-1.0.1-jar-with-dependencies.jar ",
        "${filename} -o ${basename}.csv -p all")))
  # m <- read.csv(str_c(basename, ".csv"),
  #               header=FALSE,
  #               stringsAsFactors=FALSE)
  # m <- read_csv(str_c(basename, ".csv"),
  #               col_names=c("X1","X2","X3","X4","X5"),
  #               col_types="ccccc", na=character()) %>% as.matrix
  m <- my_read_csv(str_c(basename, ".csv"))
  m %>% formatTable(filename) %>%
    communityTable(filename) %>% tidy
  # if (nrow(m) == 0)
  #   warning("\"", filename, "\" is empty")
  # else if (ncol(m) != 5) NULL
  #   # warning("\"", filename, "\" contains ",
  #   #         ncol(m), " instead of 5 columns")
  # else {
  #   m %>% formatTable %>% communityTable(filename) %>% tidy
  # }
}

read_pdf_pdftohtml <- function(filename) {
  # system("pdftohtml")
}

read_pdf_pdftotext <- function(filename) {
  # system("pdftotext")
}

skip <- c(
  "elqui-2010-gualliguaica",
  "limari-2010-carcamo",
  "limari-2010-castillo_mal_paso_y_otros",
  "limari-2010-el_espinal",
  "limari-2010-el_potrero_de_huatulame",
  "limari-2010-fundina_norte",
  "limari-2010-rinconada_punitaqui",
  "limari-2010-cerro_blanco"
)

manually_extracted_csvs <-
  str_c("manually-extracted/", skip, "-M.csv")

col_names <- c("comunidad", "right_id",
               "firstname", "surname", "shares")

read_community_csv <- function(filename) {
  read_csv(filename, col_names = col_names,
           col_types = "ciccd")
}

remove_accents <- function(df) {
  for (col in c("comunidad", "firstname", "surname")) {
    if (any(str_detect(df[[col]], "_")))
      stop("\"_\" found in column \"", col,"\"")
    df[[col]] <- df[[col]] %>%
      str_replace_all("Ñ", "_") %>%
      stri_trans_general("latin-ascii") %>%
      str_replace_all("_", "Ñ")
  }
  df
}

remove_nonpeople <- function(df) {
  df %>% filter(firstname != "ELIMINADO")
}

fix_shares <- function(df) {
  # fix <- function(df, i, j) {
  #   if (j < i) stop(j, " < ", i, "!")
  #   if (sum(df$shares[i:j]) != df$shares[[i]])
  #     stop("sum(df$shares[", i, ":", j, "]) == ",
  #          sum(df$shares[i:j]), " != ", df$shares[[i]],
  #          " == df$shares[[", i, "]]")
  #   df$shares[i:j] <- df$shares[[i]] / (j-i+1)
  #   df
  # }
  # by_row <- function(df, i, last=NA, from=NA) {
  #   if (i >= nrow(df)) return(df)
  #   else if (is.na(last) ||
  #            df$comunidad[[last]] != df$comunidad[[i]]) {
  #     if (df$shares[[i]] == 0)
  #       stop("community \"", df$comunidad[[i]],
  #            "\" starts with commoner with zero shares")
  #     if (!is.na(from))
  #       df <- fix(df, from, i-1)
  #     from <- NA
  #   } else if (is.na(from) && df$shares[[i]] == 0) {
  #     from <- last
  #   } else if (!is.na(from) && df$shares[[i]] > 0) {
  #     df <- fix(df, from, i-1)
  #     from <- NA
  #   }
  #   by_row(df, i+1, i, from)
  # }
  # by_row(df, 1)
  check <- function(i, j) {
    if (j < i) stop(j, " < ", i, "!")
    if (sum(df$shares[i:j]) != df$shares[[i]])
      stop("sum(df$shares[", i, ":", j, "]) == ",
           sum(df$shares[i:j]), " != ", df$shares[[i]],
           " == df$shares[[", i, "]]")
  }
  last <- NA
  from <- NA
  for (i in seq_along(df$shares)) {
    if (is.na(last) ||
        df$comunidad[[last]] != df$comunidad[[i]]) {
      if (df$shares[[i]] == 0)
        warning("community \"", df$comunidad[[i]],
                "\" starts with commoner with zero shares ",
                "at row ", i)
      if (!is.na(from)) {
        check(from, i-1)
        df$shares[from:i-1] <- df$shares[[from]] / (i-from)
        # df <- fix(df, from, i-1)
      }
      from <- NA
    } else if (is.na(from) && df$shares[[i]] == 0) {
      from <- last
    } else if (!is.na(from) && df$shares[[i]] > 0) {
      check(from, i-1)
      df$shares[from:i-1] <- df$shares[[from]] / (i-from)
      # df <- fix(df, from, i-1)
      from <- NA
    }
    last <- i
  }
  df
}

fix_repeated <- function(df) {
  df %>%
    group_by(comunidad, firstname, surname) %>%
    summarise(right_id = first(right_id),
              shares = sum(shares)) %>%
              # count = n())
    .[,col_names] # reorder columns
}

# fix_surnames <- function(df) {
#   w <- which(df$surname == "")
# }

words_in_col <- function(df, col, ind) {
  df[[col]][ind] %>%
    str_split("[ ]+") %>%
    unlist %>%
    unique %>%
    sort
}

firstname_set <- function(df) {
  words_in_col(df, "firstname", which(df$surname != ""))
}

surname_set <- function(df) {
  words_in_col(df, "surname", which(df$surname != ""))
}

check_table <- function(df) {
  all(df$comunidad != "") &&
    all(df$right_id > 0) &&
    all(df$firstname != "") &&
    all(df$shares > 0)
}

comuneros <- function() {
  files <- list.files(path="pdfs", pattern="\\.pdf$",
                      full.names=TRUE)
  # skip files in skip
  skipr <- str_c(skip, collapse="|")
  files <- files[!str_detect(files, skipr)]
  rbind(
    files[1] %>% map_dfr(~read_pdf(print(.))),
    manually_extracted_csvs %>%
      map_dfr(~read_community_csv(print(.))) %>%
      mutate_if(is.character, ~str_replace_na(., "")),
    stringsAsFactors=FALSE)
}

main <- function() {
  comuneros() %>%
    remove_accents %>%
    remove_nonpeople %>%
    fix_shares %>%
    fix_repeated %>%
    arrange(comunidad, right_id) %>%
    write_csv("comuneros.csv")
}
