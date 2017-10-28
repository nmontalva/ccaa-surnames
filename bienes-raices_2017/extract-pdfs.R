library(stringr)
library(stringi)
library(tidyverse)
library(tabulizer)

on <- function(f, .x) { f(.x); .x }
view <- function(x) on(View, x)

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

is_complete <- function(df, fname, flexible=TRUE) {
  check_ids <- function(col, i, j) {
    if (i >= length(col)) TRUE
    else if (col[[i]] == as.character(j))
      check_ids(col, i+1, j)
    else if (col[[i]] == as.character(j+1))
      check_ids(col, i+1, j+1)
    else if (flexible && col[[i]] == as.character(j+2))
      check_ids(col, i+1, j+2)
    else FALSE
  }
  length(df) == 5 && !is_empty(df) &&
    if (str_detect(fname, "limari-2010-oruro"))
      check_ids(df$right_id, 1, 7)
    else if (str_detect(fname, "limari-2010-alcones"))
      check_ids(df$right_id, 1, 4)
    else
      check_ids(df$right_id, 1, 1)
}

community_table <- function(m, fname) {
  if (is_empty(m) || ncol(m) != 5)
    return(NULL)
  tibble(community=m[,1],
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
    filter(right_id != "") %>%
    mutate(right_id=as.integer(right_id)) %>%
    arrange(right_id)
}

format_table <- function(m, filename) {
  if (str_detect(filename, "totoralillo")) {
    m <- m[,-c(1,2)]
    m[,1] <- rep("TOTORALILLO", nrow(m))
  } else if (str_detect(filename, "el_potrero_alto"))
    m[,1] <- rep("EL POTRERO ALTO", nrow(m))
  else if (str_detect(filename, "los_clonquis"))
    m <- cbind(rep("LOS CLONQUIS", nrow(m)-1), m[-1,1:4])
  else if (str_detect(filename, "los_morales"))
    m <- cbind(rep("LOS MORALES", nrow(m)-1), m[-1,1:4])
  else if (str_detect(filename, "diaz_y_ocaranza"))
    m[,1] <- rep("DIAZ Y OCARANZA", nrow(m))
  else if (str_detect(filename, "atelcura|carrizo"))
    m <- cbind(m[,1:4], rep("UN DERECHO", nrow(m)))
  for (i in c(1,2,3)) {
    if (nrow(m) >= i && endsWith(m[i,2], "NOMBRES"))
      m <- m[-(1:i),]
  }
  for (i in c(1,2,3)) {
    if (nrow(m) >= i && startsWith(m[i,3], "NOMBRES"))
      m <- m[-(1:i),]
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
  if (ncol(m) >= 4 && any(str_detect(m[,4], "DERECHOS?$"))) {
    t <- str_match(m[,4], "^(.* )?([^ ]* DERECHOS?)$")
    w <- which(!is.na(t[,1]))
    m[,4][w] <- str_replace_na(t[,2][w], "")
    if (ncol(m) == 4)
      m <- cbind(m, str_replace_na(t[,3], ""))
      # m[,5] <- str_replace_na(m[,5], "")
      # m[,5][which(is.na(m[,5]))] <- ""
    else
      m[,5][w] <- t[,3][w]
  }
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
    if (is.data.frame(df) && is_complete(df, filename)) df
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
  if (is.character(tbls)) tbls # when tbls contains the warning
  else if (length(tbls) > 0) {
    tbls %>% map_dfr(function(m) {
      m %>% format_table(filename) %>%
        community_table(filename)
    }) %>% tidy
  } else NULL
}

read_pdf_tabula <- function(filename) {
  basename <- str_split(filename, "\\.")[[1]][[1]]
  if (!file.exists(str_c(basename, ".csv")))
    system(str_interp(
      c("java -jar tabula-1.0.1-jar-with-dependencies.jar ",
        "${filename} -o ${basename}.csv -p all")))
  m <- my_read_csv(str_c(basename, ".csv"))
  m %>% format_table(filename) %>%
    community_table(filename) %>% tidy
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

col_names <- c("community", "right_id",
               "firstname", "surname", "shares")

read_community_csv <- function(filename) {
  read_csv(filename, col_names = col_names,
           col_types = "ciccd")
}

remove_accents <- function(df, cols=which(
  sapply(df, class) == 'character')) {
  for (col in cols) {
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
  # df %>% filter(firstname != "ELIMINADO",
  #               !startsWith(firstname, "COMUNIDAD AGRICOLA"),
  #               !str_detect(firstname, "CIA\.|LTDA\."))
}

fix_compound_names <- function(df) {
  for (col in c("firstname", "surname")) {
    df[[col]] <- df[[col]] %>%
      str_replace_all("DEL BUEN PASTOR", "DEL_BUEN_PASTOR") %>%
      str_replace_all("CASAS DEL VALLE", "CASAS_DEL_VALLE") %>%
      str_replace_all("(?<![A-Z])SAN +(MARTIN|LUIS|JUAN)",
                      "SAN_\\1") %>%
      str_replace_all("(?<![A-Z])DE +(LAS?|LOS) +", "DE_\\1_") %>%
      # negative lookbehind (?<![A-Z]) so that DEL is not
      # the end of a word (like in "FINDEL")
      str_replace_all("(?<![A-Z])(DEL?) +", "\\1_") %>%
      # str_replace_all("(?<![A-Z])(DEL?) ([A-Z]([A-Z]+|.))",
      #                 "\\1_\\2") %>%
      # replace twice to catch cases like "VDA DE DEL CAMPO"
      str_replace_all("(?<![A-Z])(DEL?) +", "\\1_") %>%
      # str_replace_all("(?<![A-Z])(DEL?) ([A-Z][A-Z]+)",
      #                 "\\1_\\2") %>%
      str_replace_all("(VIUDA|VDA.?) +", "VDA_")
  }
  df
}

remove_abbrev <- function(df) {
  for (col in c("firstname", "surname")) {
    df[[col]] <- df[[col]] %>%
      str_replace_all(" [A-Z]\\.", "") %>%
      # if the name starts by an initial
      str_replace_all("(?<![A-Z])[A-Z]\\.", "")
  }
  df
}

words_in_col <- function(df, col, ind=seq_along(df[[col]])) {
  df[[col]][ind] %>%
    str_split(" +") %>%
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

alt_regex <- function(options, complete_match=TRUE) {
  options %>%
    str_c(collapse="|") %>%
    str_c(if (complete_match) "^(", .,
          if (complete_match) ")$") %>%
    # for some weird reason the replacement string
    # is not "\\." but "\\\\."
    str_replace_all("\\.", "\\\\.")
}

show_row <- function(df, i)
  str_c(df[i,], collapse=",")

fix_names <- function(df, surnames=surname_set(df),
                      not_surnames=c(),
                      firstnames=firstname_set(df),
                      not_firstnames=c()) {
  sns <- alt_regex(surnames)
  # fill empty surnames
  ind <- which(df$surname == "")
  for (i in ind) {
    xs <- str_split(df$firstname[[i]], " +")[[1]]
    ys <- head(xs, -2)
    zs <- tail(xs, 2)
    df$firstname[[i]] <- str_c("", ys, collapse=" ")
    if (length(xs) >= 4) {
      df$surname[[i]] <- str_c("", zs, collapse=" ")
    } else {
      for (z in zs) {
        if (str_detect(z, sns)) {
          df$surname[[i]] <- str_c(
            df$surname[[i]], if (df$surname[[i]] != "") " ", z)
        } else {
          df$firstname[[i]] <- str_c(
            df$firstname[[i]], if (df$firstname[[i]] != "") " ", z)
        }
      }
    }
  }
  # move nots from column "from" to column "to"
  move_nots_from_to <- function(df, from, to, nots,
                                in_reverse) {
    rgx <- alt_regex(nots, FALSE)
    ind <- which(str_detect(df[[from]], rgx))
    r <- function(...)
      if (in_reverse) rev(c(...)) else c(...)
    for (i in ind) {
      xs <- r(str_split(df[[i, from]], " +")[[1]])
      js <- c()
      # for firstnames one has to go in reverse order
      for (j in seq_along(xs)) {
        if (xs[[j]] %in% nots)
          js <- c(js, j)
      }
      j <- max(c(0, js))
      if (j == 0)
        next
      else if (length(js) != j)
        warning("\"", xs[[j]],
                "\" found in position ", j,
                " in row \"", show_row(df, i), "\"")
      keep <- seq_along(xs)[-js]
      df[[i, from]] <- str_c(r(xs[keep]), collapse=" ")
      # add xs[js] (possibly in reverse order)
      # to the column "to"
      # either at the end (when in_reverse=FALSE)
      # or at the beginning (when in_reverse=TRUE)
      df[[i, to]] <- str_c(r(xs[js]), collapse=" ") %>%
        r(df[[i, to]], .) %>%
        str_c(collapse=" ")
    }
    df
  }
  df %>%
    # move not_surnames from surnames to firstnames
    move_nots_from_to("surname", "firstname", not_surnames, FALSE) %>%
    # move not_firstnames from firstnames to surnames
    move_nots_from_to("firstname", "surname", not_firstnames, TRUE)
}

remove_sucesiones <- function(df) {
  ind <- which(str_detect(df$firstname, "^SUC(ESION|.)"))
  # fix surnames in sucesiones first
  for (i in ind) {
    from <- i + 1
    last <- from
    while (df$shares[[last]] == 0)
      last <- last + 1
    last <- last - 1
    empties <- df$surname[from:last] == ""
    if (all(empties)) {
      surname <- df$surname[[i]]
      if (surname == "") {
        print(df[i:last,])
        stop("empty surname")
      } else if (str_detect(surname, " "))
        df$surname[from:last] <- surname
      else {
        warning("sucesion only has one surname: ", surname)
        df$surname[from:last] <- surname
      }
    } else if (any(empties)) {
      surname <- df$surname[from:last][!empties] %>% unique
      if (length(surname) > 1)
        stop("too many surnames to choose from (",
             str_c(surname, collapse=","),
             ") between rows ", from, " and ", last)
      df$surname[from:last] <- surname[[1]]
    }
  }
  sucesiones_ocultas <- c()
  for (i in ind) {
    if (df$shares[[i+1]] == 0)
      df$shares[[i+1]] <- df$shares[[i]]
    else {
      df$firstname[[i]] <- df$firstname[[i]] %>%
        str_replace("^SUC(ESION|.) ", "")
      sucesiones_ocultas <- c(sucesiones_ocultas, i)
      # warning("the share of the first commoner of a ",
      #         "succesion isn't zero but ", df$shares[[i+1]],
      #         "(", df[i+1,], ")")
    }
  }
  ind <- ind[!(ind %in% sucesiones_ocultas)]
  df[-ind,]
}

fix_shares <- function(df) {
  check <- function(i, j) {
    if (j < i) stop(j, " < ", i, "!")
    if (sum(df$shares[i:j]) != df$shares[[i]])
      stop("sum(df$shares[", i, ":", j, "]) == ",
           sum(df$shares[i:j]), " != ", df$shares[[i]],
           " == df$shares[[", i, "]]")
  }
  from <- NA
  for (i in seq_along(df$shares)) {
    if (i == 1 ||
        df$community[[i-1]] != df$community[[i]]) {
      if (df$shares[[i]] == 0)
        warning("community \"", df$community[[i]],
                "\" starts with commoner with zero shares ",
                "at row ", i)
      if (!is.na(from)) {
        check(from, i-1)
        divided_share <- df$shares[[from]] / (i-from)
        df$shares[from:i-1] <- divided_share
      }
      from <- NA
    } else if (is.na(from) && df$shares[[i]] == 0) {
      from <- i-1
    } else if (!is.na(from) && df$shares[[i]] > 0) {
      check(from, i-1)
      divided_share <- df$shares[[from]] / (i-from)
      df$shares[from:(i-1)] <- divided_share
      from <- NA
    }
  }
  df
}

fix_repeated <- function(df) {
  df %>%
    mutate(i = row_number()) %>%
    group_by(community, firstname, surname) %>%
    summarise(right_id = first(right_id),
              shares = sum(shares),
              i = first(i)) %>%
    arrange(i) %>% # reorder rows
    .[,col_names] # reorder columns and get rid of "i"
}

check_table <- function(df) {
  all(df$community != "") &&
    all(df$right_id > 0) &&
    all(df$firstname != "") &&
    all(df$shares > 0)
}

communities_csv <- "communities.csv"
add_commune <- function(df) {
  col_names <- c("commune", "id", "community",
                 "area", "num_houses", "num_commoners",
                 "population", "num_male", "num_female")
  communities <- read_csv(communities_csv, skip=1,
                          col_names=col_names,
                          col_types="ciciiiiii") %>%
    mutate_if(is.character, toupper) %>%
    select(community, commune) %>%
    remove_accents %>%
    rbind(c("TOTORAL", "COPIAPO"),
          c("HUASCOALTINOS", "ALTO DEL CARMEN"))
  left_join(df, communities, by="community")
}

communes_csv <- "communes.csv"
read_commune_csv <- function(filename=communes_csv) {
  col_names <- c("commune", "province", "region",
                 "area", "population", "density", "IDH")
  communes <- read_csv(communes_csv,
                       skip=1,
                       col_names=col_names,
                       col_types="cccdidd") %>%
    remove_accents %>%
    mutate_if(is.character, toupper) %>%
    select(commune, province, region)
}

add_province <- function(df) {
  left_join(df, read_commune_csv(), by="commune")
}

get <- function(v, i, default=NA) {
  if (i <= length(v)) v[[i]]
  else default
}

split_names <- function(df, col, fst, snd) {
  names <- str_split(df[[col]], " +")
  df[[fst]] <- map_chr(names, ~get(., 1, ""))
  df[[snd]] <- map_chr(names, ~get(., 2, ""))
  df
}

split_surnames <- function(df)
  split_names(df, "surname",
              "surname_father",
              "surname_mother")

split_firstnames <- function(df)
  split_names(df, "firstname",
              "firstname1",
              "firstname2")

assign_sex <- function(df) {
  df$sex <- rep("A", n=nrow(df))
  b1 <- endsWith(df$firstname1, "A")
  df$sex[b1] <- "F"
  b2 <- endsWith(df$firstname1, "O")
  df$sex[b2] <- "M"
  b3 <- !b1 & !b2
  b4 <- endsWith(df$firstname2[b3], "A")
  df$sex[b3][b4] <- "F"
  b5 <- endsWith(df$firstname2[b3], "O")
  df$sex[b3][b5] <- "M"
  df
}

read_pdfs <- function() {
  files <- list.files(path="pdfs", pattern="\\.pdf$",
                      full.names=TRUE)
  # skip files in skip
  skipr <- str_c(skip, collapse="|")
  files <- files[!str_detect(files, skipr)]
  rbind(
    files %>% map_dfr(~read_pdf(print(.))),
    manually_extracted_csvs %>%
      map_dfr(~read_community_csv(print(.))) %>%
      mutate_if(is.character, ~str_replace_na(., "")),
    stringsAsFactors=FALSE)
}

extract_pdfs <- function(commoners_df) {
  df <- commoners_df %>%
    mutate_if(is.character, trimws) %>%
    remove_accents %>%
    remove_nonpeople %>%
    fix_compound_names %>%
    remove_abbrev
  # fix specific entries
  df$firstname <- df$firstname %>%
    str_replace("DEL_TRANSITOESPINOZA",
                "DEL_TRANSITO ESPINOZA") %>%
    str_replace("GUILLERMO ENRIQUE DEL_ROSARVICENCIO CORTES",
                "GUILLERMO ENRIQUE DEL_ROSARIO VICENCIO CORTES") %>%
    str_replace("CASIMIRA DEL_CARMEN CLAUDICORTES",
                "CASIMIRA DEL_CARMEN CORTES")
  df$surname <- df$surname %>%
    str_replace("HENRIQUE HENRIQUEZ",
                "HENRIQUEZ HENRIQUEZ")
  w <- which(str_detect(df$firstname, "(?<![A-Z])DEL$"))
  if (w != c(12450))
    stop("missing entry: ", str_c(w, collapse=" "))
  df$firstname[[12450]] <- "ELISA DEL_ROSARIO"
  df$surname[[12450]] <- "ARAYA"
  # this makes the splitting of surnames easier
  # for one case in community "canelilla"
  # but I'm not sure it might have unintended consequences
  df$surname <- df$surname %>%
    str_replace_all(" \\(SUCESION\\)", "")
  # fix surnames and firstnames
  surnames <- surname_set(df)
  firstnames <- firstname_set(df)
  fis <- firstnames[firstnames %in% surnames]
  # sif <- surnames[surnames %in% firstnames]
  # all(fis == sif) # == TRUE
  # "ALFONSO" is the surname of "JACOBITA DE LOURDES"
  # in "LA CHACARILLA"
  not_surnames <- c("AIRES", "ANTONIO", "AQUILES",
                    "ASCENCIO", "BENEDICTO", "BAUTISTA",
                    "BENITO",
                    "DE_LA_CRUZ", "DE_LA_ROSA",
                    "DEL_TRANSITO", "DIOGENES",
                    "ENRIQUE", "ESTER", "GENERAL",
                    "HUMBERTO",
                    "MERCEDES", "MERY", "MIGUEL",
                    "SAN_JUAN", "TOMAS")
  not_firstnames <- c("ALFARO", "ALUCEMA", "ARAYA", "ANAIS",
                      "AVILES", "BACHO", "BERNALES",
                      "BUGUEÑO", "CASTRO", "CATALDO",
                      "CIELO", "GODOY", "HONORES",
                      "LEYTON", "MEDINA", "MOLINA",
                      "OLIVARES", "OLMOS", "ORO",
                      "ORREGO", "PIZARRO", "RIVERA",
                      "ROJAS", "TAPIA", "VICENCIO")
  w1 <- which(surnames %in% not_surnames)
  surnames <- c(surnames[-w1],
                "APEY", "ARGALUZA", "PRYOR", "MILTON")
  w2 <- which(firstnames %in% not_firstnames)
  firstnames <- firstnames[-w2]
  # the first commoner in ALCONES, "YAMILET VANESSA",
  # has no annotated shares. give her 1.
  w <- which(str_detect(df$firstname, "YAMILET VANESSA"))
  if (w[[1]] != 6685)
    stop("missing entry: ", str_c(w, collapse=" "))
  df$shares[[6685]] <- 1L
  # correct community names
  df[df$community == "ALHUEMILLAS LAS PALMAS", "community"] <-
    "ALHUEMILLA LAS PALMAS"
  df[df$community == "CUZ CUZ", "community"] <-
    "CUZCUZ"
  df[df$community == "UCHUMI - DIAGUITA", "community"] <-
    "UCHUMI-DIAGUITAS"
  df[df$community == "ATUNHUAICO", "community"] <-
    "ATUHUAICO"
  df[df$community == "CANELILLA", "community"] <-
    "CANELILLA OVALLE"
  df[df$community == "COIPO O COYUNCAVI", "community"] <-
    "COIPO Y CUYUNCAVI"
  df[df$community == "DAIN Y CORTADERILLA", "community"] <-
    "DAIN Y CORTADERA"
  df[df$community == "FUNDINA SUR", "community"] <-
    "FUNDIDA SUR"
  df[df$community == "FUNDINA NORTE", "community"] <-
    "FUNDIDA NORTE"
  df[df$community == "HUAPI LAS MOLLAQUITAS", "community"] <-
    "HUAPILLAS MOLLAQUITAS"
  df[df$community == "LA CHACARILLA", "community"] <-
    "CHACARILLAS"
  df[df$community == "MACANO", "community"] <-
    "EL MACANO"
  df[df$community == "QUEBRADA DE COLLIGUAYCITO", "community"] <-
    "QUEBRADA DE COLLIGUACITO"
  df[df$community == "CASTILLO MAL PASO Y OTROS", "community"] <-
    "CASTILLO MAL PASO"
  df[df$community == "RINCONADA DE PUNITAQUI", "community"] <-
    "LA RINCONADA DE PUNITAQUI"
  df[df$community == "CARRIZO MENDOZA Y ROMERO", "community"] <-
    "CARRIZO, MENDOZA Y ROMERO"
  df %>%
    # remove communities that appear as commoners
    filter(!startsWith(firstname, "COMUNIDAD")) %>%
    # remove companies that appear as commoners
    filter(!str_detect(surname, "(?<![A-Z])LTDA")) %>%
    fix_names(surnames, not_surnames,
              firstnames, not_firstnames) %>%
    remove_sucesiones %>%
    fix_shares %>%
    fix_repeated %>%
    add_commune %>%
    add_province %>%
    split_firstnames %>%
    split_surnames %>%
    assign_sex %>%
    write_csv("commoners.csv")
}
