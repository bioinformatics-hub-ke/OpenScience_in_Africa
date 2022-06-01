

#function for extracting pubmed information
extract_pubMed_data <- function (pubMed_query, 
                                 batch_size = 1000, 
                                 getKeywords = FALSE, 
                                 affi_regex_exclude = NULL) 
{
  ptm <- proc.time()
  my.idlist <- get_pubmed_ids(pubMed_query)
  record.num <- my.idlist$Count
  my.seq <- seq(1, as.numeric(my.idlist$Count), by = batch_size)
  pubmed.data <- lapply(my.seq, (function(ret.start) {
    batch.xml <- NULL
    message(paste("round #", which(my.seq == ret.start), 
                  " of ", length(my.seq), " ", sep = ""), appendLF = FALSE)
    while (is.null(batch.xml)) {
      batch.xml <- tryCatch({
        tmp.idlist <- get_pubmed_ids(pubMed_query)
        fetch_pubmed_data(tmp.idlist, retstart = ret.start, 
                          retmax = batch_size)
      }, error = function(e) {
        NULL
      })
    }
    record.list <- easyPubMed::articles_to_list(batch.xml)
    xtracted.data <- lapply(1:length(record.list), (function(i) {
      if (length(record.list) > 60) {
        custom.seq <- as.integer(seq(1, length(record.list), 
                                     length.out = 50))
        if (i %in% custom.seq) {
          message(".", appendLF = FALSE)
        }
      } else {
        message(".", appendLF = FALSE)
      }
      if (as.numeric(installed.packages()["easyPubMed","Version"]) > 2.5){
        tmp.record <- tryCatch(article_to_df(pubmedArticle = record.list[[i]], 
                                             getKeywords = getKeywords, 
                                             autofill = TRUE, max_chars = 10), 
                               error = function(e) {
                                 NULL
                               })
        required.cols <- c("title", "year", "journal", 
                           "lastname", "firstname", "address", "email", "pmid", "month","day")
        
        
      } else {
        tmp.record <- tryCatch(article_to_df(pubmedArticle = record.list[[i]], 
                                             autofill = TRUE, max_chars = 10), 
                               error = function(e) {
                                 NULL
                               })
        required.cols <- c("title", "year", "journal", 
                           "lastname", "firstname", "address", "email","pmid", "month", "day")
        
      }
      
      if (!is.null(tmp.record)) {
        out.record <- data.frame(matrix(NA, nrow = nrow(tmp.record), 
                                        ncol = length(required.cols)))
        colnames(out.record) <- required.cols
        match.cols <- colnames(tmp.record)[colnames(tmp.record) %in% 
                                             required.cols]
        out.record[, match.cols] <- tmp.record[, match.cols]
      } else {
        out.record <- NULL
      }
      out.record
    }))
    xtracted.data <- do.call(rbind, xtracted.data)
    message(" Filtering... ", appendLF = FALSE)
    xtracted.data <- xtracted.data[!is.na(xtracted.data$address), 
    ]
    if (!is.null(affi_regex_exclude)) {
      xtracted.data <- xtracted.data[regexpr(affi_regex_exclude, 
                                             xtracted.data$address,
                                             ignore.case = TRUE) < 0, ]
    }
    message("done!", appendLF = TRUE)
    xtracted.data
  }))
  stop.watch <- proc.time() - ptm
  pubmed.data <- do.call(rbind, pubmed.data)
  out.data <- list()
  out.data$params <- list()
  out.data$params$query_string <- pubMed_query
  out.data$params$pubMed_id_list <- my.idlist
  out.data$params$batch_size <- batch_size
  out.data$params$affi_regex_exclude <- affi_regex_exclude
  out.data$params$timing <- stop.watch
  out.data$data <- pubmed.data
  return(out.data)
}



#DESIGNING FUNCTION 
#fetch Pmid from entrez
fetch_PMID_data <- function(pmids, format = "xml", encoding = "UTF-8", verbose = TRUE) 
{
  
  # Hardcoded
  pmids.per.batch = 100
  
  # custom f(x)
  query_and_fetch <- function(x, format = format, encoding = encoding) {
    q0 <- paste(paste0(x, "[PMID]"), collapse = " OR ")
    
    q1 <- easyPubMed::get_pubmed_ids(pubmed_query_string = q0)
    r1 <- easyPubMed::fetch_pubmed_data(pubmed_id_list = q1, format = format, encoding = encoding)
    r2 <- easyPubMed::articles_to_list(pubmed_data = r1, encoding = encoding, simplify = FALSE)
    
    nm2 <- lapply(r2, function(xx) { 
      y <- xx[1]; y <- gsub("</PMID>.+$", "</PMID>", y); custom_grep(y, "PMID", format = "list")[[1]]
    })
    
    
    nm2 <- as.character(do.call(c, nm2))
    
    # Order
    nuORD <- as.numeric(sapply(x, function(z) {
      which(nm2 == z)
    }))      
    r3 <- r2[nuORD]
    return(r3)
  }
  
  # standardize input
  if(is.character(pmids)) {
    pmids <- pmids[!is.na(pmids)]
  } else if (is.list(pmids)) {
    pmids <- do.call(c, pmids)
    pmids <- as.character(pmids)
    pmids <- pmids[!is.na(pmids)]
  }
  
  # split ID in small batches
  out <- list()
  ti <- Sys.time() - 2
  
  if (length(pmids) < pmids.per.batch) {
    
    out <- list(query_and_fetch(pmids, format = format, encoding = encoding))
    
  } else {
    all_i <- seq(1, (length(pmids) - 1), by = pmids.per.batch)
    all_i <- c(all_i, (length(pmids) + 1))
    
    for(j in 1:(length(all_i) - 1)) {
      
      if (verbose)
        message(".", appendLF = FALSE)
      
      i0 <- all_i[j]
      i1 <- (all_i[(j+1)] - 1)
      tdf <- as.numeric(difftime(time1 = Sys.time(), time2 =  ti, units = "sec"))
      
      # check time
      if (tdf < 1) {
        Sys.sleep(time = (1 - tdf))
      } 
      
      # generate query string
      tmp <- pmids[i0:i1]
      tmp2 <- query_and_fetch(x = tmp, format = format, encoding = encoding)
      names(tmp2) <- NULL
      
      ti <- Sys.time()
      out[[length(out) + 1]] <- tmp2
    }
  }
  
  # loop over and recompose
  OUT <- list()
  for(l1 in 1:length(out)) {
    TMP <- out[[l1]]
    for(l2 in 1:length(TMP)) {
      OUT[[length(OUT) + 1]] <- TMP[[l2]]
    }
  }
  
  if(verbose) {
    message("", appendLF = TRUE)
    message("Done!", appendLF = TRUE)
  }
  
  return(OUT)
}   



#extract ids(doi,pmcid,pmid) from pubmed 
extract_article_ids <- function(pubmed_data_list) {
  
  x <- pubmed_data_list
  
  # Tmp f(x)
  my_grep <- function(X, idtype = "pubmed") {
    
    myPAT <- paste0("^.*", idtype)
    out <- tryCatch({
      if (!grepl(myPAT, X)) {
        Y <- NA
      } else {
        Y <- sub(myPAT, "", X)
        Y <- sub("</ArticleId>.*$", "", Y)
        Y <- sub("^.*>", "", Y)
      }
      Y
    }, error = function(e) NA)
    
    return(out)  
  }
  
  # dddd
  out <- lapply(x, function(xx) {
    y <- sub("</ArticleIdList>.+", "</ArticleIdList>", xx)
    y <- sub("^.+(<ArticleIdList.+$)", "\\1", y)
    
    data.frame(
      PMID = my_grep(X = y, idtype = "pubmed"),
      DOI = my_grep(X = y, idtype = "doi"), 
      PMCID = my_grep(X = y, idtype = "pmc"), 
      stringsAsFactors = FALSE)
  })
  
  # Return
  out <- do.call(rbind, out)
  return(out)
} 

