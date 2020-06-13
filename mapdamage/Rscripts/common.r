options(error = function(...) {
  traceback(2)
  quit(status = 1)
})


abort <- function(...) {
  cat("Error: ", sprintf(...), "\n", sep = "", file = stderr())
  quit(status = 2)
}


getArgument <- function(key, type = as.character) {
  if (!exists(".KWARGS")) {
    .KWARGS <<- list()

    for (.arg in commandArgs(TRUE)) {
      .idx <- regexpr("=", .arg, fixed = TRUE)
      if (.idx != -1) {
        .key <- toupper(sub("=.*", "", .arg))
        .value <- sub("[^=]*=", "", .arg)

        .KWARGS[[.key]] <<- .value
      }
    }
  }

  key <- toupper(key)
  if (key %in% names(.KWARGS)) {
    value <- .KWARGS[[key]]
    parsed <- type(value)
    if (is.na(parsed)) {
      abort("Invalid value for '%s': '%s'", key, value)
    }

    return(parsed)
  } else {
    abort("Missing argument '%s'", key)
  }
}


read.mapDamage.table <- function(file) {
  return(read.table(file = file, sep = "\t", header = TRUE, as.is = TRUE, check.names = FALSE))
}


iterLibraries <- function(tbl, min_libraries=0) {
    result <- list()
    libraries <- unique(tbl[, c("Sample", "Library")])
    for (row.i in rownames(libraries)) {
        row <- libraries[row.i,]
        result <- append(result, list(as.list(row)))
    }

    if (length(result) < min_libraries) {
        return(list())
    }

    return(result)
}
