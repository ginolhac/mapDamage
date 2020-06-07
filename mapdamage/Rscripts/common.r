on_error <- function(...) {
  traceback(2)
  quit(status = 1)
}

options(error = on_error)


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
