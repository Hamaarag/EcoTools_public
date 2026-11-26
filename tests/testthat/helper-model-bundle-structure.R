expect_unique_bundle_names <- function(
  bundle
) {
  duplicated_names <- unique(
    names(bundle)[
      duplicated(
        names(bundle)
      )
    ]
  )

  testthat::expect_length(
    duplicated_names,
    0L,
    info = if (length(duplicated_names) > 0L) {
      paste0(
        "Duplicated bundle element name(s): ",
        paste(
          duplicated_names,
          collapse = ", "
        )
      )
    } else {
      NULL
    }
  )
}
