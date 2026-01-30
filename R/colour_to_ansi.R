#' Convert plain text colour to ANSI code
#'
#' @param colour colour in plain text ("red", "green", etc.) to convert to ANSI
#'
#' @return string representing provided colour as ANSI encoding
#'
#' @examples
#' colour_to_ansi("red") # gives: "\033[31m"
#' @export
#'
colour_to_ansi <- function(colour) {
  # Note ANSI colour codes
  colour_codes <- list("black" = 30,
                       "red" = 31,
                       "green" = 32,
                       "yellow" = 33,
                       "blue" = 34,
                       "magenta" = 35,
                       "cyan" = 36,
                       "white" = 37)

  # Create ANSI version of colour
  ansi_colour <- paste0("\033[", colour_codes[[colour]], "m")
  return(ansi_colour)
}
