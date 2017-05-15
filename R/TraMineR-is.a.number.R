# Check if an object is a number.
#
# @param x Object to test.
#
# @return TRUE if the object is a number. FALSE if not.
#
# @author Pierre-Alexandre Fonta (2016-2017)

is.a.number <- function(x) {
  return(length(x) == 1 && is.numeric(x))
}
