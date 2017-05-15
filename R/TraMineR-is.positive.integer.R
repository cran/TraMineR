# Check if an object is a positive integer.
#
# @param x Object to test.
#
# @return TRUE if the object is a positive integer. FALSE if not.
#
# @author Pierre-Alexandre Fonta (2016-2017)

is.positive.integer <- function(x) {
  return(is.a.number(x) && x >= 0 && x %% 1 == 0)
}
