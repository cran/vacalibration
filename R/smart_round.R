#' Round and maintain a target sum
#'
#' Rounds a vector to the specified number of decimal places and maintains the sum it had before rounding.
#'
#' @param x Numeric vector.
#'
#' @param target_sum Numeric. The target sum to be maintained after rounding. Default is \code{NULL} which sets \code{target_sum=sum(x)}.
#'
#' @param digits Positive integer. Indicates the number of decimal places to be used.
#'
#' @return Numeric vector.
#'
#' @examples
#'
#' x = rep(1/3, 3)
#' round(x, 2)
#' smart_round(x, 1, 2)
#'
#' @export
smart_round <- function(x, target_sum, digits = 0) {

  if(is.null(target_sum)){target_sum = sum(x)}

  # 1. Normalize the input vector into proportions (so they sum exactly to 1)
  # This allows the user to input raw frequencies, weights, or probabilities.
  proportions <- x / sum(x)

  # 2. Determine the multiplier based on the desired decimal precision
  multiplier <- 10^digits

  # 3. Calculate the scaled target_sum as an integer
  # (e.g., target_sum 1 with 2 digits becomes 100 integer units)
  target_sum_scaled <- round(target_sum * multiplier)

  # 4. Scale proportions to the integer target_sum
  scaled <- proportions * target_sum_scaled
  floored <- floor(scaled)
  remainders <- scaled - floored

  # 5. Calculate the shortfall in integer units
  shortfall <- target_sum_scaled - sum(floored)

  # 6. Distribute the shortfall to the items with the largest remainders
  if (shortfall > 0) {
    # order() handles ties by maintaining original vector order
    top_indices <- order(remainders, decreasing = TRUE)[1:shortfall]
    floored[top_indices] <- floored[top_indices] + 1
  }

  # 7. Divide by the multiplier to restore the decimal places
  return(floored / multiplier)
}
