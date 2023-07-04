#' Availability calculation for a point transect
#'
#' @param delta bin width
#' @param midpoint bin midpoint
#' @param max_distance maximum distance of survey (truncation distance)
#'
#' @return numeric
availability_fn <- function(delta, midpoint, max_distance) {
  (2 * delta * midpoint)/max_distance^2
}

#' Probability closest animal is in a given bin, given a provided bin and number of individuals in a photo
#'
#' @param bin_start start point of bin
#' @param bin_end end point of bin
#' @param max_distance maximum distance of survey (truncation distance)
#' @param n number of individuals in photo
#'
#' @return numeric
pr_closest <- function(bin_start, bin_end, max_distance, n) {
  if(bin_start == 0) {
    prob_closer <- 0
  } else {
    closer_midpoint <- (bin_start)/2
    prob_closer <- 1 - (1 - availability_fn(bin_start, closer_midpoint, max_distance = max_distance))^n
  }
  if(bin_end == max_distance) {
    prob_further <- 0
  } else {
    further_midpoint <- (bin_end+max_distance)/2
    further_delta <- max_distance-bin_end
    prob_further <- availability_fn(further_delta, further_midpoint, max_distance = max_distance)
  }

  # Combined probability
  # probability that the closest indiv is in this bin is 1 minus the probability that the closest was in a closer bin, minus the probability that all individuals are in further bins
  pi_m <- 1 - (prob_closer + prob_further^n)

  return(pi_m)
}

get.beta.param<- function(mu,sdev){
  # get parameters of beta by method of moments
  a<- mu*((mu*(1-mu))/sdev^2 - 1)
  b<- (1-mu)*((mu*(1-mu))/sdev^2 - 1)
  if(is.na(a) || is.na(b)) stop("Invalid Beta parameters")
  else return(list(a=a,b=b))
}
