# functions for optimising the alpha and beta parameters

library(ggplot2)
library(dplyr)
library(tidyr)
library(optimx)

# function for generating threshold locations, taken from Selker paper
selker_thresh <- function(a, b, n_resp=7) {
  C <- n_resp
  prop <- 1:(C-1)/C
  Yc <- log(prop / (1-prop))
  a * Yc + b
}

# function for returning MSE of a given alpha and beta for describing threshold locations (actual)
selker_mse <- function(pars = c(a=1, b=0), actual) {
  est <- selker_thresh(a=pars["a"], b=pars["b"], n_resp=length(actual)+1)
  mean((actual - est)**2)
}

# function for optimising alpha and beta for a given vector of threshold locations
# (returns a dataframe based on the output of optimx())
selker_optim <- function(thresh, method="Nelder-Mead") {
  optimx::optimx(
    par = c(a=1, b=0),  # start at equidistant thresholds
    fn = selker_mse,
    method = method,
    actual = thresh
  ) %>%
    as_tibble(rownames="method") %>%
    mutate(
      # if multiple methods, will space them out in the plot
      method = factor(method),
      y = 0 - (as.numeric(method)-1) * 0.01,
      x = list(selker_thresh(a, b, n_resp = length(thresh)+1))
    ) %>%
    group_by(method) %>%
    mutate(actual = list(thresh)) %>%
    unnest(cols=c(x, actual))
}

# selker_optim(sort(selker_thresh(0.5, 0.25) + rnorm(6, 0, 0.5)))

# function for plotting the output of selker_optim() on a latent distribution
plot_latent <- function(opt_res) {
  xlims <- range(c(opt_res$x), c(opt_res$actual), -5, 5)
  
  ggplot() +
    geom_function(aes(alpha="Latent Density"), fun = dnorm, size=1.5, colour="darkgrey") +
    geom_vline(aes(alpha="Actual Locations", xintercept = actual), data.frame(actual=unique(opt_res$actual)), size=1, show.legend=FALSE) +
    geom_point(aes(x, y, alpha="Estimated Locations"), colour="red", size=3, data=opt_res) +
    lims(x = xlims) +
    labs(x = "Value", y = "Density", colour=NULL) +
    scale_alpha_manual(
      name = NULL,
      values = c(1, 1, 1),
      breaks = c("Latent Density", "Actual Locations", "Estimated Locations"),
      guide = guide_legend(override.aes = list(
        linetype = c(1, 1, 0),
        size = c(1.5, 0.75, 3),
        shape = c(NA, NA, 16),
        color = c("darkgrey", "black", "red")
      ))
    ) +
    theme(legend.position = "bottom")
}

# sort(selker_thresh(0.5, 0.25) + rnorm(6, 0, 0.5)) %>%
#   selker_optim() %>%
#   plot_latent()

# function for plotting the distortion that would occur based on the selker_optim() results
plot_distort <- function(opt_res) {
  opt_res %>%
    ggplot(aes(actual, x, group=method)) +
    geom_abline(aes(alpha="Ideal", intercept=i, slope=s), data=data.frame(i = 0, s = 1), linetype = "dashed", colour = "darkgrey", size=1, show.legend = FALSE) +
    geom_line(aes(alpha="Observed"), size=1.5) +
    geom_point(size=5) +
    labs(x = "Actual Location", y = "Estimated Location") +
    scale_alpha_manual(
      name = NULL,
      values = c(1, 1),
      breaks = c("Observed", "Ideal"),
      guide = guide_legend(override.aes = list(
        linetype = c("solid", "dashed"),
        color = c("black", "darkgrey")
      ))
    ) +
    theme(legend.position = "bottom")
}

# sort(selker_thresh(0.5, 0.25) + rnorm(6, 0, 0.5)) %>%
#   selker_optim() %>%
#   plot_distort()
