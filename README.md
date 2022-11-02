# selker

This Shiny app demonstrates [Selker et al.'s (2019)](https://doi.org/10.3758/s13428-019-01231-3) method of describing an arbitrary number of threshold locations with just two parameters, showing how well the method can describe thresholds of known locations.

## Running the App

You can run the app locally or remotely:

* Run it locally in R, with `shiny::runGitHub(repo="selker", username="JackEdTaylor")`
* Go to https://jackt.shinyapps.io/selker/
* Go to https://rstudio-connect.psy.gla.ac.uk/selker/
* Go to https://shiny.psy.gla.ac.uk/jackt/selker/

## More on the Method

Models estimating threshold locations on latent distributions are very useful for analyses of ordinal ratings. A common way to account for participant or item variablility in response biases is to estimate thresholds separately for each participant or item, but for large rating scales this requires a large number of parameters to be estimated. [Selker et al.](https://doi.org/10.3758/s13428-019-01231-3) suggest a method for estimating the locations of any number of thresholds which only requires the estimation of two parameters: $a$ (describing scale) and $b$ (descrbing shift). Applying this to estimate per-participant and per-item thresholds can be more efficient.

Threshold $\lambda_c$ in position $c$, when there are $C$ ordered regions in the latent distribution, has location:

$$\gamma_c = log ( \frac{c/C}{1-c/C} )$$

$$\lambda_c = a \gamma_c + b$$
