Test code
================

## Fitting the models

As detailed in the manuscript, we derive five nested versions of the
model:

-   *B* arbitrary, which can be selected using `model = "full"`
-   *B = D(d) + vw^T*, which can be selected using `model = "diag_vwt"`
-   *B = D(d) + vv^T*, which can be selected using `model = "diag_vvt"`
-   *B = D(d) + ^T*, which can be selected using `model = "diag_a11t"`

We also have two goal functions, attempting to minimize either the sum
of squared deviations (OLS, selected using `goal = "SSQ"`), or the
weighted sum of the squared deviations (WLS, selected using
`goalf = "WLS"`). Options also allow to :

-   `pars = "NULL` by default, the parameters are initialized with *B =
    I* for all models. The user can provide alternative parameters
    instead.
-   `skipEM = FALSE` skip the iterative algorithm and perform the
    numerical search directly.
-   `plot_results = TRUE` plot the results at the end of the
    calculation.

For example, to fit the data by Kuebbing et al. (2015), non-native
species, using OLS and the simplified model *B = D(d) + vw^T*, one can
invoke:

``` r
source("general.R")
run_model(datafile = "../data/kuebbing_2015_non_natives.csv", # location of the data
          model = "diag_vwt", # one of "diag_a11t", "diag_vvt", "diag_vwt", or "full"
          goalf = "SSQ", # one of SSQ, WLS
          pars = NULL, # pre-computed parameters; otherwise use the identity matrix
          skipEM = FALSE, # skip the iterative algorithm
          plot_results = TRUE # plot the results as boxplots
          )
```

    ## [1] "iterative step 1 -> 67.0062462911419"
    ## [1] "iterative step 2 -> 58.1125868909348"
    ## [1] "iterative step 3 -> 52.831579796576"
    ## [1] "iterative step 4 -> 48.1614107930199"
    ## [1] "iterative step 5 -> 43.7510941997679"
    ## [1] "iterative step 6 -> 40.7931051934239"
    ## [1] "iterative step 7 -> 38.7693436324155"
    ## [1] "iterative step 8 -> 37.3940743308755"
    ## [1] "iterative step 9 -> 36.3404144951679"
    ## [1] "iterative step 10 -> 35.9586264383472"
    ## [1] "iterative step 11 -> 34.918386145092"
    ## [1] "iterative step 12 -> 34.495479516673"
    ## [1] "iterative step 13 -> 34.0632653666747"
    ## [1] "iterative step 14 -> 33.7667109953849"
    ## [1] "iterative step 15 -> 33.3999522502816"
    ## [1] "iterative step 16 -> 33.23408137218"
    ## [1] "iterative step 17 -> 32.9573872515676"
    ## [1] "iterative step 18 -> 32.7704882115572"
    ## [1] "iterative step 19 -> 32.5549861467456"
    ## [1] "iterative step 20 -> 32.3779777912159"
    ## [1] "iterative step 21 -> 32.2017947003788"
    ## [1] "iterative step 22 -> 32.0677438786954"
    ## [1] "iterative step 23 -> 32.0020931017848"
    ## [1] "iterative step 24 -> 31.9375212099102"
    ## [1] "iterative step 25 -> 31.7738898224771"
    ## [1] "numerical search"
    ## [1] 29.90851
    ## [1] 29.90851
    ## [1] 29.90851
    ## [1] 29.90851
    ## [1] 29.90851
    ## [1] 29.90851

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

## Out-of-fit predictions

The model can be fitted using part of the data, and then the parameters
used to predict the whole data set, including the data that was not used
to fit the model. To implement a simple Leave-One-Out approach, one can
call `run_model_LOO`, with an extra parameter, `LOO_row_num`. This is
the row number of the data to be excluded. All communities of the same
type will be excluded as well. For example, in the file
`kuebbing_2015_non_natives.csv`, row 22 contains a replicate of the
community `as + fa`. By selecting this row, all communities of the same
type are excluded from the fit. When plotting, the out-of-fit prediction
has colored boxplots, while the in-fit data is plotted in white:

``` r
source("general.R")
run_model_LOO(datafile = "../data/kuebbing_2015_non_natives.csv", # location of the data
          model = "diag_vwt", # one of "diag_a11t", "diag_vvt", "diag_vwt", or "full"
          goalf = "WLS", # one of SSQ, WLS
          LOO_row_num = 22, # exclude all replicates of as + fa from the fit
          pars = NULL, # pre-computed parameters; otherwise use the identity matrix
          skipEM = FALSE, # skip the iterative algorithm
          plot_results = TRUE # plot the results as boxplots
          )
```

    ## [1] "iterative step 1 -> 644.094376223836"
    ## [1] "iterative step 2 -> 576.48941099368"
    ## [1] "iterative step 3 -> 422.402096274447"
    ## [1] "iterative step 4 -> 422.402096274447"
    ## [1] "iterative step 5 -> 313.974067213791"
    ## [1] "iterative step 6 -> 313.974067213791"
    ## [1] "iterative step 7 -> 313.974067213791"
    ## [1] "iterative step 8 -> 313.974067213791"
    ## [1] "iterative step 9 -> 313.974067213791"
    ## [1] "iterative step 10 -> 313.974067213791"
    ## [1] "iterative step 11 -> 313.974067213791"
    ## [1] "iterative step 12 -> 313.974067213791"
    ## [1] "iterative step 13 -> 313.974067213791"
    ## [1] "iterative step 14 -> 313.974067213791"
    ## [1] "iterative step 15 -> 313.974067213791"
    ## [1] "iterative step 16 -> 313.974067213791"
    ## [1] "iterative step 17 -> 313.974067213791"
    ## [1] "iterative step 18 -> 313.974067213791"
    ## [1] "iterative step 19 -> 313.974067213791"
    ## [1] "iterative step 20 -> 313.974067213791"
    ## [1] "iterative step 21 -> 313.974067213791"
    ## [1] "iterative step 22 -> 313.974067213791"
    ## [1] "iterative step 23 -> 313.974067213791"
    ## [1] "iterative step 24 -> 313.974067213791"
    ## [1] "iterative step 25 -> 313.974067213791"
    ## [1] "numerical search"
    ## [1] 256.6098
    ## [1] 256.6098
    ## [1] 256.6098
    ## [1] 256.6098
    ## [1] 256.6098
    ## [1] 256.6098

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
