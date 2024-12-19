Test code
================

## Fitting the models

As detailed in the manuscript, we derive five nested versions of the
model:

-   *B* arbitrary, which can be selected using `model = "full"`
-   *B = D(d) + vw^T*, which can be selected using `model = "diag_vwt"`
-   *B = D(d) + vv^T*, which can be selected using `model = "diag_vvt"`
-   *B = D(d) + alpha 11^T*, which can be selected using
    `model = "diag_a11t"`

We can also choose between two goal functions, attempting to minimize
either the sum of squared deviations (OLS, selected using
`goal = "SSQ"`), or the weighted sum of the squared deviations (WLS,
selected using `goalf = "WLS"`).

Other options allow to :

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
res <- run_model(datafile = "../data/kuebbing_2015_non_natives.csv", # location of the data
          model = "diag_vwt", # one of "diag_a11t", "diag_vvt", "diag_vwt", or "full"
          goalf = "SSQ", # one of SSQ, WLS
          pars = NULL, # pre-computed parameters; otherwise use the identity matrix
          skipEM = FALSE, # skip the iterative algorithm
          plot_results = TRUE # plot the results as boxplots
          )
```

    ## [1] "iterative step 1 -> 990.888949636642"
    ## [1] "iterative step 2 -> 629.642704358945"
    ## [1] "iterative step 3 -> 504.921563372439"
    ## [1] "iterative step 4 -> 418.84452081321"
    ## [1] "iterative step 5 -> 376.117486397484"
    ## [1] "iterative step 6 -> 354.227370774127"
    ## [1] "iterative step 7 -> 335.004319289932"
    ## [1] "iterative step 8 -> 315.039958827847"
    ## [1] "iterative step 9 -> 299.682234330772"
    ## [1] "iterative step 10 -> 287.845924104116"
    ## [1] "iterative step 11 -> 273.910497474283"
    ## [1] "iterative step 12 -> 267.045959348568"
    ## [1] "iterative step 13 -> 263.598770503308"
    ## [1] "iterative step 14 -> 254.381567707855"
    ## [1] "iterative step 15 -> 254.27602348795"
    ## [1] "iterative step 16 -> 248.602828117017"
    ## [1] "iterative step 17 -> 245.998972247379"
    ## [1] "iterative step 18 -> 244.141377327376"
    ## [1] "iterative step 19 -> 240.220213084485"
    ## [1] "iterative step 20 -> 240.061305491651"
    ## [1] "iterative step 21 -> 239.159615031142"
    ## [1] "iterative step 22 -> 237.378566693952"
    ## [1] "iterative step 23 -> 237.109223060851"
    ## [1] "iterative step 24 -> 234.948355492571"
    ## [1] "iterative step 25 -> 234.42701216238"
    ## [1] "numerical search"
    ## [1] 212.7557
    ## [1] 212.7556
    ## [1] 212.7556
    ## [1] 212.7556
    ## [1] 212.7556
    ## [1] 212.7556

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
str(res)
```

    ## List of 8
    ##  $ data_name    : chr "kuebbing_2015_non_natives"
    ##  $ observed     : num [1:140, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : NULL
    ##   .. ..$ : chr [1:4] "as" "fa" "la" "po"
    ##  $ predicted    : num [1:140, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : NULL
    ##   .. ..$ : chr [1:4] "as" "fa" "la" "po"
    ##  $ variances    : num [1:140, 1:4] 0 0 0 0 0 0 0 0 0 0 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : NULL
    ##   .. ..$ : chr [1:4] "as" "fa" "la" "po"
    ##  $ B            : num [1:4, 1:4] 0.3984 0.314 0.3255 0.3134 0.0758 ...
    ##  $ pars         : num [1:12] 11.95 54.39 13.3 8.3 -1.51 ...
    ##  $ goal_function: num 213
    ##  $ goal_type    : chr "SSQ"

## Out-of-fit predictions

The model can be fitted using part of the data, and then the parameters
used to predict the whole data set, including the data that was not used
to fit the model. To implement a simple Leave-One-Out approach, one can
call `run_model_LOO`, with an extra parameter, `LOO_row_num`. This is
the row number of the data to be excluded. All communities of the same
type will be excluded as well. For example, in the file
`kuebbing_2015_non_natives.csv`, row 22 contains a replicate of the
community `as + fa`. By selecting this row, all communities of the same
type are excluded from the fit.

When plotting, the out-of-fit prediction is reported using colored
boxplots, while the in-fit data is plotted in white:

``` r
source("general.R")
res <- run_model_LOO(datafile = "../data/kuebbing_2015_non_natives.csv", # location of the data
          model = "diag_vwt", # one of "diag_a11t", "diag_vvt", "diag_vwt", or "full"
          goalf = "WLS", # one of SSQ, WLS
          LOO_row_num = 22, # exclude all replicates of as + fa from the fit
          pars = NULL, # pre-computed parameters; otherwise use the identity matrix
          skipEM = FALSE, # skip the iterative algorithm
          plot_results = TRUE # plot the results as boxplots
          )
```

    ## [1] "iterative step 1 -> 621.808110694175"
    ## [1] "iterative step 2 -> 621.808110694175"
    ## [1] "iterative step 3 -> 599.913139473817"
    ## [1] "iterative step 4 -> 384.2634314538"
    ## [1] "iterative step 5 -> 357.417389231598"
    ## [1] "iterative step 6 -> 311.892906469798"
    ## [1] "iterative step 7 -> 311.892906469798"
    ## [1] "iterative step 8 -> 289.171194810318"
    ## [1] "iterative step 9 -> 285.761608937196"
    ## [1] "iterative step 10 -> 285.574791134673"
    ## [1] "iterative step 11 -> 285.574791134673"
    ## [1] "iterative step 12 -> 285.574791134673"
    ## [1] "iterative step 13 -> 285.574791134673"
    ## [1] "iterative step 14 -> 285.574791134673"
    ## [1] "iterative step 15 -> 285.574791134673"
    ## [1] "iterative step 16 -> 285.574791134673"
    ## [1] "iterative step 17 -> 285.574791134673"
    ## [1] "iterative step 18 -> 285.574791134673"
    ## [1] "iterative step 19 -> 285.574791134673"
    ## [1] "iterative step 20 -> 285.574791134673"
    ## [1] "iterative step 21 -> 285.574791134673"
    ## [1] "iterative step 22 -> 285.574791134673"
    ## [1] "iterative step 23 -> 285.574791134673"
    ## [1] "iterative step 24 -> 285.574791134673"
    ## [1] "iterative step 25 -> 285.574791134673"
    ## [1] "numerical search"
    ## [1] 260.4811
    ## [1] 260.4811
    ## [1] 260.4811
    ## [1] 260.4811
    ## [1] 260.4811
    ## [1] 260.4811

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
