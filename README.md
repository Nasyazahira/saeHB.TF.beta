
<!-- README.md is generated from README.Rmd. Please edit that file -->

# saeHB.TF.beta

<!-- badges: start -->
<!-- badges: end -->

`saeHB.TF.beta` provides several functions for area and subarea level of
small area estimation under Twofold Subarea Level Model using
hierarchical Bayesian (HB) method with Beta distribution for variables
of interest. Some dataset simulated by a data generation are also
provided. The ‘rstan’ package is employed to obtain parameter estimates
using STAN.

## Function

## Installation

You can install the development version of saeHB.TF.beta from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
# devtools::install.github("Nasyazahira/saeHB.TF.beta")
```

## Example

Here is a basic example of using the **betaTF** function to make
estimates based on sample data in this package

``` r
library(saeHB.TF.beta)

#Load Dataset
data(dataBeta) #for dataset with nonsampled subarea use dataBetaNS

#Fitting model
fit <- betaTF(y~X1+X2, area="codearea", weight="w", data=dataBeta)
#> 
#> SAMPLING FOR MODEL 'saeHB_TF_beta' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0.000157 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 1.57 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 1.515 seconds (Warm-up)
#> Chain 1:                1.403 seconds (Sampling)
#> Chain 1:                2.918 seconds (Total)
#> Chain 1:
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> 
#> SAMPLING FOR MODEL 'saeHB_TF_beta' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 8.7e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.87 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 1.066 seconds (Warm-up)
#> Chain 1:                0.783 seconds (Sampling)
#> Chain 1:                1.849 seconds (Total)
#> Chain 1:
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> 
#> SAMPLING FOR MODEL 'saeHB_TF_beta' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 6.6e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.66 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:   1 / 1000 [  0%]  (Warmup)
#> Chain 1: Iteration: 100 / 1000 [ 10%]  (Warmup)
#> Chain 1: Iteration: 200 / 1000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 300 / 1000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 400 / 1000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 500 / 1000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 501 / 1000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 600 / 1000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 700 / 1000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 800 / 1000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 900 / 1000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 1000 / 1000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.896 seconds (Warm-up)
#> Chain 1:                1.056 seconds (Sampling)
#> Chain 1:                1.952 seconds (Total)
#> Chain 1:
```

Extract subarea mean estimation

``` r
fit$Est_sub
#>             Mean         SD       2.5%       25%       50%       75%     97.5%
#> mu[1]  0.9028999 0.09360662 0.63963221 0.8807558 0.9355760 0.9609809 0.9897710
#> mu[2]  0.8884830 0.08929324 0.66040068 0.8544526 0.9092362 0.9536080 0.9870566
#> mu[3]  0.8877850 0.08808811 0.64227430 0.8418483 0.9132770 0.9518893 0.9853221
#> mu[4]  0.8821074 0.08887543 0.67103155 0.8437779 0.9071352 0.9438826 0.9812516
#> mu[5]  0.6248557 0.19787880 0.19979730 0.5048811 0.6407873 0.7761079 0.9422180
#> mu[6]  0.8254342 0.12453690 0.49820070 0.7659647 0.8580730 0.9116334 0.9765322
#> mu[7]  0.8791863 0.09831153 0.60761287 0.8433229 0.9052431 0.9484217 0.9838472
#> mu[8]  0.8574591 0.10633678 0.56813340 0.8017131 0.8886431 0.9344778 0.9740133
#> mu[9]  0.8822929 0.08823396 0.64495522 0.8516709 0.9043668 0.9437214 0.9840131
#> mu[10] 0.9013078 0.08716471 0.68704143 0.8707152 0.9270385 0.9615326 0.9899438
#> mu[11] 0.8684039 0.10602093 0.57582753 0.8274586 0.8974754 0.9428251 0.9849724
#> mu[12] 0.8902321 0.09286219 0.64568821 0.8512685 0.9206774 0.9548377 0.9894241
#> mu[13] 0.6999023 0.16182222 0.35419889 0.6023424 0.7278719 0.8119006 0.9444410
#> mu[14] 0.8225585 0.12738247 0.48741259 0.7616193 0.8608365 0.9170292 0.9749569
#> mu[15] 0.4426269 0.21761025 0.08868534 0.2622981 0.4338759 0.5981972 0.8746677
#> mu[16] 0.8359049 0.11167703 0.58295991 0.7824194 0.8604786 0.9178787 0.9763436
#> mu[17] 0.9123215 0.08181373 0.69769857 0.8851462 0.9376697 0.9667158 0.9898391
#> mu[18] 0.7875051 0.14256993 0.45464683 0.7119661 0.8170070 0.8968464 0.9689241
#> mu[19] 0.3879337 0.22018559 0.08128812 0.2032067 0.3517219 0.5474827 0.8256371
#> mu[20] 0.8134362 0.13450009 0.44528979 0.7570132 0.8531842 0.9116585 0.9741450
#> mu[21] 0.6038188 0.18574257 0.20891634 0.4743283 0.6123170 0.7514118 0.9105432
#> mu[22] 0.7897103 0.14503904 0.43212812 0.7245095 0.8278949 0.8953567 0.9656462
#> mu[23] 0.7237086 0.16307546 0.33727326 0.6354550 0.7512128 0.8443797 0.9535892
#> mu[24] 0.6439534 0.17701854 0.25475546 0.5240780 0.6603727 0.7811855 0.9237107
#> mu[25] 0.4795057 0.21697879 0.12890988 0.3088672 0.4501044 0.6483116 0.9022565
#> mu[26] 0.4392388 0.19763481 0.11815828 0.2857771 0.4134936 0.5744910 0.8547213
#> mu[27] 0.6399657 0.18246131 0.22906604 0.5252883 0.6600119 0.7851923 0.9191572
#> mu[28] 0.9237816 0.07174340 0.71731590 0.9051032 0.9447943 0.9697162 0.9933175
#> mu[29] 0.8309203 0.11819816 0.52677042 0.7733554 0.8560810 0.9195697 0.9776610
#> mu[30] 0.8625773 0.10988535 0.55764913 0.8178639 0.8935973 0.9379587 0.9857123
```

Extract area mean estimation

``` r
fit$Est_area
#>         Mean         SD      2.5%       25%       50%       75%     97.5%
#> 1  1.2145310 0.08925347 0.9771087 1.1756451 1.2362475 1.2767032 1.3263293
#> 2  0.9695564 0.11076941 0.7203668 0.9026264 0.9783977 1.0493926 1.1463753
#> 3  1.2794284 0.10808427 0.9980827 1.2240069 1.3001245 1.3577017 1.4223757
#> 4  0.9667121 0.07859691 0.7615088 0.9362999 0.9828508 1.0203278 1.0667598
#> 5  0.9795126 0.16143157 0.6336765 0.8885102 0.9898502 1.1007650 1.2618834
#> 6  1.4801560 0.13070754 1.1868860 1.4072872 1.5058130 1.5761023 1.6666746
#> 7  0.8259752 0.17305504 0.4850795 0.7020690 0.8357205 0.9479483 1.1477195
#> 8  0.8467616 0.14943671 0.5375826 0.7422109 0.8594054 0.9565852 1.0941389
#> 9  0.7162518 0.17877442 0.3512586 0.5979816 0.7225286 0.8371337 1.0451359
#> 10 0.7955220 0.06550355 0.6376557 0.7630145 0.8092711 0.8418919 0.8846159
```

Extract coefficient estimation $\beta$

``` r
fit$coefficient
#>           Mean        SD       2.5%        25%       50%       75%    97.5%
#> b[0] 0.8430376 0.3242408  0.2675138 0.61893948 0.8137375 1.0512616 1.544077
#> b[1] 0.3519589 0.4019291 -0.3677773 0.07786376 0.3534125 0.6169805 1.143463
#> b[2] 1.1169256 0.5110980  0.1081647 0.76438571 1.1033670 1.4624942 2.222093
```

Extract estimation of subarea random effect variance $\sigma^2_u$ and
area random effect variance $\sigma^2_v$

``` r
fit$refVar
#>      b_var    a_var
#> 1 1.017496 1.077574
```

Calculate Relative Standard Error (RSE)

``` r
RSE <- (fit$Est_sub$SD)/(fit$Est_sub$Mean)*100
summary(RSE)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   7.766  10.383  14.656  20.240  26.397  56.759
```
