<img width="160" height="160" src="https://raw.githubusercontent.com/jsilve24/fido/master/inst/fido.png" />

  <!-- badges: start -->
  [![Travis build status](https://travis-ci.org/jsilve24/fido.svg?branch=master)](https://travis-ci.org/jsilve24/fido)
  [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) 
  <!-- badges: end -->

# fido (formerly stray)
Multinomial Logistic-Normal Models (really fast) <br>
*its a little **tar**-ball of joy*

## Citation ##
Silverman, JD, Roche, K, Holmes, ZC, David, LA, and Mukherjee, S. Bayesian Multinomial Logistic Normal Models through Marginally Latent Matrix-T Processes. 2019, arXiv e-prints, arXiv:1903.11695

## License ##
All source code freely availale under [GPL-3 License](https://www.gnu.org/licenses/gpl-3.0.en.html). 

## Installation ##

``` r
devtools::install_github("jsilve24/fido")
```
Or to download the development version

``` r
devtools::install_github("jsilve24/fido", ref="develop")
```

A few notes:

* There are a few installation options that can greatly speed fido up (often by as much as 10-50 fold). For a more detailed description of installation, take a look at [the installation page](https://github.com/jsilve24/fido/wiki/Installation-Details). 
* Vignettes are prebuilt on the [*fido* webpage](https://jsilve24.github.io/fido/). If you 
want vignettes to build locally during package installation you must also pass the `build=TRUE` and `build_opts = c("--no-resave-data", "--no-manual")` options to `install_github`. 


## Bugs/Feature requests ##
Have you checked out [our FAQ](https://github.com/jsilve24/fido/wiki/Frequently-Asked-Questions)? 

I appreciate bug reports and feature requests. Please post to the github issue tracker [here](https://github.com/jsilve24/fido/issues). 


