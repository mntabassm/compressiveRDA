# R-package | _compressiveRDA_  
* __Version__: 1.1 
* __GitHub Link__: https://github.com/mntabassm/compressiveRDA

* __Title__: Compressive regularized discriminant analysis for simultaneous feature selection and classification of high-dimensional data, with applications to genomic studies. 	 
* __Short Title__: _Compressive Regularized Discriminant Analysis (CRDA)_
* __Authors__: Muhammad Naveed Tabassum and Esa Ollila
* __Maintainer__: Muhammad Naveed Tabassum
* __Language__: R
* __Date__: 26.09.2018
* __Date (Last update)__: 15.10.2019

## Introduction 

The <tt>compressiveRDA</tt> pacakge implements the CRDA approach whose goal is to address three facets of high-dimensional classification: namely, accuracy, computational complexity, and interpretability. The currently available competitors of CRDA method present a weak spot for at least one of the aforementioned criteria of an HD classifier.  

## Installation

The <tt>compressiveRDA</tt> pacakge can be installed from GitHub, using the <tt>devtools</tt> pacakge as:

``` r
devtools::install_github("mntabassm/compressiveRDA")
library(compressiveRDA)
```

NOTE: If there is some problem coming then, do as:
``` r
devtools::install_github("mntabassm/compressiveRDA", force = TRUE)
library(compressiveRDA)
```

<!-- or by the <tt>githubinstall</tt> package that provides a function _githubinstall_. It does not need developer’s name and work as:

``` r
install.packages('githubinstall', dependencies = TRUE)
library(githubinstall)
githubinstall::githubinstall("compressiveRDA")
library(compressiveRDA)
```
-->
## Example

As an example, just run the function 'crda.demo()' that does the classification for one split of a real genomic dataset, Khan'2001.

## basic example code

> * crda.demo() : It does classification using a uniform prior.

> * crda.demo(prior = 'estimated') : It does classification using a empirically estimated prior.
