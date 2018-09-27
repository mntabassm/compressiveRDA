# R-package | _compressiveRDA_  
* __Version__: 1.0.0 
* __GitHub Link__: https://github.com/mntabassm/compressiveRDA

* __Title__: Compressive regularized discriminant analysis for simultaneous feature selection and classification of high-dimensional data, with applications to genomic studies. 	 
* __Short Title__: _Compressive Regularized Discriminant Analysis (CRDA)_
* __Authors__: Muhammad Naveed Tabassum and Esa Ollila
* __Maintainer__: Muhammad Naveed Tabassum
* __Language__: R
* __Date__: 26.09.2018
* __Date (Last update)__: 26.09.2018

## Introduction 

The compressiveRDA implements the CRDA approach whose goal is to address three facets of high-dimensional classification: namely, accuracy, computational complexity, and interpretability. The currently available competitors of CRDA method present a weak spot for at least one of the aforementioned criteria of an HD classifier.  

## Installation

The <tt> 'compressiveRDA' </tt> pacakge can be installed from GitHub, using the <tt> 'devtools' </tt> pacakge as:

``` r
devtools::install_github("mntabassm/compressiveRDA")
library(compressiveRDA)
```

## Example

As an example, just run the function 'crda.classify()' that does the classification for setup 3 of the paper.

## basic example code

> * crda.classify() : It does classification using defualt settings (e.g., L=10 training and test splits)

> * crda.classify(L=1) : It does classification for one split of data.
