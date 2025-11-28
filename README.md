# Bayesian Nonparametric (BNP) Regression for Insurance Claims Data

This repository contains R code for the Dirichlet process mixture model (DPMM) and the Pitman–Yor process mixture model (PYMM) applied to the insurance claims data in our manuscript *"Modeling insurance claims using Bayesian nonparametric regression"* by Mostafa Shams and Kaushik Ghosh.

Each script includes:

- The model’s MCMC sampling algorithm  
- Computation of the posterior predictive distribution  
- Evaluation of clustering performance


## How to Install and Load the Data

The French motor insurance claims data used in this study are third-party data originating from the book Computational Actuarial Science with R (edited by Arthur Charpentier) and are distributed through the R package CASdatasets[^1]. These datasets are publicly and freely available via R by installing the CASdatasets package and loading the datasets `freMTPLfreq` and `freMTPLsev` using the following commands:

```r
install.packages(
	"CASdatasets", 
	repos = "https://dutangc.perso.math.cnrs.fr/RRepository/pub/", 
	type="source"
)

library(CASdatasets) 
data("freMTPLfreq")
data("freMTPLsev")
```


Further installation instructions and documentation for CASdatasets are provided at its official repository:  https://dutangc.github.io/CASdatasets/


Note. The data are not redistributed in this repository; users should obtain them directly via the CASdatasets package. 

**Reference**:  
[^1]: Christophe Dutang and Arthur Charpentier (2025). CASdatasets: Insurance datasets, R package version 1.2-1, DOI: 10.57745/P0KHAG.


