# FastBayesFMMs

FastBayesFMMs is a Julia package that implements the Fast Longitudinal Function on Scalar Regression method developed by Thomas Sun and Dan Kowal which was published at the Journal of Computational and Graphical Statistics in 2024. 

 Consider having $i = 1,2,..., N$ study subjects, which are visted on $M_i$ occasions, this quantity may vary throught subjects, where evidence associated to functional observations $Y_{ij}(t)$'s is recorded. This is an example of a longitudinal data context as the functional recordings asociated to a single subject share common information, making methods which relly on independence assumptions innapropiate. We wish to study the efect of covariates $X_{\ell, ij}$ $\ell = 1,2,..., L$ on the $Y_{ij}(t)$, while seeking to incorporate the dependence information. To accomplish this we may use the following Functional Mixed Model: 


 $$Y_{ij}(t) = \alpha_0(t) + \sum_{\ell = 1}^L X_{\ell, ij}\alpha_\ell(t) + \gamma_{i}(t) + \omega_{ij}(t) + \varepsilon_{ij}(t)$$

In this model $\alpha_{0}(t)$ is interpreted pointwise as the avarage value that $Y_{ij}(t)$ would take if the values $X_{ij, \ell}$ were equal to zero. For $\alpha_\ell(t)$'s and fixing $t^\prime \in \mathcal{T}$, $\alpha_\ell(t^\prime)$ is interpreted as the change in $Y_{ij}(t^\prime)$ for a unit increase in $X_{ij, \ell}$ at time $t^\prime$ if everyting else remains constant. For a fixed subject $i^\prime$, the term $\gamma_{i^\prime}(t)$ *absorbs* the dependence that $i^\prime$ carries across its $M_i$ observations. Because of this, $\gamma_i(t)$ is called a *subject level random effect* and is used to explain the variation between the $N$ subjects. These are assumed to be independent among subjects. $\omega_{ij}(t)$ is called a *visit specific random effects* and, as its name suggests, it has the purpose of measuring the deviation of $Y_{ij}(t)$ which is specific to the $i$-th subjects $j$-th visit. These are independent among subjects and visits.  Because of this, the $\omega_{ij}(t)$'s are used to measure within-subject variation. 


[![Build Status](https://github.com/LarsDanielJohaN/FastBayesFMMs.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/LarsDanielJohaN/FastBayesFMMs.jl/actions/workflows/CI.yml?query=branch%3Amaster)
