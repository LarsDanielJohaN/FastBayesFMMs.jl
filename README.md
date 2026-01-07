# FastBayesFMMs

FastBayesFMMs is a Julia package that implements the Fast Longitudinal Function on Scalar Regression method developed by Thomas Sun and Dan Kowal which was published at the Journal of Computational and Graphical Statistics in 2024. 

 Consider having $i = 1,2,..., N$ study subjects, which are visted on $M_i$, this quantity may vary throught subjects, where evidence associated to functional observations $Y_{ij}(t)$'s is recorded. This is an example of a longitudinal data context as the functional recordings asociated to a single subject share common information, making methods which relly on independence assumptions innapropiate. We wish to study the efect of covariates $X_{\ell, ij}$ $\ell = 1,2,..., L$ on the $Y_{ij}(t)$, while seeking to incorporate the dependence information. To accomplish this we may use the following Functional Mixed Model: 


 $$Y_{ij}(t) =  \alpha_0(t) + \sum_{\ell = 1}^L X_{\ell, i}\alpha_l(t) + \gamma_{i}(t) + \delta_{ij}(t) $$



[![Build Status](https://github.com/LarsDanielJohaN/FastBayesFMMs.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/LarsDanielJohaN/FastBayesFMMs.jl/actions/workflows/CI.yml?query=branch%3Amaster)
