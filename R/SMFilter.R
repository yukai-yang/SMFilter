################################################################################
## package name: SMFilter
## author: Yukai Yang
## Statistiska Inst., Uppsala Universitet
## Nov 2018
#################################################################################


#' SMFilter: a package implementing the filtering algorithms for the state-space models on the Stiefel manifold.
#'
#' The package implements the filtering algorithms for the state-space models on the Stiefel manifold.
#' It also implements sampling algorithms for uniform, vector Langevin-Bingham and matrix Langevin-Bingham distributions on the Stiefel manifold.
#'
#' Two types of the state-space models on the Stiefel manifold are considered.
#'
#' The type one model on Stiefel manifold takes the form:
#' \deqn{\boldsymbol{y}_t \quad = \quad \boldsymbol{\alpha}_t \boldsymbol{\beta} ' \boldsymbol{x}_t + \boldsymbol{B} \boldsymbol{z}_t + \boldsymbol{\varepsilon}_t}
#' \deqn{\boldsymbol{\alpha}_{t+1} | \boldsymbol{\alpha}_{t} \quad \sim \quad ML (p, r, \boldsymbol{\alpha}_{t} \boldsymbol{D})}
#' where \eqn{\boldsymbol{y}_t} is a \eqn{p}-vector of the dependent variable,
#' \eqn{\boldsymbol{x}_t} and \eqn{\boldsymbol{z}_t} are explanatory variables wit dimension \eqn{q_1} and \eqn{q_2},
#' \eqn{\boldsymbol{x}_t} and \eqn{\boldsymbol{z}_t} have no overlap,
#' matrix \eqn{\boldsymbol{B}} is the coefficients for \eqn{\boldsymbol{z}_t},
#' \eqn{\boldsymbol{\varepsilon}_t} is the error vector.
#'
#' The matrices \eqn{\boldsymbol{\alpha}_t} and \eqn{\boldsymbol{\beta}} have dimensions \eqn{p \times r} and \eqn{q_1 \times r}, respectively.
#' Note that \eqn{r} is strictly smaller than both \eqn{p} and \eqn{q_1}.
#' \eqn{\boldsymbol{\alpha}_t} and \eqn{\boldsymbol{\beta}} are both non-singular matrices.
#' \eqn{\boldsymbol{\alpha}_t} is time-varying while \eqn{\boldsymbol{\beta}} is time-invariant.
#'
#' Furthermore, \eqn{\boldsymbol{\alpha}_t} fulfills the condition \eqn{\boldsymbol{\alpha}_t' \boldsymbol{\alpha}_t = \boldsymbol{I}_r},
#' and therefor it evolves on the Stiefel manifold.
#'
#' \eqn{ML (p, r, \boldsymbol{\alpha}_{t} \boldsymbol{D})} denotes the Matrix Langevin distribution or matrix von Mises-Fisher distribution on the Stiefel manifold.
#' Its density function takes the form
#' \deqn{f(\boldsymbol{\alpha_{t+1}} ) = \frac{ \mathrm{etr} \left\{ \boldsymbol{D} \boldsymbol{\alpha}_{t}' \boldsymbol{\alpha_{t+1}} \right\} }{ _{0}F_1 (\frac{p}{2}; \frac{1}{4}\boldsymbol{D}^2 ) }}
#' where \eqn{\mathrm{etr}} denotes \eqn{\mathrm{exp}(\mathrm{tr}())},
#' and \eqn{_{0}F_1 (\frac{p}{2}; \frac{1}{4}\boldsymbol{D}^2 )} is the (0,1)-type hypergeometric function for matrix.
#'
#' The type two model on Stiefel manifold takes the form:
#' \deqn{\boldsymbol{y}_t \quad = \quad \boldsymbol{\alpha} \boldsymbol{\beta}_t ' \boldsymbol{x}_t + \boldsymbol{B}' \boldsymbol{z}_t + \boldsymbol{\varepsilon}_t}
#' \deqn{\boldsymbol{\beta}_{t+1} | \boldsymbol{\beta}_{t} \quad \sim \quad ML (q_1, r, \boldsymbol{\beta}_{t} \boldsymbol{D})}
#' where \eqn{\boldsymbol{y}_t} is a \eqn{p}-vector of the dependent variable,
#' \eqn{\boldsymbol{x}_t} and \eqn{\boldsymbol{z}_t} are explanatory variables wit dimension \eqn{q_1} and \eqn{q_2},
#' \eqn{\boldsymbol{x}_t} and \eqn{\boldsymbol{z}_t} have no overlap,
#' matrix \eqn{\boldsymbol{B}} is the coefficients for \eqn{\boldsymbol{z}_t},
#' \eqn{\boldsymbol{\varepsilon}_t} is the error vector.
#'
#' The matrices \eqn{\boldsymbol{\alpha}} and \eqn{\boldsymbol{\beta}_t} have dimensions \eqn{p \times r} and \eqn{q_1 \times r}, respectively.
#' Note that \eqn{r} is strictly smaller than both \eqn{p} and \eqn{q_1}.
#' \eqn{\boldsymbol{\alpha}} and \eqn{\boldsymbol{\beta}_t} are both non-singular matrices.
#' \eqn{\boldsymbol{\beta}_t} is time-varying while \eqn{\boldsymbol{\alpha}} is time-invariant.
#'
#' Furthermore, \eqn{\boldsymbol{\beta}_t} fulfills the condition \eqn{\boldsymbol{\beta}_t' \boldsymbol{\beta}_t = \boldsymbol{I}_r},
#' and therefor it evolves on the Stiefel manifold.
#'
#' \eqn{ML (p, r, \boldsymbol{\beta}_t \boldsymbol{D})} denotes the Matrix Langevin distribution or matrix von Mises-Fisher distribution on the Stiefel manifold.
#' Its density function takes the form
#' \deqn{f(\boldsymbol{\beta_{t+1}} ) = \frac{ \mathrm{etr} \left\{ \boldsymbol{D} \boldsymbol{\beta}_{t}' \boldsymbol{\beta_{t+1}} \right\} }{ _{0}F_1 (\frac{p}{2}; \frac{1}{4}\boldsymbol{D}^2 ) }}
#' where \eqn{\mathrm{etr}} denotes \eqn{\mathrm{exp}(\mathrm{tr}())},
#' and \eqn{_{0}F_1 (\frac{p}{2}; \frac{1}{4}\boldsymbol{D}^2 )} is the (0,1)-type hypergeometric function for matrix.
#'
#'
#' @section Author and Maintainer:
#' Yukai Yang
#'
#' Department of Statistics, Uppsala University
#'
#' \email{yukai.yang@@statistik.uu.se}
#'
#' @section References:
#' Yang, Yukai and Bauwens, Luc. (2018) "\href{https://www.mdpi.com/2225-1146/6/4/48}{State-Space Models on the Stiefel Manifold with a New Approach to Nonlinear Filtering}", Econometrics, 6(4), 48.
#'
#' @section Simulation:
#' \code{\link{SimModel1}} simulate from the type one state-space model on the Stiefel manifold.
#'
#' \code{\link{SimModel2}} simulate from the type two state-space model on the Stiefel manifold.
#'
#' @section Filtering:
#' \code{\link{FilterModel1}} filtering algorithm for the type one model.
#'
#' \code{\link{FilterModel2}} filtering algorithm for the type two model.
#'
#' @section Sampling:
#' \code{\link{runif_sm}} sample from the uniform distribution on the Stiefel manifold.
#'
#' \code{\link{rvlb_sm}} sample from the vector Langevin-Bingham distribution on the Stiefel manifold.
#'
#' \code{\link{rmLB_sm}} sample from the matrix Langevin-Bingham distribution on the Stiefel manifold.
#'
#' @section Other Functions:
#' \code{\link{version}} shows the version number and some information of the package.
#'
#' @docType package
#' @name SMFilter
NULL

#' @importFrom stats optim rbeta rnorm runif
NULL
