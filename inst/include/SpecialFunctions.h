#ifndef MONGREL_SPECIALFUNC_H
#define MONGREL_SPECIALFUNC_H

#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

//' Log of Multivarate Gamma Function
//' Gamma_p(a) - https://en.wikipedia.org/wiki/Multivariate_gamma_function
double lmvgamma(double a, int p);
//' Derivative of Log of Multivariate Gamma Function
//' https://en.wikipedia.org/wiki/Multivariate_gamma_function
//' Gamma_p(a)
double lmvgamma_deriv(double a, int p);
//' Efficient, stable calculation of repeated matrix multiplication (for a fixed matrix)
Eigen::MatrixXd power_G(Eigen::MatrixXd G, int it_begin, int it_end);
//' Calculate (marginal) matrix normal mean for DLM
Eigen::MatrixXd dlm_B(Eigen::MatrixXd F, Eigen::MatrixXd G, Eigen::MatrixXd M0, Eigen::VectorXd observations);
//' Calculate (marginal) matrix normal covariance for DLM
Eigen::MatrixXd dlm_U(Eigen::VectorXd F, Eigen::MatrixXd G, Eigen::MatrixXd W, Eigen::MatrixXd C0, Eigen::VectorXd observations);

#endif