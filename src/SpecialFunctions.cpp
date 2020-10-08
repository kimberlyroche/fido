#include <Rcpp.h>
#include <RcppEigen.h>
#include "MatDist.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::MatrixXcd;
using Eigen::VectorXd;
using Eigen::VectorXcd;
using Eigen::ArrayXXd;
using Eigen::Map;

//' Log of Multivarate Gamma Function - Gamma_p(a)
//' @param a defined by Gamma_p(a)
//' @param p defined by Gamma_p(a)
//' @references https://en.wikipedia.org/wiki/Multivariate_gamma_function
// [[Rcpp::export]]
double lmvgamma(double a, int p){
  static const double pi = log(3.14159265); 
  double s=0;
  double x = pi*(p*(p-1.0))/2.0;
  for (int i=1; i<=p; i++){
    s += lgamma(a+(1.0-i)/2);
  }
  return(x+s);
}

//' Derivative of Log of Multivariate Gamma Function - Gamma_p(a)
//' @param a defined by Gamma_p(a)
//' @param p defined by Gamma_p(a)
//' @references https://en.wikipedia.org/wiki/Multivariate_gamma_function
// [[Rcpp::export]]
double lmvgamma_deriv(double a, int p){
  double s=0;
  for (int i=1; i<=p; i++){
    s += R::digamma(a + 0.5*(1-i));
  }
  return s*lmvgamma(a,p);
}


// explicit references to standard namespace here are motivated by the observation that
// these calls don't properly resolve on some Windows installs
// [[Rcpp::export]]
Eigen::MatrixXd power_G(Eigen::MatrixXd G, 
                        int it_begin,
                        int it_end){
  if(it_begin == it_end) { return(G); }
  int p = G.rows();
  int power_it = std::abs(it_end - it_begin) + 1;
  // powering must explicitly be carried out on complex components (if they exist) using
  // complex exponents; doesn't look like these a vector wrapper for this; perform
  // element-wise powering by hand
  std::complex<double> power_complex(power_it, 0);
  Eigen::ComplexEigenSolver<MatrixXcd> Gdec(G);
  VectorXcd powered_eigenvalues = Gdec.eigenvalues();
  for(int j=0; j<p; j++) {
    powered_eigenvalues(j) = std::pow(powered_eigenvalues(j), power_complex);
  }
  MatrixXcd Lambda = MatrixXcd::Zero(p,p);
  // use of .diagonal() in assignment doesn't work here for type reasons
  // troubleshoot this
  for(int j=0; j<p; j++) {
    Lambda(j,j) = powered_eigenvalues(j,1);
  }
  MatrixXcd reconstituted = Gdec.eigenvectors()*Lambda*(Gdec.eigenvectors().inverse());
  // in R the imaginary part was occasionally vanishingly non-zero
  // doesn't seem to happen here but chop it off just in case
  return reconstituted.real();
}

// [[Rcpp::export]]
Eigen::MatrixXd dlm_B(Eigen::MatrixXd F, Eigen::MatrixXd G, Eigen::MatrixXd M0, Eigen::VectorXd observations) {
  int D = M0.cols()+1;
  int N = observations.size();
  int T = observations.maxCoeff();
  MatrixXd B(D-1, N);
  MatrixXd alpha(1, D);
  MatrixXd Ft = F.transpose();
  int t;
  for(int t_incr=0; t_incr<observations.size(); t_incr++) {
    // column-wise
    t = observations(t_incr);
    alpha = Ft*power_G(G, t, 1)*M0;
    B.block(0,t_incr,D-1,1) = alpha;
  }
  return(B);
}

// build U matrix (unscaled covariance over samples) assuming time-invariant parameters F, G, W
// [[Rcpp::export]]
Eigen::MatrixXd dlm_U(Eigen::VectorXd F, Eigen::MatrixXd G, Eigen::MatrixXd W, Eigen::MatrixXd C0, Eigen::VectorXd observations) {
  // check T >= 1
  int N = observations.size();
  int T = observations.maxCoeff();
  MatrixXd res = MatrixXd::Zero(N,N);
  MatrixXd Ft = F.transpose();
  MatrixXd Gt = G.transpose();
  int system_dim = G.rows();
  int t;
  int tk;
  for(int t1_incr=0; t1_incr<observations.size(); t1_incr++) { // cov rows
    for(int t2_incr=0; t2_incr<observations.size(); t2_incr++) { // cov cols
      // check this observation exists, else skip this row/col in covariance matrix
      if(t1_incr == t2_incr) {
        t = observations(t1_incr);
        // diagonal
        res(t1_incr, t1_incr) += (Ft*W*F)(0,0); // 1x1 vector product
        for(int ell=t; ell>=2; ell--) {
          res(t1_incr, t1_incr) += (Ft*power_G(G, t, ell)*W*power_G(Gt, ell, t)*F)(0,0); // 1x1 vector
        }
        res(t1_incr, t1_incr) += (Ft*power_G(G, t, 1)*C0*power_G(Gt, 1, t)*F)(0,0); // 1x1 vector
      } else {
        // off-diagonal
        if(t1_incr > t2_incr) {
          // symmetric: if you're in the lower triangular portion, just copy from the (already computed) upper triangular
          res(t1_incr, t2_incr) = res(t2_incr, t1_incr);
        } else {
          t = observations(t2_incr); // col
          tk = observations(t1_incr); // row
          res(t1_incr, t2_incr) += (Ft*power_G(G, t, tk+1)*W*F)(0,0);
          MatrixXd temp = MatrixXd::Zero(system_dim, system_dim);
          for(int ell=tk; ell>=2; ell--) {
            temp += power_G(G, t, ell)*W*power_G(Gt, ell, tk);
          }
          res(t1_incr, t2_incr) += (Ft*temp*F)(0,0);
          res(t1_incr, t2_incr) += (Ft*power_G(G, t, 1)*C0*power_G(Gt, 1, tk)*F)(0,0);
        }
      }
    }
  }
  return res;
}
