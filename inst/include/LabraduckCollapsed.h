#ifndef LABRADUCK_MMTC_H
#define LABRADUCK_MMTC_H

#include <RcppNumerical.h>
#include <MatrixAlgebra.h>
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::ArrayXXd;
using Eigen::VectorXd;
using Eigen::Ref;

class LabraduckCollapsed : public Numer::MFuncGrad
{
  private:
    const ArrayXXd Y;
    const double upsilon;
    const MatrixXd B;
    const MatrixXd KInv;
    const MatrixXd U;
    MatrixXd A;
    MatrixXd AInv; // AInv = (gamma*I + W_scale*U)^{-1}
    // computed quantities 
    int D;
    int N;
    bool optimize_gamma_scale;
    bool optimize_W_scale;
    double delta;
    double phi;
    Eigen::ArrayXd m;
    Eigen::RowVectorXd n;
    MatrixXd S;  // I_D-1 + KEAE'
    //Eigen::HouseholderQR<MatrixXd> Sdec;
    Eigen::PartialPivLU<MatrixXd> Sdec;
    Eigen::PartialPivLU<MatrixXd> Adec;
    MatrixXd E;  // eta-B
    ArrayXXd O;  // exp{eta}
    // only needed for gradient and hessian
    MatrixXd rhomat;
    VectorXd rho; 
    MatrixXd C;
    MatrixXd R;
    
    // testing
    bool sylv;
    
    
  public:
    LabraduckCollapsed(const ArrayXXd Y_,          // constructor
                        const double upsilon_,
                        const MatrixXd B_,
                        const MatrixXd KInv_,
                        const MatrixXd U_,
                        bool optimize_gamma_scale_,
                        bool optimize_W_scale_,
                        bool sylv=false) :
    Y(Y_), upsilon(upsilon_), B(B_), KInv(KInv_), U(U_), optimize_gamma_scale(optimize_gamma_scale_), optimize_W_scale(optimize_W_scale_)
    {
      D = Y.rows();           // number of multinomial categories
      N = Y.cols();           // number of samples
      n = Y.colwise().sum();  // total number of counts per sample
      delta = 0.5*(upsilon + N + D - 2.0);
      phi = 0.5*(D-1);
      // the assumption here will be that zero values for gamma_scale or W_scale indicate we should optimize them
      this->sylv = sylv;
    }
    ~LabraduckCollapsed(){}        
    
    // Update with Eta when it comes in as a vector
    void updateWithEtaLL(const Ref<const VectorXd>& etavec, const Ref<const VectorXd>& scale_pars){
      // contents of scale_pars should be (0) gamma (1) W
      const Map<const MatrixXd> eta(etavec.data(), D-1, N);
      E = eta - B;

      double e_gamma = exp(scale_pars(0));
      double e_W = exp(scale_pars(1));
      A = e_gamma*MatrixXd::Identity(N,N) + e_W*U;
      Adec.compute(A);
      AInv = Adec.inverse();
      if (sylv & (N < (D-1))){
        S.noalias() = AInv*E.transpose()*KInv*E;
        S.diagonal() += VectorXd::Ones(N);
      } else {
        S.noalias() = KInv*E*AInv*E.transpose();
        S.diagonal() += VectorXd::Ones(D-1);  
      }
      Sdec.compute(S);
      O = eta.array().exp();
      m = O.colwise().sum();
      m += Eigen::ArrayXd::Ones(N);
    }
    
    // Must be called after updateWithEtaLL 
    void updateWithEtaGH(){
      rhomat = (O.rowwise()/m.transpose()).matrix();
      Map<VectorXd> rhovec(rhomat.data() , rhomat.size());
      rho = rhovec; // probably could be done in one line rather than 2 (above)
      if (sylv & (N < (D-1))){
        C.noalias() = KInv*E;
        R.noalias() = Sdec.solve(AInv); // S^{-1}AInv
      } else {
        C.noalias() = AInv*E.transpose();
        R.noalias() = Sdec.solve(KInv); // S^{-1}KInv    
      }
    }
    
    // Must have called updateWithEtaLL first 
    double calcLogLik(const Ref<const VectorXd>& etavec){
      const Map<const MatrixXd> eta(etavec.data(), D-1, N);
      double ll=0.0;
      // start with multinomial ll
      ll += (Y.topRows(D-1)*eta.array()).sum() - n*m.log().matrix();
      // Now compute collapsed prior ll
      //ll -= delta*Sdec.logAbsDeterminant();
      // Following was adapted from : 
      //   https://gist.github.com/redpony/fc8a0db6b20f7b1a3f23
      double ld = 0.0;
      double c = Sdec.permutationP().determinant();
      VectorXd diagLU = Sdec.matrixLU().diagonal();
      for (unsigned i = 0; i < diagLU.rows(); ++i) {
        const double& lii = diagLU(i);
        if (lii < 0.0) c *= -1;
        ld += log(std::abs(lii));
      }
      ld += log(c);
      ll -= delta*ld;
      // repeat for term -P/2 log |A|
      ld = 0.0;
      c = Adec.permutationP().determinant();
      diagLU = Adec.matrixLU().diagonal();
      for (unsigned i = 0; i < diagLU.rows(); ++i) {
        const double& lii = diagLU(i);
        if (lii < 0.0) c *= -1;
        ld += log(std::abs(lii));
      }
      ld += log(c);
      ll -= phi*ld; // repeated can speed up in future
      return ll;
    }
    
    // Must have called updateWithEtaLL and then updateWithEtaGH first 
    VectorXd calcGrad(const Ref<const VectorXd>& scale_pars){
      // For Multinomial
      MatrixXd g = (Y.topRows(D-1) - (rhomat.array().rowwise()*n.array())).matrix();
      // For MatrixVariate T
      if (sylv & (N < (D-1))){
        g.noalias() += -delta*C*(R+R.transpose());
      } else {
        g.noalias() += -delta*(R + R.transpose())*C.transpose();        
      }
      Map<VectorXd> grad_eta(g.data(), g.size());
      VectorXd grad(grad_eta.size() + 2);
      VectorXd g2(1);
      VectorXd g3(1);
      MatrixXd M;
      double e_gamma = exp(scale_pars(0));
      double e_W = exp(scale_pars(1));
      if (sylv & (N < (D-1))){
        if(optimize_gamma_scale || optimize_W_scale) {
          M = AInv*(E.transpose())*C*R;
          if(optimize_gamma_scale) {
            // gradient for gamma scale
            // Frobenius inner product may work if I can convince myself this whole product is symmetric (TODO)
            g2(0) = delta*M.diagonal().sum();
            g2(0) -= phi*AInv.diagonal().sum();
            g2(0) = e_gamma*g2(0);
          } else {
            g2(0) = 0;
          }
          if(optimize_W_scale) {
            // gradient for W scale
            g3(0) = delta*(M*U).diagonal().sum();
            g3(0) -= phi*(AInv*U).diagonal().sum();
            g3(0) = e_W*g3(0);
          } else {
            g3(0) = 0;
          }
        } else {
          g2(0) = 0;
          g3(0) = 0;
        }
      } else {
        if(optimize_gamma_scale || optimize_W_scale) {
          M = C*R*E*AInv;
          if(optimize_gamma_scale) {
            // gradient for gamma scale
            // Frobenius inner product
            g2(0) = delta*M.diagonal().sum();
            g2(0) -= phi*AInv.diagonal().sum();
            g2(0) = e_gamma*g2(0);
          } else {
            g2(0) = 0;
          }
          if(optimize_W_scale) {
            // gradient for W scale
            g3(0) = delta*(M*U).diagonal().sum();
            g3(0) -= phi*(AInv*U).diagonal().sum();
            g3(0) = e_W*g3(0);
          } else {
            g3(0);
          }
        } else {
          g2(0) = 0;
          g3(0) = 0;
        }
      }
      grad << grad_eta, g2, g3;
      return grad; // not transposing (leaving as vector)
    }
    
    
    // Must have called updateWithEtaLL and then updateWithEtaGH first 
    MatrixXd calcHess(const Ref<const VectorXd>& etavec, const Ref<const VectorXd>& scale_pars){
      bool tmp_sylv = sylv;
      if (sylv & (N < (D-1))){
        this->sylv=false;
        updateWithEtaLL(etavec, scale_pars);
        updateWithEtaGH();
      }
      // for MatrixVariate T
      MatrixXd H(N*(D-1), N*(D-1));
      MatrixXd RCT(D-1, N);
      MatrixXd CR(N, D-1);
      MatrixXd L(N*(D-1), N*(D-1));
      RCT.noalias() = R*C.transpose();
      CR.noalias() = C*R;
      krondense_inplace(L, C*RCT, R.transpose());
      //L.noalias() = krondense(C*RCT, R.transpose());
      krondense_inplace(H, AInv, R+R.transpose());
      //H.noalias() = krondense(AInv, R+R.transpose());
      H.noalias() -= L+L.transpose();
      krondense_inplace(L, RCT, RCT.transpose());
      //L.noalias() = krondense(RCT, RCT.transpose()); // reuse L
      krondense_inplace_add(L, CR.transpose(), CR);
      //L.noalias() += krondense(CR.transpose(), CR); // reuse L
      tveclmult_minus(N, D-1, L, H);
      //H.noalias() -= tveclmult(N, D-1, L);
      H.noalias() = -delta * H;
      
      // For Multinomial
      VectorXd rho_parallel;
      VectorXd n_parallel;
      rho_parallel = rho; 
      n_parallel = n;
      
      #pragma omp parallel shared(rho_parallel, n_parallel)
      {
      MatrixXd W(D-1, D-1);
      //VectorXd rhoseg(D-1);
      #pragma omp for 
      for (int j=0; j<N; j++){
        //rhoseg = rho.segment(j*(D-1), D-1);
        Eigen::Ref<VectorXd> rhoseg = rho_parallel.segment(j*(D-1), D-1);
        W.noalias() = rhoseg*rhoseg.transpose();
        W.diagonal() -= rhoseg;
        H.block(j*(D-1), j*(D-1), D-1, D-1).noalias()  += n_parallel(j)*W;
      }
      }
      // Turn back on sylv option if it was wanted:
      this->sylv = tmp_sylv;
      return H;
    }
    
    // function for use by ADAMOptimizer wrapper (and for RcppNumeric L-BFGS)
    virtual double f_grad(Numer::Constvec& pars, Numer::Refvec grad){
      const Map<const VectorXd> eta(pars.head(N*(D-1)).data(), N*D-1);
      const Map<const VectorXd> scale_pars(pars.tail(2).data(), 2); // may want to scale blocks of W separately in future
      updateWithEtaLL(eta, scale_pars);  // precompute things needed for LogLik
      updateWithEtaGH();                 // precompute things needed for gradient and hessian
      grad = -calcGrad(scale_pars);      // negative because wraper minimizes
      return -calcLogLik(eta);           // negative because wraper minimizes
    }
};


#endif