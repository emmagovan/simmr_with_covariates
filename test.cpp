#include <RcppArmadillo.h>
#include <RcppDist.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

static double const log2pi = std::log(2.0 * M_PI);


void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}


// [[Rcpp::export]]
arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = true) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}


// [[Rcpp::export]]
NumericMatrix crossprod(NumericMatrix X){
  NumericMatrix ans(X.nrow(), X.ncol());
  
  for(int i = 0; i<X.ncol(); i++){
    for(int j=0; j<X.ncol(); j++){
      for(int n =0; n<X.nrow(); n++){
        
        ans(i,j) += X(n,i) * X(n,j);
      }
    }}
  return(ans);
}

// [[Rcpp::export]]
NumericMatrix matmult(NumericMatrix x, NumericMatrix y) {
  NumericMatrix ans(x.nrow(), y.ncol());
  
  for(int i = 0; i<x.nrow(); i++){
    for(int j=0; j<y.ncol(); j++){
      for(int k =0; k<y.nrow(); k++){
        
        ans(i,j) += x(i,k) * y(k,j);
      }
    }}
  
  return ans;
}


// [[Rcpp::export]]
NumericMatrix rMVNormCpp(const double n,
                         const arma::vec mu,
                         const NumericMatrix U) {
  
  
  // Dimension of MVN
  int p = mu.size();
  
  // Simulate iid standard normals
  arma::mat Z(p, n);
  Z.imbue(norm_rand);
  
  // Now backsolve and add back on the means
  arma::mat X = solve(as<arma::mat>(U), Z);
  for ( int i = 0; i < n; ++i ) {
    X.col(i) += mu;
  }
  
  return Rcpp::wrap(X.t());
}



// [[Rcpp::export]]
NumericMatrix solvearma(const NumericMatrix X) {

  arma::mat b = arma::eye(X.nrow(), X.ncol());


  // Now backsolve and add back on the means
  arma::mat ans = solve(as<arma::mat>(X), b);


  return Rcpp::wrap(ans.t());
}


//[[Rcpp::export]]
NumericMatrix sim_thetacpp(int S, NumericVector lambda, int n_sources,
                           int n_tracers, int n_cov){
  NumericMatrix theta(S, (n_cov*n_sources + n_tracers));
  NumericMatrix mean_beta((n_cov), n_sources);
  int mat_size = n_sources * (n_sources+1) /2;

  for(int i=0; i<n_cov; i++){
    for(int k=0; k<n_sources;k++){
      mean_beta(i,k) = lambda(i * mat_size + i * n_sources + k);
    }
  }

  NumericMatrix sig_beta(n_cov, mat_size);

  for(int m = 0; m<mat_size; m++){
    for(int i =0; i<n_cov; i++){
    sig_beta(i,m) = lambda(i* mat_size + (i+1) * n_sources + m);

  }
  }

  NumericVector count(n_cov);

  for(int i =0; i<n_cov; i++){
    count(i) = 0;
  }

  arma::cube chol_prec(n_sources, n_sources, n_cov);

  for(int j = 0; j< n_sources; j++){
    for(int i = 0; i<n_sources; i++){
      for(int m = 0; m<n_cov; m++){
        if (i <= j){
          count(m) +=1;
          chol_prec(i,j,m) = sig_beta(m, count(m)-1);


        }

        else{
          chol_prec(i,j,m) = 0;
        }
      }
    }
  }




arma::mat theta_arma(S, (n_cov*n_sources + n_tracers)); //Want to go from chol_prec array to
  // A matrix of thetas generated using rMVNormCpp

  for(int i=0; i<n_cov; i++){
  theta_arma.submat(0, (i)*n_sources, S-1, (i+1)*n_sources - 1) = as<arma::mat>(rMVNormCpp(S, mean_beta(i,_), Rcpp::wrap(chol_prec.slice(i))));
  }

theta = Rcpp::wrap(theta_arma);



  for(int i = 0; i<n_tracers; i++){
    theta(_,i+n_sources*n_cov) = (Rcpp::rgamma(S,  lambda(n_cov * mat_size + n_cov *n_sources +i),
          1/lambda(n_cov * mat_size + n_cov *n_sources +i + n_tracers)));
  }


  return theta;
}


//[[Rcpp::export]]
NumericMatrix hfn(NumericVector theta, int n_sources, int n, int n_cov, NumericMatrix x_scaled){
  NumericMatrix p(n, n_sources);
  NumericMatrix exptheta(n_sources, n_sources);
  NumericMatrix f(n, n_sources);
  NumericMatrix beta(n_cov, n_sources);
  

    for(int i = 0; i<n_cov; i++){
      for(int j=0; j<n_sources; j++){
      beta(i,j) = theta((i)*n_sources +j);
    }
    }
    
    f = matmult(x_scaled, beta);
    
    
    NumericMatrix expf(n, n_sources); 
    
    for(int i =0; i<n; i++){
      for(int j = 0; j<n_sources; j++){
        expf(i,j) = exp(f(i,j));
      }
    }
    
    NumericVector sumexpf(n);
    
    for(int i = 0; i<n; i++){
      for(int j=0; j<n_sources; j++){
        sumexpf(i) +=expf(i,j);
      }
    }
    
    for(int i=0; i<n; i++){
      for(int j =0; j<n_sources; j++){
        p(i,j) = expf(i,j)/sumexpf(i);
      }
    }


  
  return p;
  
}


// Next step is to do the actual h function itself - just a case of changing the maths and looping over isotopes I think
// This will also make it so it can use any number of isotopes, not just 2



//[[Rcpp::export]]
double hcpp(int n_sources, int n_isotopes, int n_covariates,
            double beta_prior,
            NumericMatrix x_scaled,
            NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
            NumericMatrix correctionmeans,
            NumericMatrix corrsds, NumericMatrix sourcesds,
            NumericVector theta, NumericMatrix y){



  Rcpp::NumericMatrix beta(n_covariates, n_sources);

  for(int i = 0; i<n_covariates; i++){
    for(int j=0; j<n_sources; j++){
      beta(i,j) = theta((i)*n_sources +j);
    }
  }

  int n = y.rows();
  
  NumericMatrix p(n, n_sources);

  p = hfn(theta, n_sources, n, n_covariates, x_scaled);



  // Setting prior values for hyper parameters
  NumericMatrix prior_means(n_covariates, n_sources);
  NumericMatrix prior_sd(n_covariates, n_sources);
  NumericVector c_0(n_isotopes);
  NumericVector d_0(n_isotopes);

  // Setting up prior values
  for(int i=0; i<n_covariates; i++){
    for(int j=0; j<n_sources; j++){

    prior_means(i,j) = 0;
    prior_sd(i,j) = 1;
  }
  }

  for (int i = 0; i<n_isotopes; i++){
    c_0(i) = 0.001;
    d_0(i) = beta_prior;
  }
  
  double hold = 0;
  
  NumericMatrix mutop(n*n_isotopes, n_sources);
  NumericMatrix mubtm(n*n_isotopes, n_sources);
  
  for(int j=0; j<n_isotopes; j++){
  for(int i=0; i<n; i++){
      mutop(i+j*n,_) = p(i,_)*concentrationmeans(_,j) * (sourcemeans(_,j) + correctionmeans(_,j));
      mubtm(i+j*n, _) = p(i,_)*concentrationmeans(_,j);
    }
  }
  
NumericVector musumt(n*n_isotopes);
  NumericVector musumb(n*n_isotopes);
  
  for(int j=0; j<n_sources; j++){
  for(int i=0; i<(n*n_isotopes); i++){
      musumt(i) += mutop(i,j);
      musumb(i) += mubtm(i,j);
    }
  }
  
  NumericVector mutotal(n*n_isotopes);
  
  for(int i=0; i<(n*n_isotopes); i++){
    mutotal(i) = musumt(i)/musumb(i);
  }
  
  // Now do the same for sigma
  
  NumericMatrix sigtop(n*n_isotopes, n_sources);
  NumericMatrix sigbtm(n*n_isotopes, n_sources);
  
  for(int j=0; j<n_isotopes; j++){
    for(int i=0; i<n; i++){
      sigtop(i+j*n,_) = pow(p(i,_), 2)*pow(concentrationmeans(_,j),2) * (pow(sourcemeans(_,j),2) + pow(correctionmeans(_,j),2));
      sigbtm(i+j*n, _) = pow(p(i,_)*concentrationmeans(_,j),2);
    }
  }
  
  NumericVector sigsumt(n*n_isotopes);
  NumericVector sigsumb(n*n_isotopes);
  
  for(int j=0; j<n_sources; j++){
    for(int i=0; i<(n*n_isotopes); i++){
      sigsumt(i) += sigtop(i,j);
      sigsumb(i) += sigbtm(i,j);
    }
  }
  
  NumericVector sigtotal(n*n_isotopes);
  
  for(int i=0; i<(n*n_isotopes); i++){
    sigtotal(i) = pow(sigsumt(i)/sigsumb(i), 0.5);
  }
  
  
  for(int i=0; i<(n); i++){
    for(int j=0; j<n_isotopes; j++){
    hold = hold -n *log(sigtotal(i+j*n)) - 0.5 *n *log(2 * M_PI) - 
      0.5 * pow((y(i,j)-mutotal(i+j*n)),2) * 1/pow(sigtotal(i+j*n),2);
  }

  }


  double betanorm = 0;



  for(int i = 0; i<n_covariates; i++){
    for(int j=0; j<n_sources; j++){
    betanorm +=  - n_sources * log(prior_sd(i,j)) - 0.5 * n_sources* log(2 * M_PI) -
       0.5 * (pow((beta(i,j) - prior_means(i,j)), 2)
      *1/pow(prior_sd(i,j), 2));
  }
  }
  double gammaprior = 0;
  for (int i=0; i <(n_isotopes); i++){
    gammaprior += c_0(i) * log(d_0(i)) - log(tgamma(c_0(i))) +
      (c_0(i) - 1) * theta((i+n_sources*n_covariates))-
      d_0(i) * theta((i+n_sources*n_covariates));

  }

  double totx = hold + gammaprior + betanorm;
  
   return (totx);

}



//[[Rcpp::export]]
double log_q_cpp(NumericVector theta, NumericVector lambda, 
                 int n_sources, int n_tracers, int S, int n_covariates){
  
  NumericMatrix mean_beta((n_covariates), n_sources);
  int mat_size = n_sources * (n_sources+1) /2;
  
  for(int i=0; i<n_covariates; i++){
    for(int k=0; k<n_sources;k++){
      mean_beta(i,k) = lambda(i * mat_size + i * n_sources + k);
    }
  }
  
  NumericVector count(n_covariates);
  
  for(int i =0; i<n_covariates; i++){
    count(i) = 0;
  }
  NumericMatrix sig_beta(n_covariates, mat_size);
  
  for(int m = 0; m<mat_size; m++){
    for(int i =0; i<n_covariates; i++){
      sig_beta(i,m) = lambda(i* mat_size + (i+1) * n_sources + m);
      
    }
  }
  
  arma::cube chol_prec(n_sources, n_sources, n_covariates);
  
  for(int j = 0; j< n_sources; j++){
    for(int i = 0; i<n_sources; i++){
      for(int m = 0; m<n_covariates; m++){
        if (i <= j){
          count(m) +=1;
          chol_prec(i,j,m) = sig_beta(m, count(m)-1);
          
          
        }
        
        else{
          chol_prec(i,j,m) = 0;
        }
      }
    }
  }
  
  Rcpp::NumericMatrix beta(n_covariates, n_sources);
  
  for(int i = 0; i<n_covariates; i++){
    for(int j=0; j<n_sources; j++){
      beta(i,j) = theta((i)*n_sources +j);
    }
  }
  
  NumericVector sigma(n_tracers);
  
  for(int i=0; i<n_tracers; i++){
    sigma(i) = theta(n_covariates*n_sources+i);
  }
  
  NumericMatrix pmat(n_covariates, n_sources);
  
  NumericMatrix y(n_covariates, beta.ncol());
  
  for(int i=0; i<n_covariates; i++){
   y(i,_) = beta(i,_) - mean_beta(i,_);
  }
  

  

  


 double thetanorm = 0;
 for(int i=0; i<n_covariates; i++){
   NumericMatrix prec(n_sources, n_sources);
   prec = crossprod(Rcpp::wrap(chol_prec.slice(i)));
   NumericMatrix solve_prec(n_sources, n_sources);
   solve_prec = solvearma(prec);
  thetanorm += *REAL(Rcpp::wrap(dmvnrm_arma_fast(as<arma::mat>(y), beta(i,_), as<arma::mat>(solve_prec))));
 }

    
    
    
    double gamman = 0;
    for (int i=0; i <(n_tracers); i++){
      gamman += lambda(n_covariates * mat_size + n_covariates *n_sources +i) * log(lambda(n_covariates * mat_size + n_covariates *n_sources +i + n_tracers))  -
      log(tgamma(lambda(n_covariates * mat_size + n_covariates *n_sources +i)))  +
        lambda((n_covariates * mat_size + n_covariates *n_sources +i) - 1) * log(theta((i+n_sources*n_covariates))) - 
      lambda(n_covariates * mat_size + n_covariates *n_sources +i + n_tracers) * theta((i+n_sources*n_covariates));
    }
    
    
    
    
    double x = thetanorm + gamman;
    
    return x;
    
}






