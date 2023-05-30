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
            NumericMatrix concentrationmeans, NumericMatrix sourcemeans,
            NumericMatrix correctionmeans,
            NumericMatrix corrsds, NumericMatrix sourcesds, NumericVector theta, NumericMatrix y ){
  
  double x =0;
  
  NumericVector p(n_sources);
  
  p = hfn(theta, n_sources);
  
  double ly = y.rows();
  
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
  
  
    
    
    // This is to get dnorm(y[,1], sum(p*q) etc)
    double mutop1 = 0;
    double mubtm1 = 0;
    double mu1 = 0;
    double sigmasq1 = 0;
    double sigmatopsq1 = 0;
    double sigmabtmsq1 = 0;
    
    
    
    // Calculate numerator and denominator of mu
    for(int i=0; i<n_sources; i++){
      for(int j=0; j<n_isotopes; j++){
        mutop1 +=  p(i)*concentrationmeans(i,j) * (sourcemeans(i,0) + correctionmeans(i,0));
        mubtm1 += p(i) * concentrationmeans(i,j);
      }
    }
    
    // Same for sigma
    for(int i=0; i<n_sources; i++){
      for(int j =0; j<n_isotopes; j++){
        sigmatopsq1 += pow(p(i),2) * pow(concentrationmeans(i,j),2) * (pow(sourcesds(i,0),2) +
          pow(corrsds(i,0),2));
        sigmabtmsq1 += pow(p(i),2) * pow(concentrationmeans(i,j),2);
      }
    }
    
    //Calculate mu and sd
    mu1 = mutop1/mubtm1;
    sigmasq1 = sigmatopsq1/sigmabtmsq1;
    double sigma1 = pow(sigmasq1 + 1/theta((n_sources)), 0.5);
    
    
    // This is to get dnorm(y[,2], sum(p*q) etc)
    double mutop2 = 0;
    double mubtm2 = 0;
    double mu2 = 0;
    double sigmasq2 = 0;
    double sigmatopsq2 = 0;
    double sigmabtmsq2 = 0;
    for(int i=0; i<n_sources; i++){
      for(int j =0; j<n_isotopes; j++){
        mutop2 += p(i) * concentrationmeans(i,j) * (sourcemeans(i,1) + correctionmeans(i,1));
        mubtm2 += p(i) * concentrationmeans(i,j);
      }
    }
    
    
    for(int i=0; i<n_sources; i++){
      for(int j=0; j<n_isotopes; j++){
        sigmatopsq2 += pow(p(i),2) * pow(concentrationmeans(i,j),2) * (pow(sourcesds(i,1),2) +
          pow(corrsds(i,1),2));
        sigmabtmsq2 += pow(p(i),2) * pow(concentrationmeans(i,j),2);
      }
    }
    
    mu2 = mutop2/mubtm2;
    sigmasq2 = sigmatopsq2/sigmabtmsq2;
    
    double sigma2 = pow(sigmasq2 + 1/theta((1+n_sources)), 0.5);
    
    double yminusmu1 = 0;
    double yminusmu2 = 0;
    
    for(int i = 0; i<ly; i++){
      yminusmu1 += pow((y(i,0) - mu1),2);
      yminusmu2 +=  pow((y(i,1) - mu2),2);
    }
    
    
    
    // This is log(dnorm(y, p*q, p^2*q^2 etc) for y1 and y2
    
    x = - ly * log(sigma1) - 0.5 * ly * log(2 * M_PI)
      - 0.5 * yminusmu1 * 1/(pow(sigma1,2))
      - ly * log(sigma2) - 0.5 * ly * log(2 * M_PI)
      - 0.5 * yminusmu2 * 1/(pow(sigma2,2));
      
      
  
  
  double thetanorm = 0;
  
  
  
  
  for(int i = 0; i<n_sources; i++){
    thetanorm +=  - n_sources * log(prior_sd(i)) - 0.5 * log(2 * M_PI) - (pow((theta(i) - prior_means(i)), 2)
                                                                            * 1/(2 * pow(prior_sd(i), 2)));
  }
  
  double gammaprior = 0;
  for (int i=0; i <(n_isotopes); i++){
    gammaprior += c_0(i) * log(d_0(i)) - log(tgamma(c_0(i))) +(c_0(i) - 1) * theta((i+n_sources)) -
      d_0(i) * theta((i+n_sources));
    
  }
  
  double totx = x + gammaprior + thetanorm;
  
  return totx;
  
}









