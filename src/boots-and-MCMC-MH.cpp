#include <Rcpp.h>
using namespace Rcpp;

//' @title Bootstrap in Rcpp
//' @name BootstrapC
//' @description The bootstrap in Rcpp which can reduce the computation time largely when comparing with R.
//' @param x Data
//' @param B the sampling times
//' @return a list of Bootstrap estimator
//' @examples
//' \dontrun{
//' data <- rnorm(100)
//' bootstrap_results <- bootsC(data, 1000)
//' }
//' @export
// [[Rcpp::export]]
NumericVector bootsC(NumericVector x,int B) {
  NumericVector thetastar(B);
  double theta = mean(x);
  int n = x.size();

  for(int b = 0; b < B; b++) {
    NumericVector xstar = Rcpp::sample(x, n, true);
    thetastar[b] = mean(xstar);
  }

  double bias = mean(thetastar) - theta;
  double se_boot = sd(thetastar);
  double se_samp = sd(x) / sqrt(n);

  return NumericVector::create(Named("bias") = bias,
                               Named("se.boot") = se_boot,
                               Named("se.samp") = se_samp);
}


double f(double x, double sigma) {
  if (x < 0) return 0;
  return (x / pow(sigma, 2)) * exp(-pow(x, 2) / (2 * pow(sigma, 2)));
}

//' @title MCMC Metropolis-Hastings sampler generating random number of Rayleigh distribution in Rcpp
//' @name MCMCMHC
//' @description Metropolis-Hastings sampler generates random numbers from the Rayleigh distribution in Rcpp which can reduce the computation time largely when comparing with R.
//' @param m random numbers size
//' @param sigma Rayleigh distribution's parameter
//' @param b burn-in size
//' @return a random sample from Rayleigh distribution
//' @examples
//' \dontrun{
//' RayRandom=RayleighC(m=2000,sigma=4,b=100)
//' }
//' @export
// [[Rcpp::export]]
NumericVector RayleighC(int m, double sigma, int b) {
  NumericVector x(m);
  x[0] = R::rchisq(1); 
  int k = 0;
  NumericVector u = runif(m);
   
  Environment stats("package:stats");
  Function dchisq = stats["dchisq"];
  Function rchisq = stats["rchisq"];
   
  for (int i = 1; i < m; ++i) {
    double xt = x[i - 1];
    double y = R::rchisq(1); 
    double num = f(y, sigma) * as<double>(dchisq(xt, y));
    double den = f(xt, sigma) * as<double>(dchisq(y, xt));
     
  if (u[i] <= num / den) {
    x[i] = y;
    } else {
    x[i] = xt;
    k++;
    }
  }
  return x[Range(b, m - 1)]; 
 }