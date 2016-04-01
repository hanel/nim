#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List cppREG(NumericVector p, NumericVector xpar, NumericMatrix xi, float g, float k){

  int n = xpar.size(); int s = xi.ncol();
  float pom1, pom2, pom3;
  NumericMatrix loc(n, s); NumericMatrix scale(n, s); NumericMatrix shape(n, s);

  for(int i = 0; i < n; ++i) {

    pom1 = p[0] + p[1] * xpar[i];
    pom2 = exp(p[2] + p[3] * xpar[i]);
    pom3 = p[4] ;

    for(int j = 0; j < s; ++j) {

      loc(i, j) = xi(i, j) * pom1;
      scale(i, j) = loc(i, j) * pom2;
      shape(i, j) = pom3;

    }}

  return List::create(Rcpp::Named("loc") = loc, Rcpp::Named("scale") = scale, Rcpp::Named("shape") = shape);
}

// [[Rcpp::export]]
List cppREG_lin(NumericVector p, NumericVector xpar, NumericMatrix xi, float g, float k){

  int n = xpar.size(); int s = xi.ncol();
  float pom1, pom2, pom3;
  NumericMatrix loc(n, s); NumericMatrix scale(n, s); NumericMatrix shape(n, s);

  for(int i = 0; i < n; ++i) {

    pom1 = p[0] + p[1] * xpar[i];
    pom2 = exp(p[2] + p[3] * xpar[i]);
    pom3 = p[4] + p[5] * xpar[i] ;

    for(int j = 0; j < s; ++j) {

      loc(i, j) = xi(i, j) * pom1;
      scale(i, j) = loc(i, j) * pom2;
      shape(i, j) = pom3;

    }}

  return List::create(Rcpp::Named("loc") = loc, Rcpp::Named("scale") = scale, Rcpp::Named("shape") = shape);
}


// [[Rcpp::export]]
List cppREG_constK(NumericVector p, NumericVector xpar, NumericMatrix xi, float g, float k){

  int n = xpar.size(); int s = xi.ncol();
  float pom1, pom2;
  NumericMatrix loc(n, s); NumericMatrix scale(n, s); NumericMatrix shape(n, s);

  for(int i = 0; i < n; ++i) {

    pom1 = p[0] + p[1] * xpar[i];
    pom2 = exp(p[2] + p[3] * xpar[i]);

    for(int j = 0; j < s; ++j) {

      loc(i, j) = xi(i, j) * pom1;
      scale(i, j) = loc(i, j) * pom2;
      shape(i, j) = k;

    }}

  return List::create(Rcpp::Named("loc") = loc, Rcpp::Named("scale") = scale, Rcpp::Named("shape") = shape);

}

// [[Rcpp::export]]
List cppREG_constGK(NumericVector p, NumericVector xpar, NumericMatrix xi, float g, float k){

  int n = xpar.size(); int s = xi.ncol();
  float pom1;
  NumericMatrix loc(n, s); NumericMatrix scale(n, s); NumericMatrix shape(n, s);

  for(int i = 0; i < n; ++i) {

    pom1 = p[0] + p[1] * xpar[i];

    for(int j = 0; j < s; ++j) {

      loc(i, j) = xi(i, j) * pom1;
      scale(i, j) = loc(i, j) * exp(g);
      shape(i, j) = k;

    }}

  return List::create(Rcpp::Named("loc") = loc, Rcpp::Named("scale") = scale, Rcpp::Named("shape") = shape);

}

// [[Rcpp::export]]
List cppREG_const(NumericVector p, NumericVector xpar, NumericMatrix xi, float g, float k){

  int n = xpar.size(); int s = xi.ncol();
  NumericMatrix loc(n, s); NumericMatrix scale(n, s); NumericMatrix shape(n, s);

  for(int i = 0; i < n; ++i) {

    for(int j = 0; j < s; ++j) {

      scale(i, j) = xi(i, j) * exp(g);
      shape(i, j) = k;

    }}

  return List::create(Rcpp::Named("loc") = loc, Rcpp::Named("scale") = scale, Rcpp::Named("shape") = shape);

}


// [[Rcpp::export]]
List cppXI0(float p, NumericVector xi, NumericVector g, NumericVector k){
  int n = xi.size();
  NumericMatrix loc(n, 1); NumericMatrix scale(n, 1);

  for (int i = 0; i < n; i++){
    loc(i, 0) = p * xi(i);
    scale(i, 0) = loc(i, 0) * exp(g(i));
  }
  k.attr("dim") = Dimension(n, 1);
  return List::create(Rcpp::Named("loc") = loc, Rcpp::Named("scale") = scale, Rcpp::Named("shape") = k);
}


// [[Rcpp::export]]
double nll_gev(NumericVector par, Function fpar, NumericMatrix data, NumericVector w){
  List pmat = fpar(par);
  double nll = 0, y, z, n = 0;

  NumericMatrix  loc = pmat["loc"];
  NumericMatrix scale = pmat["scale"];
  NumericMatrix shape = pmat["shape"];

  for (int i = 0; i < data.nrow(); i++){
    for (int j = 0; j < data.ncol(); j++){
      
      if (NumericMatrix::is_na(data(i, j))) {n = 0;} else {
        if (scale(i, j) <= 0) {return 1e20;}
        y = (data(i, j) - loc(i, j)) / scale(i, j);
        z = 1 + shape(i, j) * y;
        if (z <= 0) {return 1e20;}
        if ((fabs(shape(i, j)) < 1e-6)) {n = y + exp(-y);}
        else {n = (1 + 1 / shape(i, j)) * log(z) + pow(z, (-1 / shape(i, j)));}
      }
      nll = nll + w[i] * (n + log(scale(i, j)));
    }
  }

  return nll;
}
