// includes from the plugin
//#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
#include <Rcpp.h>


#ifndef BEGIN_RCPP
#define BEGIN_RCPP
#endif

#ifndef END_RCPP
#define END_RCPP
#endif

using namespace Rcpp;
using namespace std;

// declarations
extern "C" {
SEXP coneACpp( SEXP y, SEXP amat, SEXP face) ;
}

// definition
SEXP coneACpp( SEXP y, SEXP amat, SEXP face){
BEGIN_RCPP

    Rcpp::NumericVector y_(y);
    Rcpp::NumericMatrix amat_(amat);
    int n = y_.size(), m = amat_.nrow();
    arma::mat namat(amat_.begin(), m, n, false);
    arma::colvec ny(y_.begin(), n, false);
    float sm = 1e-8;
    arma::colvec h(m); h.fill(0);
    arma::colvec obs = arma::linspace(0, m-1, m);
    int check = 0;
    arma::mat amat_in = namat;
//new:
    arma::colvec face_ = Rcpp::as<arma::colvec>(face);
    int nf = face_.n_rows;

    for(int im = 0; im < m; im ++){
        arma::mat nnamat = namat.row(im) * namat.row(im).t();
        amat_in.row(im) = namat.row(im) / sqrt(nnamat(0,0));
    }

    arma::mat delta = -amat_in;
    arma::colvec b2 = delta * ny;
    arma::colvec theta(n); theta.fill(0);

    if(max(b2) > 2 * sm){
        int i = min(obs.elem(find(b2 == max(b2))));
        h(i) = 1;
//new:
        if(!face_.is_empty()) {
            for (int i = 0; i < nf; i ++) {
                int posi = face_(i) - 1;
                h(posi) = 1;
            }
        }
    }
    else{check = 1;}

    int nrep = 0;
//new:
    arma::mat xmat_use;
    while(check == 0 & nrep < (n * n)){
        nrep ++ ;
        arma::colvec indice = arma::linspace(0, delta.n_rows - 1, delta.n_rows);
        indice = indice.elem(find(h == 1));
        arma::mat xmat(indice.n_elem, delta.n_cols); xmat.fill(0);

        for(int k = 0; k < indice.n_elem; k ++){
            xmat.row(k) = delta.row(indice(k));
        }

        arma::colvec a = solve(xmat * xmat.t(), xmat * ny);
        arma::colvec avec(m); avec.fill(0);

        if(min(a) < (-sm)){
            avec.elem(find(h == 1)) = a;
            int i = min(obs.elem(find(avec == min(avec))));
            h(i) = 0;
            check = 0;
        }

        else{
            check = 1;
            theta = xmat.t() * a;
            b2 = delta * (ny - theta)/n;

            if(max(b2) > 2 * sm){
                int  i = min(obs.elem(find(b2 == max(b2))));
                h(i) = 1;
                check = 0;
            }
        }
//new:
       xmat_use = xmat;
    }

    //if(nrep > (n * n - 1)){Rcpp::Rcout << "Fail to converge in coneproj!Too many steps! Number of steps:" << nrep << std::endl;}

    return wrap(Rcpp::List::create(Rcpp::Named("thetahat") = ny - theta, Named("xmat") = xmat_use, Named("dim") = n - sum(h), Named("nrep") = nrep, Named("h") = h));

END_RCPP
}

// declarations
extern "C" {
SEXP coneBCpp( SEXP y, SEXP delta, SEXP vmat, SEXP face) ;
}

// definition
SEXP coneBCpp( SEXP y, SEXP delta, SEXP vmat, SEXP face){
BEGIN_RCPP

    Rcpp::NumericVector y_(y);
    Rcpp::NumericMatrix delta_(delta);
    arma::mat nvmat = Rcpp::as<arma::mat>(vmat);
    int n = y_.size(), m = delta_.nrow(), p = nvmat.n_cols;
    arma::colvec ny(y_.begin(), n, false);
    arma::mat ndelta(delta_.begin(), m, n, false);
    arma::mat a;
    arma::mat sigma;
    arma::colvec h;
//new:
    arma::colvec face_ = Rcpp::as<arma::colvec>(face);
    int nf = face_.n_rows;
    arma::colvec obs;
    arma::mat theta(n, 1);

    float sm = 1e-8;
//new: test!
    //float sm = 1e-5;
    int check = 0;

    arma::colvec scalar(m);
    arma::mat delta_in = ndelta;

    for(int im = 0; im < m; im ++){
        arma::mat nndelta = ndelta.row(im) * ndelta.row(im).t();
        scalar(im) = sqrt(nndelta(0,0));
        delta_in.row(im) = ndelta.row(im) / scalar(im);
    }

    if(nvmat.is_empty()){
       p = p - 1;
       sigma.set_size(m, n);
       sigma = delta_in;
       h.set_size(m);
       h.fill(0);
       obs.set_size(m);
       obs = arma::linspace(0, m - 1, m);
       theta.fill(0);
    }

    if(!nvmat.is_empty()){
        sigma.set_size(m + p, n);
        sigma.rows(0, p - 1) = nvmat.t(); sigma.rows(p, m + p - 1) = delta_in;
        h.set_size(m + p);
        h.fill(0);
        for(int i = 0; i < p; i ++){
          h(i) = 1;
        }
        obs.set_size(m + p);
        obs = arma::linspace(0, m + p - 1, m + p);
        theta = nvmat * solve(nvmat.t() * nvmat, nvmat.t() * ny);
    }

    arma::colvec b2 = sigma * (ny - theta) / n;

    if(max(b2) > 2 * sm){
        int i = min(obs.elem(find(b2 == max(b2))));
        h(i) = 1;
//new:
        if(!face_.is_empty()) {
            for (int i = 0; i < nf; i ++) {
                int posi = face_(i) - 1;
                h(posi) = 1;
            }
        }
    }

    int nrep = 0;

    if(max(b2) <= 2 * sm){
        check = 1;
        theta.fill(0);

        if(nvmat.is_empty()){
           a.set_size(m, 1); a.fill(0);
        }

        if(!nvmat.is_empty()){
           a.set_size(p, 1);
           a = solve(nvmat.t() * nvmat, nvmat.t() * ny);
           theta = nvmat * solve(nvmat.t() * nvmat, nvmat.t() * ny);
        }
        arma::colvec avec(m + p); avec.fill(0);
        if(!nvmat.is_empty()){
          avec.elem(find(h == 1)) = a;
        }
        return wrap(Rcpp::List::create(Named("yhat") = theta, Named("coefs") = avec, Named("nrep") = nrep, Named("dim") = sum(h)));
    }
//new:
    //int upper = n*n - 1;
    //int upper = 1000;
   // while(check == 0 & nrep < (n * n)){
//double sc = 0;
   while(check == 0 & nrep < 1e+6){
        nrep ++;
        //if(nrep > (n * n)){
          // throw (Rcpp::exception("Fail to converge in coneproj! nrep > n^2 !"));
        //}
        arma::colvec indice = arma::linspace(0, sigma.n_rows-1, sigma.n_rows);
        indice = indice.elem(find(h == 1));
        arma::mat xmat(indice.n_elem, sigma.n_cols); xmat.fill(0);

        for(int k = 0; k < indice.n_elem; k ++){
            xmat.row(k) = sigma.row(indice(k));
        }
 	//double sc = arma::norm(xmat * xmat.t(), 2);
        a = solve(xmat * xmat.t(), xmat * ny);
//new:
       if (a.n_elem > p) {
            arma::colvec a_sub(a.n_elem - p);

            for(int i = p; i <= a.n_elem - 1; i ++){
                a_sub(i-p) = a(i);
            }

            if(min(a_sub) < (- sm)){
                arma::colvec avec(m + p); avec.fill(0);
                avec.elem(find(h == 1)) = a;
                arma::colvec avec_sub(m);

                for(int i = p; i <= p + m - 1; i ++){
                    avec_sub(i-p) = avec(i);
                }

                int i = max(obs.elem(find(avec == min(avec_sub))));
                h(i) = 0;
                check = 0;
            }

            if(min(a_sub) > (-sm)){
                check = 1;
                theta = xmat.t() * a;
                b2 = sigma * (ny - theta) / n;
//arma::mat sc0 = sqrt(b2.t() * b2);
//sc = sqrt(sc0(0,0));
//sc = arma::as_scalar(b2.t() * b2);
                //if(max(b2) > 2 * sc * sm){
		//if((max(b2) * sc) > 2 * sm){
                if(max(b2) > 2 * sm){
                    int i = min(obs.elem(find(b2 == max(b2))));
                    check = 0;
                    h(i) = 1;
                }
            }
        } else {
            check = 1;
       }
//new: avoid the mismatch problem
       if (nrep == 1e+6) {
            arma::colvec indiceEnd = arma::linspace(0, sigma.n_rows-1, sigma.n_rows);
       	    indiceEnd = indiceEnd.elem(find(h == 1));
            arma::mat xmat(indiceEnd.n_elem, sigma.n_cols); xmat.fill(0);
            for(int k = 0; k < indiceEnd.n_elem; k ++){
                xmat.row(k) = sigma.row(indiceEnd(k));
            }
            //sc = norm(xmat * xmat.t(), 2);
            a = solve(xmat * xmat.t(), xmat * ny);
            theta = xmat.t() * a;
       }
   }

   arma::colvec avec(m + p); avec.fill(0);
   avec.elem(find(h == 1)) = a;
   arma::colvec avec_orig(m + p); avec_orig.fill(0);

   for(int i = 0; i < p; i ++){
      avec_orig(i) = avec(i);
   }

   for(int i = p; i < (m + p); i ++){
      avec_orig(i) = avec(i) / scalar(i - p);
   }
   // if(nrep > (n * n - 1)){Rcpp::Rcout << "Fail to converge in coneproj!Too many steps! Number of steps:" << nrep << std::endl;}
   return wrap(Rcpp::List::create(Named("yhat") = theta, Named("coefs") = avec_orig, Named("nrep") = nrep, Named("dim") = sum(h)));

END_RCPP
}

// declarations
extern "C" {
SEXP qprogCpp( SEXP q, SEXP c, SEXP amat, SEXP b, SEXP face) ;
}

// [[Rcpp::export]]
SEXP qprogCpp( SEXP q, SEXP c, SEXP amat, SEXP b){
BEGIN_RCPP

    Rcpp::NumericVector c_(c);
    Rcpp::NumericMatrix q_(q);
    Rcpp::NumericMatrix amat_(amat);
    Rcpp::NumericVector nb(b);
    int n = c_.size(), m = amat_.nrow();
    arma::colvec nc(c_.begin(), n, false);
    arma::mat namat(amat_.begin(), m, n, false);
    arma::mat nq(q_.begin(), n, n, false);
    bool constr = is_true(any( nb != 0 ));
    arma::colvec theta0(n);
    arma::colvec nnc(n);

    //arma::colvec face_ = Rcpp::as<arma::colvec>(face);
    //int nf = face_.n_rows;

    if(constr){
        arma::colvec b_(nb.begin(), m, false);
        theta0 = solve(namat, b_);
        nnc = nc - nq * theta0;
    }

    else{nnc = nc;}

    arma::mat preu = chol(nq);
    arma::mat u = trimatu(preu);
    arma::colvec z = inv(u).t() * nnc;
    arma::mat atil = namat * inv(u);

    float sm = 1e-8;
    arma::colvec h(m); h.fill(0);
    arma::colvec obs = arma::linspace(0, m-1, m);
    int check = 0;

    for(int im = 0; im < m; im ++){
        arma::mat atilnorm = atil.row(im) * atil.row(im).t();
        atil.row(im) = atil.row(im) / sqrt(atilnorm(0,0));
    }

    arma::mat delta = -atil;
    arma::colvec b2 = delta * z;
    arma::colvec phi(n); phi.fill(0);

    if(max(b2) > 2 * sm){
        int i = min(obs.elem(find(b2 == max(b2))));
        h(i) = 1;
//new:
        //if(!face_.is_empty()) {
        //    for (int i = 0; i < nf; i ++) {
        //        int posi = face_(i) - 1;
        //        h(posi) = 1;
        //    }
        //}
    }

    else{check = 1;}

    int nrep = 0;
//new:
    arma::mat xmat_use;
    while(check == 0 & nrep < (n * n)){
        nrep ++ ;
       // if(nrep > (n * n)){
         //  throw (Rcpp::exception("Fail to converge in coneproj! nrep > n^2 !"));}
        arma::colvec indice = arma::linspace(0, delta.n_rows - 1, delta.n_rows);
        indice = indice.elem(find(h == 1));
        arma::mat xmat(indice.n_elem, delta.n_cols); xmat.fill(0);

        for(int k = 0; k < indice.n_elem; k ++){
        xmat.row(k) = delta.row(indice(k));
        }

        arma:: colvec a = solve(xmat * xmat.t(), xmat * z);
        arma:: colvec avec(m); avec.fill(0);

        if(min(a) < (-sm)){
            avec.elem(find(h == 1)) = a;
            int i = min(obs.elem(find(avec == min(avec))));
            h(i) = 0;
            check = 0;
        }

        else{
            check = 1;
            phi = xmat.t() * a;
            b2 = delta * (z - phi)/n;

            if(max(b2) > 2 * sm){
                int  i = min(obs.elem(find(b2 == max(b2))));
                h(i) = 1;
                check = 0;
            }
        }
	//new:
    	xmat_use = xmat;
    }

    arma::colvec thetahat = solve(u, z - phi);

    if(constr){
        thetahat = thetahat + theta0;
    }

    // if(nrep > (n * n - 1)){Rcpp::Rcout << "Fail to converge in qprog!Too many steps! Number of steps:" << nrep << std::endl;}

    return wrap(Rcpp::List::create(Rcpp::Named("thetahat") = thetahat, Named("xmat") = xmat_use, Named("dim") = n - sum(h), Named("nrep") = nrep, Named("h") = h));

END_RCPP
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat multiplyXtXBySparse_c(int n, arma::uvec ind, arma::mat val, arma::vec phi, arma::vec w){
  int p = val.n_cols;
  int as = ind.size();
  arma::uvec indrev = arma::zeros<arma::uvec>(n-1);
  for (int i=0; i<(n-1); ++i){
    indrev(i) = n-2-i;
  }
  arma::uvec o = sort_index(ind);
  ind = sort(ind);
  val = val.rows(o);
  arma::vec phicum = cumsum(phi%phi);
  arma::mat r = val.each_col() % w(ind);
  arma::mat s1 = r.each_col() % phicum(ind);
  arma::mat cumsums1 = cumsum(s1);
  arma::rowvec s = cumsums1.row(as-1)/n;
  arma::mat matrixT = arma::zeros<arma::mat>(n-1, p);
  matrixT.rows(indrev(ind)) = r;
  matrixT = cumsum(matrixT);
  matrixT = matrixT.rows(indrev);
  arma::mat u = matrixT.each_row()-s;
  arma::vec phisq = phi%phi;
  u.each_col() %= phisq.subvec(0,n-2);
  arma::mat U = cumsum(u);
  arma::mat C = U.each_col() % w;
  return C;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat leftMultiplyByXt_c(arma::mat Y, arma::vec phi, arma::vec w){
  int n = Y.n_rows;
  int p = Y.n_cols;
  Y.each_col() %= phi;
  arma::mat u = cumsum(Y);
  arma::mat c = arma::zeros<arma::mat>(n-1,p);
  arma::vec phicum = cumsum(phi%phi);
  for (int i=0; i<p; ++i){
    arma::vec x = u.col(i);
    c.col(i) = w % (phicum.subvec(0,n-2) * (x(n-1)/n) - x.subvec(0,n-2));
  }
  return c;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat leftMultiplyByInvXAtXA_c(double n, arma::uvec ind, arma::mat val, arma::vec phi, arma::vec w){
  int a = val.n_rows;
  int p = val.n_cols;
  int as = ind.size();
  arma::mat r = arma::zeros<arma::mat>(a, p);
  if (a!=0){
    arma::vec v1 = arma::zeros<arma::vec>(2+as);
    v1(0) = 0;
    arma::vec phicum = cumsum(phi%phi);
    v1.subvec(1,as)=phicum(ind);
    v1(as+1) = n;
    arma::vec v = diff(v1);
    arma::vec d = w(ind);
    arma::mat R = arma::zeros<arma::mat>(a+2, p);
    val.each_col() /= d;
    R.row(0).fill(0);
    R.rows(1,a) = val;
    R.row(a+1).fill(0);
    arma::mat gamma = diff(R);
    arma::mat delta = gamma;
    delta.each_col() /= v;
    r = -diff(delta);
    r.each_col() /= d;
  }
  return r;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec ic_c(int k, arma::mat Y, arma::vec phi, arma::rowvec xi, arma::mat theta, int p, int n){
  arma::mat dis = arma::zeros<arma::mat>(p, n);
  int np = n*p;
  double rss = 0;
  for(int i = 0; i < p; ++i){
    for(int j = 0; j < n; ++j){
      dis(i, j) = Y(i, j) - (theta(i, j) + xi(j))*phi(i);
      rss += dis(i, j)*dis(i, j);
    }
  }
  rss = rss/(2*np);
  const double bic_p = log(p) * ((k+1)*n + (p-k-1) + 2*k)/np;
  const double temp1 = 2*((k+1)*n+(p-k-1)+3*k);
  arma::vec out = arma::zeros<arma::vec>(3);
  out[0] = rss;
  out[1] = log(rss) + temp1/np;
  out[2] = log(rss) + bic_p;
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec defaultWeights_c(int p){
  arma::vec out = arma::zeros<arma::vec>(p-1);
  for(int i = 0; i < (p-1); ++i){
    double temp = (i+1) * (p-i-1);
    out[i] = (sqrt(p/temp));
  }
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List doLars_c(arma::mat Y, int K, arma::vec phi, arma::vec wts, int p, int n, double epsilon = 1e-09){
  arma::vec res_lambda = arma::zeros<arma::vec>(K);
  arma::uvec res_bkp = arma::zeros<arma::uvec>(K);
  arma::mat res_value = arma::zeros<arma::mat>(K,n);
  arma::mat res_value_old = arma::zeros<arma::mat>(K,n);
  arma::mat result = arma::zeros<arma::mat>(K, n+2);
  arma::uvec AS(1);
  arma::uvec AS0(1);
  arma::mat c = leftMultiplyByXt_c(Y, phi, wts);
  for (int ii=0; ii<K; ++ii){
    arma::mat cc = c%c; //381 x 32
    arma::mat ccrowsum = cumsum(cc, 1);
    arma::vec cNorm = ccrowsum.col(n-1); //381x1
    double bigcHat = cNorm.max();
    if (ii==0){
      AS(0) = cNorm.index_max();
      res_bkp[ii] = cNorm.index_max();
    }
    arma::uvec I = sort_index(AS);
    AS0 = AS;
    AS = AS(I);
    arma::mat w = leftMultiplyByInvXAtXA_c(p, AS, c.rows(AS), phi, wts);
    arma::mat a = multiplyXtXBySparse_c(p, AS, w, phi, wts);
    arma::mat aa = a%a;
    arma::mat cumsumaa = cumsum(aa, 1);
    arma::vec rowsumsa2 = cumsumaa.col(aa.n_cols-1);
    arma::vec a1 = bigcHat - rowsumsa2;
    arma::mat u = a%c;
    arma::mat cumsumu = cumsum(u, 1);
    arma::vec rowsumsu = cumsumu.col(u.n_cols-1);
    arma::vec a2 = bigcHat-rowsumsu;
    arma::vec a3 = bigcHat-cNorm;
    arma::mat gammaTemp = arma::zeros<arma::mat>(p-1, 2);
    gammaTemp.fill(NA_REAL);
    arma::uvec subset = find(a1 > epsilon);
    arma::vec delta = arma::zeros<arma::vec>(subset.size());
    arma::uvec onetemp = arma::zeros<arma::uvec>(1);
    arma::uvec twotemp = arma::ones<arma::uvec>(1);
    delta = a2(subset)%a2(subset)-a1(subset)%a3(subset);
    arma::uvec deltaneg = find(delta<0);
    arma::uvec deltapos = find(delta>=0);
    arma::uvec delta_neg = subset(deltaneg);
    gammaTemp.submat(delta_neg, onetemp).fill(NA_REAL);
    gammaTemp.submat(delta_neg, twotemp).fill(NA_REAL);
    arma::uvec delta_pos = subset(deltapos);
    gammaTemp.submat(delta_pos, onetemp) = (a2(delta_pos) +
      sqrt(delta(deltapos)))/a1(delta_pos);
    gammaTemp.submat(delta_pos, twotemp) = (a2(delta_pos) -
      sqrt(delta(deltapos)))/a1(delta_pos);
    subset = find((a1 <= epsilon) && (a2 > epsilon));
    arma::vec whattofill = a3(subset)/(2*a2(subset));
    for (int i=0; i<subset.size(); ++i){
      gammaTemp.row(subset[i]).fill(whattofill(i));
    }
    double maxg = gammaTemp.max() + 1;
    subset = find((a1 <= epsilon) && (a2 <= epsilon));
    gammaTemp.rows(subset).fill(maxg);
    gammaTemp.rows(AS).fill(maxg);
    gammaTemp(find(gammaTemp<=0)).fill(maxg);
    double gamma = gammaTemp.min();
    int idx = gammaTemp.index_min();
    arma::uvec nexttoadd(1);
    nexttoadd(0) = (idx)%(p-1);
    double nexttoadd_double = nexttoadd(0);
    res_lambda[ii] = sqrt(bigcHat);
    res_value.fill(0);
    res_value.rows(I) = gamma*w;
    if(ii > 0){
      res_value += res_value_old;
    }
    if (ii < (K-1)) {
      AS = join_cols(AS0, nexttoadd);
      res_bkp(ii+1) = nexttoadd_double;
      c = c - gamma * a;
    }
    res_value_old = res_value;
  }
  arma::uvec one = arma::ones<arma::uvec>(res_bkp.size());
  res_bkp += one;
  return Rcpp::List::create(
    Rcpp::Named("bkp") = res_bkp,
    Rcpp::Named("lambda") = res_lambda,
    Rcpp::Named("value") = res_value
  );
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List cnvJoint(arma::mat Y, arma::vec wts, int steps, int maxloop=10, bool verbose=false){
  /* Initialize phi*/
  const int p = Y.n_rows;/*number of probes*/
  const int n = Y.n_cols; /*number of samples*/
  arma::vec phi = arma::ones<arma::vec>(p);
  arma::rowvec xi = arma::ones<arma::rowvec>(n);
  int loop = 0;
  arma::vec phi_old = phi + 1e+10;
  double error = 1e+10;
  //double error2 = 1e+10;
  while (error > 1e-5 & loop < 10){
    arma::mat theta1 = (phi.t() * Y / p).t();
    phi = Y * theta1 / accu(theta1 % theta1);
    phi = phi / sqrt(accu(phi%phi)/p);
    arma::vec diff = phi_old - phi;
    error = sqrt(accu(diff%diff)/p);
    phi_old = phi;
    loop += 1;
  }
  /*Initialize variables to save*/
  arma::vec aic_list = arma::zeros<arma::vec>(steps-1);
  arma::vec bic_list = arma::zeros<arma::vec>(steps-1);
  arma::vec rss_list = arma::zeros<arma::vec>(steps-1);
  arma::vec loop_list = arma::zeros<arma::vec>(steps-1);
  arma::mat bkp_list = arma::zeros<arma::mat>(steps-1, steps-1);
  arma::mat phi_list = arma::zeros<arma::mat>(p, steps-1);
  arma::mat xi_list = arma::zeros<arma::mat>(steps-1, n);
  arma::cube theta_list = arma::zeros<arma::cube>(p, n, steps-1);
  /* start computing the result for each k*/
  for (int k = 1; k < (steps); ++k){
    if(verbose==TRUE) {cout <<"looking for "<< "k = "<< k << " changepoints.. \n";}
    /* For each given k, we compute bkp, aic, etc*/
    loop=0; error=1e5;
    /* Initialize variable for next loop */
    arma::mat theta = arma::zeros<arma::mat>(p, n);
    arma::vec bkp = arma::zeros<arma::vec>(k);
    arma::vec lambda = arma::zeros<arma::vec>(k);
    arma::mat value = arma::zeros<arma::mat>(k, n);
    arma::vec sorted_bkp = arma::zeros<arma::vec>(k);
    arma::vec phi_new = arma::zeros<arma::vec>(p);
    arma::rowvec xi_new = arma::zeros<arma::rowvec>(n);
    arma::vec q = arma::zeros<arma::vec>(maxloop);
    bool alternating = TRUE;
    while (error > 1e-3 & loop<maxloop & alternating){
      arma::mat Ynew = Y - phi*xi;
      Rcpp::List lars = doLars_c(Ynew, k, phi, wts, p, n);
      bkp = as<arma::vec>(lars["bkp"]);
      lambda = as<arma::vec>(lars["lambda"]);
      value = as<arma::mat>(lars["value"]);
      arma::uvec ord = sort_index(bkp);
      sorted_bkp = sort(bkp);
      arma::uvec sorted(k);
      for (int i=0; i<k; ++i){
        sorted(i) = sorted_bkp(i)-1;
      }
      arma::mat delta2 = arma::zeros<arma::mat>(p-1, n);
      delta2.rows(sorted) = value.rows(ord);
      arma::mat delta3 = arma::zeros<arma::mat>(p-1, n);
      delta3 = delta2.each_col() % wts;
      arma::mat delta1 = phi.t() * Ynew;
      arma::vec phisq = phi%phi;
      arma::vec cumsumphi = flipud(cumsum(flipud(phisq)));
      cumsumphi.shed_row(0);
      delta1 = (delta1-(cumsumphi.t() * delta3))/accu(phisq);
      arma::mat delta = join_cols(delta1, delta3);
      theta = delta;
      for (int i=1; i<p; ++i){
        theta.row(i) = theta.row(i-1) + delta.row(i);
      }
      //partition
      arma::mat partition(k+1, p);
      partition.fill(NA_REAL);
      if (k==1){
        for (int i=0; i < bkp(0); ++i){
          partition(0, i) = i+1;
        }
        for (int i=bkp(0); i < p; ++i){
          partition(1, i-bkp(0)) = i+1;
        }
      }
      else{
        for (int i=0; i < sorted(0)+1; ++i){
          partition(0, i) = i+1;
        }
        for (int i=sorted(k-1)+1; i<p; ++i){
          partition(k, i-sorted(k-1)-1) = i+1;
        }
        for (int i=0; i < (k-1); ++i){
          for (int j=sorted(i)+1; j<sorted(i+1)+1; ++j){
            partition(i+1, j-sorted(i)-1) = j+1;
          }
        }
      }

      /*//Update phi : this was before the phi>0 constraint
      arma::mat thetaY = theta%Y + Y.each_row()%xi;
      arma::mat thetaYcumsum = cumsum(thetaY, 1);
      arma::vec thetaYrowsum = thetaYcumsum.col(n-1);
      arma::mat newtheta = theta.each_row() + xi;
      arma::mat thetathetacumsum = cumsum(newtheta%newtheta, 1);
      arma::vec thetasqrowsum = thetathetacumsum.col(n-1);
      phi_new = thetaYrowsum / thetasqrowsum;
      for (int i=0; i<(k+1); ++i){
        arma::rowvec ith_partition = partition.row(i);
        arma::vec where = ith_partition(find_finite(ith_partition));
        arma::uvec uwhere = arma::zeros<arma::uvec>(where.size());
        for (int j = 0; j < where.size(); ++j){
          uwhere(j) = where(j)-1;
        }
        arma::vec phi_temp = phi_new(uwhere);
        double phitempsqmean = mean(phi_temp%phi_temp);
        phi_new(uwhere) = phi_temp/sqrt(phitempsqmean);
      }
      arma::vec diff = phi-phi_new;
      error = accu(diff%diff) / p;
      error = sqrt(error);
      phi = phi_new;*/

      arma::mat newtheta = theta.each_row() + xi;
      arma::mat H = arma::zeros<arma::mat>(p,p);
      arma::mat Thetasqcumsum = cumsum(newtheta%newtheta, 1);
      arma::mat Ytheta = Y%newtheta;
      arma::mat Ythetacumsum = cumsum(Ytheta, 1);
      arma::vec C = arma::zeros<arma::vec>(p);
      arma::mat identity = arma::eye<arma::mat>(p,p);
      SEXP id = wrap(identity);
      arma::mat zeros = arma::zeros<arma::vec>(p);
      SEXP zer = wrap(zeros);
      for (int i=0; i<p; ++i){
        H(i,i) = Thetasqcumsum(i, n-1);
        C(i) = Ythetacumsum(i, n-1);
      }
      SEXP H2 = wrap(H);
      SEXP C2 = wrap(C);
      Rcpp::List sphi = qprogCpp(H2, C2, id, zer);
      phi_new = as<arma::vec>(sphi["thetahat"]);
      for (int i=0; i<(k+1); ++i){
        arma::rowvec ith_partition = partition.row(i);
        arma::vec where = ith_partition(find_finite(ith_partition));
        arma::uvec uwhere = arma::zeros<arma::uvec>(where.size());
        for (int j = 0; j < where.size(); ++j){
          uwhere(j) = where(j)-1;
        }
        arma::vec phi_temp = phi_new(uwhere);
        if(mean(phi_temp%phi_temp)!=0){
          double phitempsqmean = mean(phi_temp%phi_temp);
          phi_new(uwhere) = phi_temp/sqrt(phitempsqmean);
        }
        if(mean(phi_temp%phi_temp)==0){
          phi_new(uwhere) = phi_temp;
        }
      }
      arma::vec diff = phi-phi_new;
      error = accu(diff%diff) / p;
      error = sqrt(error);
      phi = phi_new;

      //Update xi
      arma::mat xitemp1 = (Y - theta.each_col()%phi);
      arma::mat xitemp2 = xitemp1.each_col() % phi;
      arma::mat xitempcumsum = cumsum(xitemp2, 0);
      arma::rowvec colsumxi = xitempcumsum.row(p-1);
      xi = colsumxi / accu(phi%phi);
      arma::rowvec meanxi = arma::zeros<arma::rowvec>(n);
      meanxi.fill(accu(xi)/n);
      xi = xi - meanxi;

      loop += 1;

      //Compute the optimization objective Q1, check convergence
      arma::vec res = ic_c(k,Y,phi,xi,theta,p,n);
      double q1 = res(0) * 2* n * p;
      for (int i = 0; i < (p-1); ++i){
        arma::rowvec temp = theta.row(i+1)-theta.row(i);
        q1 += lambda(k-1) * sqrt(accu(temp%temp)) / wts(i);
      }
      q(loop-1) = q1;
      if(loop>5){
        alternating = (abs(q(loop-1)-q(loop-3)) + abs(q(loop-2) - q(loop-4)) > 1e-3);
      }
    }
    arma::vec res = ic_c(k, Y, phi, xi, theta, p, n);
    rss_list(k-1) = res(0);
    aic_list(k-1) = res(1);
    bic_list(k-1) = res(2);
    loop_list(k-1) = loop;
    for (int i=0; i<k; ++i){
      bkp_list(k-1, i) = sorted_bkp(i);
    }
    theta_list.slice(k-1) = theta;
    phi_list.col(k-1) = phi;
    xi_list.row(k-1) = xi;
  }
  int aic_k = aic_list.index_min();
  int bic_k = bic_list.index_min();

  /*Organize result*/
  arma::rowvec tempaicbkp = bkp_list.row(aic_k);
  arma::rowvec tempbicbkp = bkp_list.row(bic_k);
  arma::vec aicbkp = arma::zeros<arma::vec>(aic_k+1);
  arma::vec bicbkp = arma::zeros<arma::vec>(bic_k+1);
  for (int i=0; i < aic_k+1; ++i){
    aicbkp(i) = tempaicbkp(i)+1;
  }
  for (int i=0; i < bic_k+1; ++i){
    bicbkp(i) = tempbicbkp(i)+1;
  }
  /*organize result*/
  arma::mat aictheta = theta_list.slice(aic_k);
  arma::vec aicphi = phi_list.col(aic_k);
  arma::rowvec aicxi = xi_list.row(aic_k);
  double aicrss = rss_list(aic_k);
  double aicaic = aic_list(aic_k);
  double aicbic = bic_list(bic_k);
  arma::mat bictheta = theta_list.slice(bic_k);
  arma::vec bicphi = phi_list.col(bic_k);
  arma::rowvec bicxi = xi_list.row(bic_k);
  double bicrss = rss_list(bic_k);
  double bicaic = aic_list(bic_k);
  double bicbic = bic_list(bic_k);
  return Rcpp::List::create(
    Rcpp::Named("aic") = Rcpp::List::create(
      Rcpp::Named("bkp") = aicbkp,
      Rcpp::Named("theta") = aictheta,
      Rcpp::Named("phi") = aicphi,
      Rcpp::Named("xi") = aicxi,
      Rcpp::Named("rss") = aicrss,
      Rcpp::Named("aic") = aicaic,
      Rcpp::Named("bic") = aicbic
    ),
    Rcpp::Named("bic") = Rcpp::List::create(
      Rcpp::Named("bkp") = bicbkp,
      Rcpp::Named("theta") = bictheta,
      Rcpp::Named("phi") = bicphi,
      Rcpp::Named("xi") = bicxi,
      Rcpp::Named("rss") = bicrss,
      Rcpp::Named("aic") = bicaic,
      Rcpp::Named("bic") = bicbic
    ),
    Rcpp::Named("bkp") = bkp_list,
    Rcpp::Named("aicerror") = aic_list,
    Rcpp::Named("bicerror") = bic_list,
    Rcpp::Named("rss") = rss_list,
    Rcpp::Named("looplist") = loop_list,
    Rcpp::Named("thetalist") = theta_list,
    Rcpp::Named("philist") = phi_list,
    Rcpp::Named("xilist") = xi_list
  );
}



/////////////////OLD VERSIONS OF CODES//////////
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat multiplyXtXBySparse_c_old(int n, arma::uvec ind, arma::mat val, arma::vec phi, arma::vec w){
  int p = val.n_cols;
  int as = ind.size();
  arma::uvec indrev = arma::zeros<arma::uvec>(n-1);
  for (int i=0; i<(n-1); ++i){
    indrev(i) = n-2-i;
  }
  arma::uvec o = sort_index(ind);
  ind = sort(ind);
  val = val.rows(o);
  arma::vec phicum = cumsum(phi%phi);
  arma::mat r = val.each_col() % w(ind);
  arma::mat s1 = r.each_col() % phicum(ind);
  arma::mat cumsums1 = cumsum(s1);
  arma::rowvec s = cumsums1.row(as-1)/n;
  arma::mat matrixT = arma::zeros<arma::mat>(n-1, p);
  matrixT.rows(indrev(ind)) = r;
  matrixT = cumsum(matrixT);
  matrixT = matrixT.rows(indrev);
  arma::mat u = matrixT.each_row()-s;
  arma::vec phisq = phi%phi;
  u.each_col() %= phisq.subvec(0,n-2);
  arma::mat U = cumsum(u);
  arma::mat C = U.each_col() % w;
  return C;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat leftMultiplyByXt_c_old(arma::mat Y, arma::vec phi, arma::vec w){
  int n = Y.n_rows;
  int p = Y.n_cols;
  Y.each_col() %= phi;
  arma::mat u = cumsum(Y);
  arma::mat c = arma::zeros<arma::mat>(n-1,p);
  arma::vec phicum = cumsum(phi%phi);
  for (int i=0; i<p; ++i){
    arma::vec x = u.col(i);
    c.col(i) = w % (phicum.subvec(0,n-2) * (x(n-1)/n) - x.subvec(0,n-2));
  }
  return c;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat leftMultiplyByInvXAtXA_c_old(double n, arma::uvec ind, arma::mat val, arma::vec phi, arma::vec w){
  int a = val.n_rows;
  int p = val.n_cols;
  int as = ind.size();
  arma::mat r = arma::zeros<arma::mat>(a, p);
  if (a!=0){
    arma::vec v1 = arma::zeros<arma::vec>(2+as);
    v1(0) = 0;
    arma::vec phicum = cumsum(phi%phi);
    v1.subvec(1,as)=phicum(ind);
    v1(as+1) = n;
    arma::vec v = diff(v1);
    arma::vec d = w(ind);
    arma::mat R = arma::zeros<arma::mat>(a+2, p);
    val.each_col() /= d;
    R.row(0).fill(0);
    R.rows(1,a) = val;
    R.row(a+1).fill(0);
    arma::mat gamma = diff(R);
    arma::mat delta = gamma;
    delta.each_col() /= v;
    r = -diff(delta);
    r.each_col() /= d;
  }
  return r;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec ic_c_old(int k, arma::mat Y, arma::vec phi, arma::mat theta, double p, double n){
  arma::mat dis = arma::zeros<arma::mat>(p, n);
  double np = n*p;
  double rss = 0;
  for(int i = 0; i < p; ++i){
    for(int j = 0; j < n; ++j){
      dis(i, j) = Y(i, j) - theta(i, j)*phi(i);
      rss += dis(i, j)*dis(i, j);
    }
  }
  rss = rss/(2*np);
  const double bic_p = log(p) * ((k+1)*n + (p-k-1) + 2*k)/np;
  const double temp1 = 2*((k+1)*n+(p-k-1)+3*k);
  arma::vec out = arma::zeros<arma::vec>(3);
  out[0] = rss;
  out[1] = log(rss) + temp1/np;
  out[2] = log(rss) + bic_p;
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec defaultWeights_c_old(int p){
  arma::vec out = arma::zeros<arma::vec>(p-1);
  for(int i = 0; i < (p-1); ++i){
    double temp = (i+1) * (p-i-1);
    out[i] = (sqrt(p/temp));
  }
  return out;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List doLars_slow_old(arma::mat Y, int K, arma::vec phi, arma::vec wts, int p, int n, double epsilon = 1e-09){
  arma::vec res_lambda = arma::zeros<arma::vec>(K);
  arma::uvec res_bkp = arma::zeros<arma::uvec>(K);
  arma::mat res_value = arma::zeros<arma::mat>(K,n);
  arma::mat res_value_old = arma::zeros<arma::mat>(K,n);
  arma::mat matrixT = phi * wts.t(); //382 x 381
  arma::mat result = arma::zeros<arma::mat>(K, n+2);
  /*empty the upper triangle*/
  for (int i = 0; i < matrixT.n_cols; ++i){
    matrixT.submat(i,i,i, matrixT.n_cols-1).fill(0);
  }
  arma::mat QT1 = phi.t() * matrixT; //1x381
  arma::mat QT2 = (phi/p) * QT1; //382 x 381
  arma::mat QT = matrixT - QT2; //382 x 381
  arma::mat QY1 = phi.t() * Y; //1x32
  arma::mat QY2 = (phi/p) * QY1; //382x32
  arma::mat QY = Y - QY2; //382x32
  arma::mat c = QT.t() * QY; //381 x 32
  for (int ii=0; ii<K; ++ii){
    arma::mat cc = c%c; //381 x 32
    arma::mat ccrowsum = cumsum(cc, 1);
    arma::vec cNorm = ccrowsum.col(n-1); //381x1
    double bigcHat = cNorm.max();
    if (ii==0){
      res_bkp[ii] = cNorm.index_max();
    }
    arma::uvec AS = res_bkp.subvec(0,ii);
    arma::uvec I = sort_index(AS);
    arma::mat w1 = inv_sympd(QT.cols(AS).t() * QT.cols(AS)); //(ii+1)x(ii+1)
    arma::mat w = w1 * c.rows(AS); //(ii+1) x n
    arma::mat tempa = QT.t() * QT.cols(AS); //(p-1)x(ii+1)
    arma::mat a = tempa * w; //(p-1)xn
    arma::mat aa = a%a;
    arma::mat cumsumaa = cumsum(aa, 1);
    arma::vec rowsumsa2 = cumsumaa.col(aa.n_cols-1);
    arma::vec a1 = bigcHat - rowsumsa2;
    arma::mat u = a%c;
    arma::mat cumsumu = cumsum(u, 1);
    arma::vec rowsumsu = cumsumu.col(u.n_cols-1);
    arma::vec a2 = bigcHat-rowsumsu;
    arma::vec a3 = bigcHat-cNorm;
    arma::mat gammaTemp = arma::zeros<arma::mat>(p-1, 2);
    gammaTemp.fill(NA_REAL);
    arma::uvec subset = find(a1 > epsilon);
    arma::vec delta = arma::zeros<arma::vec>(subset.size());
    arma::uvec onetemp = arma::zeros<arma::uvec>(1);
    arma::uvec twotemp = arma::ones<arma::uvec>(1);
    delta = a2(subset)%a2(subset)-a1(subset)%a3(subset);
    arma::uvec deltaneg = find(delta<0);
    arma::uvec deltapos = find(delta>=0);
    if(deltaneg.size()>0){
      arma::uvec delta_neg = subset(deltaneg);
      gammaTemp.submat(delta_neg, onetemp).fill(NA_REAL);
      gammaTemp.submat(delta_neg, twotemp).fill(NA_REAL);
    }
    if (deltapos.size()>0){
      arma::uvec delta_pos = subset(deltapos);
      gammaTemp.submat(delta_pos, onetemp) = (a2(delta_pos) +
        sqrt(delta(deltapos)))/a1(delta_pos);
      gammaTemp.submat(delta_pos, twotemp) = (a2(delta_pos) -
        sqrt(delta(deltapos)))/a1(delta_pos);
    }
    subset = find((a1 <= epsilon) && (a2 > epsilon));
    if(subset.size()>0){
      arma::vec whattofill = a3(subset)/(2*a2(subset));
      for (int i=0; i<subset.size(); ++i){
        gammaTemp.row(subset[i]).fill(whattofill(i));
      }
    }
    double maxg = gammaTemp.max() + 1;
    subset = find((a1 <= epsilon) && (a2 <= epsilon));
    if(subset.size()>0){
      gammaTemp.rows(subset).fill(maxg);
    }
    gammaTemp.rows(AS).fill(maxg);
    gammaTemp(find(gammaTemp<=0)).fill(maxg);
    double gamma = gammaTemp.min();
    int idx = gammaTemp.index_min();
    int nexttoadd = (idx)%(p-1);
    res_lambda[ii] = sqrt(bigcHat);
    res_value.rows(0,ii) = gamma*w;
    if(ii > 0){
      res_value += res_value_old;
    }
    if (ii < (K-1)) {
      res_bkp(ii+1) = nexttoadd;
      c = c - gamma * a;
    }
    res_value_old = res_value;
  }
  arma::uvec one = arma::ones<arma::uvec>(res_bkp.size());
  res_bkp += one;
  return Rcpp::List::create(
    Rcpp::Named("bkp") = res_bkp,
    Rcpp::Named("lambda") = res_lambda,
    Rcpp::Named("value") = res_value
  );
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List doLars_c_old(arma::mat Y, int K, arma::vec phi, arma::vec wts, int p, int n, double epsilon = 1e-09){
  arma::vec res_lambda = arma::zeros<arma::vec>(K);
  arma::uvec res_bkp = arma::zeros<arma::uvec>(K);
  arma::mat res_value = arma::zeros<arma::mat>(K,n);
  arma::mat res_value_old = arma::zeros<arma::mat>(K,n);
  //arma::mat res_value;
  //arma::mat res_value_old;
  arma::mat result = arma::zeros<arma::mat>(K, n+2);
  arma::uvec AS(1);
  arma::uvec AS0(1);
  /*empty the upper triangle*/
  arma::mat c = leftMultiplyByXt_c_old(Y, phi, wts);
  for (int ii=0; ii<K; ++ii){
    arma::mat cc = c%c; //381 x 32
    arma::mat ccrowsum = cumsum(cc, 1);
    arma::vec cNorm = ccrowsum.col(n-1); //381x1
    double bigcHat = cNorm.max();
    if (ii==0){
      AS(0) = cNorm.index_max();
      res_bkp[ii] = cNorm.index_max();
    }
    //arma::uvec AS = res_bkp.subvec(0,ii);
    arma::uvec I = sort_index(AS);
    AS0 = AS;
    AS = AS(I);
    arma::mat w = leftMultiplyByInvXAtXA_c_old(p,AS,c.rows(AS),phi,wts);
    arma::mat a = multiplyXtXBySparse_c_old(p, AS, w, phi, wts);
    arma::mat aa = a%a;
    arma::mat cumsumaa = cumsum(aa, 1);
    arma::vec rowsumsa2 = cumsumaa.col(aa.n_cols-1);
    arma::vec a1 = bigcHat - rowsumsa2;
    arma::mat u = a%c;
    arma::mat cumsumu = cumsum(u, 1);
    arma::vec rowsumsu = cumsumu.col(u.n_cols-1);
    arma::vec a2 = bigcHat-rowsumsu;
    arma::vec a3 = bigcHat-cNorm;
    arma::mat gammaTemp = arma::zeros<arma::mat>(p-1, 2);
    gammaTemp.fill(NA_REAL);
    arma::uvec subset = find(a1 > epsilon);
    arma::vec delta = arma::zeros<arma::vec>(subset.size());
    arma::uvec onetemp = arma::zeros<arma::uvec>(1);
    arma::uvec twotemp = arma::ones<arma::uvec>(1);
    delta = a2(subset)%a2(subset)-a1(subset)%a3(subset);
    arma::uvec deltaneg = find(delta<0);
    arma::uvec deltapos = find(delta>=0);
    //if(deltaneg.size()>0){
    arma::uvec delta_neg = subset(deltaneg);
    gammaTemp.submat(delta_neg, onetemp).fill(NA_REAL);
    gammaTemp.submat(delta_neg, twotemp).fill(NA_REAL);
    //}
    //if (deltapos.size()>0){
    arma::uvec delta_pos = subset(deltapos);
    gammaTemp.submat(delta_pos, onetemp) = (a2(delta_pos) +
      sqrt(delta(deltapos)))/a1(delta_pos);
    gammaTemp.submat(delta_pos, twotemp) = (a2(delta_pos) -
      sqrt(delta(deltapos)))/a1(delta_pos);
    //}
    subset = find((a1 <= epsilon) && (a2 > epsilon));
    //if(subset.size()>0){
    arma::vec whattofill = a3(subset)/(2*a2(subset));
    for (int i=0; i<subset.size(); ++i){
      gammaTemp.row(subset[i]).fill(whattofill(i));
    }
    //}
    double maxg = gammaTemp.max() + 1;
    subset = find((a1 <= epsilon) && (a2 <= epsilon));
    //if(subset.size()>0){
    gammaTemp.rows(subset).fill(maxg);
    //}
    gammaTemp.rows(AS).fill(maxg);
    gammaTemp(find(gammaTemp<=0)).fill(maxg);
    double gamma = gammaTemp.min();
    int idx = gammaTemp.index_min();
    arma::uvec nexttoadd(1);
    nexttoadd(0) = (idx)%(p-1);
    double nexttoadd_double = nexttoadd(0);
    res_lambda[ii] = sqrt(bigcHat);
    res_value.fill(0);
    res_value.rows(I) = gamma*w;
    if(ii > 0){
      res_value += res_value_old;
    }
    if (ii < (K-1)) {
      AS = join_cols(AS0, nexttoadd);
      res_bkp(ii+1) = nexttoadd_double;
      c = c - gamma * a;
    }
    res_value_old = res_value;
  }
  arma::uvec one = arma::ones<arma::uvec>(res_bkp.size());
  res_bkp += one;
  return Rcpp::List::create(
    Rcpp::Named("bkp") = res_bkp,
    Rcpp::Named("lambda") = res_lambda,
    Rcpp::Named("value") = res_value
  );
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List cnv_c_old(arma::mat Y, arma::vec wts, int steps, int maxloop=10){
  /* Initialize phi*/
  const int p = Y.n_rows;/*number of probes*/
  const int n = Y.n_cols; /*number of samples*/
  arma::vec phi = arma::ones<arma::vec>(p);
  int loop = 0;
  arma::vec phi_old = phi + 1e+10;
  double error = 1e+10;
  //double error2 = 1e+10;
  while (error > 1e-5 & loop < 10){
    arma::mat theta1 = (phi.t() * Y / p).t();
    phi = Y * theta1 / accu(theta1 % theta1);
    phi = phi / sqrt(accu(phi%phi)/p);
    arma::vec diff = phi_old - phi;
    error = sqrt(accu(diff%diff)/p);
    phi_old = phi;
    loop += 1;
  }
  /*Initialize variables to save*/
  arma::vec aic_list = arma::zeros<arma::vec>(steps-1);
  arma::vec bic_list = arma::zeros<arma::vec>(steps-1);
  arma::vec rss_list = arma::zeros<arma::vec>(steps-1);
  arma::vec loop_list = arma::zeros<arma::vec>(steps-1);
  arma::mat bkp_list = arma::zeros<arma::mat>(steps-1, steps-1);
  arma::mat phi_list = arma::zeros<arma::mat>(p, steps-1);
  arma::cube theta_list = arma::zeros<arma::cube>(p, n, steps-1);

  /* start computing the result for each k*/
  for (int k = 1; k < (steps); ++k){
    /* For each given k, we compute bkp, aic, etc*/
    loop=0; error=1e5;
    /* Initialize variable for next loop */
    arma::mat theta = arma::zeros<arma::mat>(p, n);
    arma::vec bkp = arma::zeros<arma::vec>(k);
    arma::vec lambda = arma::zeros<arma::vec>(k);
    arma::mat value = arma::zeros<arma::mat>(k, n);
    arma::vec sorted_bkp = arma::zeros<arma::vec>(k);
    arma::vec phi_new = arma::zeros<arma::vec>(p);
    arma::vec q = arma::zeros<arma::vec>(maxloop);
    bool alternating = TRUE;
    while (error > 1e-5 & loop<maxloop & alternating){
      Rcpp::List lars = doLars_c_old(Y, k, phi, wts, p, n);
      bkp = as<arma::vec>(lars["bkp"]);
      lambda = as<arma::vec>(lars["lambda"]);
      value = as<arma::mat>(lars["value"]);
      arma::uvec ord = sort_index(bkp);
      sorted_bkp = sort(bkp);
      arma::uvec sorted(k);
      for (int i=0; i<k; ++i){
        sorted(i) = sorted_bkp(i)-1;
      }
      arma::mat delta2 = arma::zeros<arma::mat>(p-1, n);
      delta2.rows(sorted) = value.rows(ord);
      arma::mat delta3 = arma::zeros<arma::mat>(p-1, n);
      delta3 = delta2.each_col() % wts;
      arma::mat delta1 = phi.t() * Y;
      arma::vec phisq = phi%phi;
      arma::vec cumsumphi = flipud(cumsum(flipud(phisq)));
      cumsumphi.shed_row(0);
      delta1 = (delta1-(cumsumphi.t() * delta3))/accu(phisq);
      arma::mat delta = join_cols(delta1, delta3);
      theta = delta;
      for (int i=1; i<p; ++i){
        theta.row(i) = theta.row(i-1) + delta.row(i);
      }
      arma::mat partition(k+1, p);
      partition.fill(NA_REAL);
      if (k==1){
        for (int i=0; i < bkp(0); ++i){
          partition(0, i) = i+1;
        }
        for (int i=bkp(0); i < p; ++i){
          partition(1, i-bkp(0)) = i+1;
        }
      }
      else{
        for (int i=0; i < sorted(0)+1; ++i){
          partition(0, i) = i+1;
        }
        for (int i=sorted(k-1)+1; i<p; ++i){
          partition(k, i-sorted(k-1)-1) = i+1;
        }
        for (int i=0; i < (k-1); ++i){
          for (int j=sorted(i)+1; j<sorted(i+1)+1; ++j){
            partition(i+1, j-sorted(i)-1) = j+1;
          }
        }
      }
      arma::mat thetaY = theta%Y;
      arma::mat thetaYcumsum = cumsum(thetaY, 1);
      arma::mat thetathetacumsum = cumsum(theta%theta, 1);
      arma::vec thetaYrowsum = thetaYcumsum.col(n-1);
      arma::vec thetasqrowsum = thetathetacumsum.col(n-1);
      phi_new = thetaYrowsum / thetasqrowsum;
      for (int i=0; i<(k+1); ++i){
        arma::rowvec ith_partition = partition.row(i);
        arma::vec where = ith_partition(find_finite(ith_partition));
        arma::uvec uwhere = arma::zeros<arma::uvec>(where.size());
        for (int j = 0; j < where.size(); ++j){
          uwhere(j) = where(j)-1;
        }
        arma::vec phi_temp = phi_new(uwhere);
        double phitempsqmean = mean(phi_temp%phi_temp);
        phi_new(uwhere) = phi_temp/sqrt(phitempsqmean);
      }
      arma::vec diff = phi-phi_new;
      error = accu(diff%diff) / p;
      error = sqrt(error);
      phi = phi_new;
      arma::vec res = ic_c_old(k,Y,phi,theta,p,n);
      loop += 1;
      //Compute the optimization objective Q1, check convergence
      double q1 = res(0) * 2* n * p;
      for (int i = 0; i < (p-1); ++i){
        arma::rowvec temp = theta.row(i+1)-theta.row(i);
        q1 += lambda(k-1) * sqrt(accu(temp%temp)) / wts(i);
      }
      q(loop-1) = q1;
      if(loop>5){
        alternating = (abs(q(loop-1)-q(loop-3)) + abs(q(loop-2) - q(loop-4)) > 1e-2);
      }
    }
    arma::vec res = ic_c_old(k, Y, phi, theta, p, n);
    rss_list(k-1) = res(0);
    aic_list(k-1) = res(1);
    bic_list(k-1) = res(2);
    loop_list(k-1) = loop;
    for (int i=0; i<k; ++i){
      bkp_list(k-1, i) = sorted_bkp(i);
    }
    theta_list.slice(k-1) = theta;
    phi_list.col(k-1) = phi;
  }
  int aic_k = aic_list.index_min();
  int bic_k = bic_list.index_min();

  /*Organize result*/
  arma::rowvec tempaicbkp = bkp_list.row(aic_k);
  arma::rowvec tempbicbkp = bkp_list.row(bic_k);
  arma::vec aicbkp = arma::zeros<arma::vec>(aic_k+1);
  arma::vec bicbkp = arma::zeros<arma::vec>(bic_k+1);
  for (int i=0; i < aic_k+1; ++i){
    aicbkp(i) = tempaicbkp(i)+1;
  }
  for (int i=0; i < bic_k+1; ++i){
    bicbkp(i) = tempbicbkp(i)+1;
  }
  /*organize result*/
  arma::mat aictheta = theta_list.slice(aic_k);
  arma::vec aicphi = phi_list.col(aic_k);
  double aicrss = rss_list(aic_k);
  double aicaic = aic_list(aic_k);
  double aicbic = bic_list(bic_k);
  arma::mat bictheta = theta_list.slice(bic_k);
  arma::vec bicphi = phi_list.col(bic_k);
  double bicrss = rss_list(bic_k);
  double bicaic = aic_list(bic_k);
  double bicbic = bic_list(bic_k);
  return Rcpp::List::create(
    Rcpp::Named("aic") = Rcpp::List::create(
      Rcpp::Named("bkp") = aicbkp,
      Rcpp::Named("theta") = aictheta,
      Rcpp::Named("phi") = aicphi,
      Rcpp::Named("rss") = aicrss,
      Rcpp::Named("aic") = aicaic,
      Rcpp::Named("bic") = aicbic
    ),
    Rcpp::Named("bic") = Rcpp::List::create(
      Rcpp::Named("bkp") = bicbkp,
      Rcpp::Named("theta") = bictheta,
      Rcpp::Named("phi") = bicphi,
      Rcpp::Named("rss") = bicrss,
      Rcpp::Named("aic") = bicaic,
      Rcpp::Named("bic") = bicbic
    ),
    Rcpp::Named("bkp") = bkp_list,
    Rcpp::Named("aicerror") = aic_list,
    Rcpp::Named("bicerror") = bic_list,
    Rcpp::Named("rss") = rss_list,
    Rcpp::Named("looplist") = loop_list
  );
}



