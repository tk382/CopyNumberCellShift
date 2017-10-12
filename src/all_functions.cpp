#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;

/*
 * We look at the model
 * Y_ij - \phi_i (\Theta_{ij} + \delta_j)
 */


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
  //arma::mat res_value;
  //arma::mat res_value_old;
  arma::mat result = arma::zeros<arma::mat>(K, n+2);
  arma::uvec AS(1);
  arma::uvec AS0(1);
  /*empty the upper triangle*/
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
    //arma::uvec AS = res_bkp.subvec(0,ii);
    arma::uvec I = sort_index(AS);
    AS0 = AS;
    AS = AS(I);
    arma::mat w = leftMultiplyByInvXAtXA_c(p,AS,c.rows(AS),phi,wts);
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
Rcpp::List cnv_c(arma::mat Y, arma::vec wts, int steps, int maxloop=10){
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
      phi = phi_new;

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
    Rcpp::Named("looplist") = loop_list
  );
}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_f(arma::mat Y, arma::vec sigma, arma::vec delta){
  double f = 0;
  int p = Y.n_rows;
  int n = Y.n_cols;
  arma::mat Ysq = Y%Y;
  arma::vec sigmasq = sigma%sigma;
  arma::vec deltasq = delta%delta;
  for (int i=0; i < p; ++i){
    for (int j=0; j < n; ++j){
      f -= log(sigmasq(i) + deltasq(j))/2;
      f -= Ysq(i,j)/(2*(sigmasq(i) + deltasq(j)));
    }
  }
  return f;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec get_g(arma::mat Y, arma::vec sigma, arma::vec delta2){
  int p = Y.n_rows;
  int n = Y.n_cols;
  arma::rowvec delta = arma::zeros<arma::rowvec>(n);
  for (int i=0; i<n; ++i){
    delta(i) = delta2(i);
  }
  arma::mat Ysq = Y%Y;
  arma::mat gsig = arma::zeros<arma::mat>(p,n);
  arma::mat gdel = arma::zeros<arma::mat>(p,n);
  arma::mat sigmat = arma::zeros<arma::mat>(p,n);
  arma::mat delmat = arma::zeros<arma::mat>(p,n);
  sigmat.each_col()+=sigma;
  delmat.each_row()+=delta;
  gsig = -sigmat%sigmat%sigmat + sigmat%(Ysq-delmat%delmat);
  arma::mat denom = sigmat%sigmat + delmat%delmat;
  denom = denom%denom;
  gsig = gsig / denom;
  gdel = -delmat%delmat%delmat + delmat%(Ysq-sigmat%sigmat);
  gdel = gdel / denom;
  arma::mat cumgsig = cumsum(gsig, 1);
  arma::mat cumgdel = cumsum(gdel.t(), 1);
  arma::vec g = arma::zeros<arma::vec>(n+p);
  g.subvec(0,p-1) = cumgsig.col(n-1);
  g.subvec(p, n+p-1) = cumgdel.col(p-1);
  return g;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List get_f_and_g(arma::vec sigma, arma::vec delta, arma::mat Y){
  double f = 0;
  int p = Y.n_rows;
  int n = Y.n_cols;
  arma::mat Ysq = Y%Y;
  arma::vec sigmasq = sigma%sigma;
  arma::vec deltasq = delta%delta;
  arma::mat gs = arma::zeros<arma::mat>(p,n);
  arma::mat gd = arma::zeros<arma::mat>(n,p);
  for (int i=0; i < p; ++i){
    for (int j=0; j < n; ++j){
      f -= log(sigmasq(i) + deltasq(j))/2;
      f-= Ysq(i,j)/(2*(sigmasq(i) + deltasq(j)));
      gs(i,j) = -pow(sigma(i), 3) + sigma(i)*(Ysq(i,j)-deltasq(j));
      gs(i,j) /= pow(sigmasq(i) + deltasq(j), 2);
      gd(j,i) = -pow(delta(j),3) + delta(j)*(Ysq(i,j)-sigmasq(i));
      gd(j,i) /= pow(sigmasq(i) + deltasq(j), 2);
    }
  }
  arma::vec g = arma::zeros<arma::vec>(p+n);
  arma::mat cumsumgs = cumsum(gs,1);
  g.subvec(0, p-1) = cumsumgs.col(n-1);
  arma::mat cumsumgd = cumsum(gd, 1);
  g.subvec(p, p+n-1) = cumsumgd.col(p-1);

  return Rcpp::List::create(Rcpp::Named("f") = f, Rcpp::Named("g") = g);
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double get_alpha(arma::mat Y, double alpha, double rho, double c, arma::vec sigma, arma::vec delta){
  int p = Y.n_rows;
  int n = Y.n_cols;
  double f1 = get_f(Y, sigma, delta);
  arma::vec g1 = get_g(Y, sigma, delta);
  double g1cross = accu(g1%g1);
  int maxit = 1000;
  int it = 1;
  double f2 = f1-100;
  while(f2 > f1 - c*alpha*g1cross && it < maxit){
    arma::vec sigma2 = sigma - alpha*g1.subvec(0,p-1);
    arma::vec delta2 = delta - alpha*g1.subvec(p, p+n-1);
    f2 = get_f(Y, sigma2, delta2);
    alpha = alpha * rho;
    it = it + 1;
  }
  return alpha;
}

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List grdesc(arma::mat Y, arma::vec sig, arma::vec del, int maxit){
  int p = Y.n_rows;
  int n = Y.n_cols;
  arma::vec newsig = sig;
  arma::vec newdel = del;
  arma::vec g = arma::zeros<arma::vec>(n+p);
  double alpha0 = 1;
  double rho = 0.9;
  double c = 0.1;
  double diff = 1;
  int it = 1;
  double alpha = 10;
  while(it < maxit && diff > 1e-4){
    alpha = get_alpha(Y, alpha0, rho, c, sig, del);
    g = get_g(Y, sig, del);
    newsig = sig - alpha * g.subvec(0,p-1);
    newdel = del - alpha * g.subvec(p, n+p-1);
    diff = accu((newsig-sig)%(newsig-sig))+accu((newdel-del)%(newdel-del));
    sig = newsig;
    del = newdel;
    it += 1;
  }
  return Rcpp::List::create(
    Rcpp::Named("sigma") = newsig, Rcpp::Named("delta") = newdel
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





