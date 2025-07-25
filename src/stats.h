/*************************
 Statistical tests and other useful tools
**************************/
// Regression models
struct REGmodel{
    // INPUTS:
    vec y;              // output data
    mat X;              // input data
    // OUTPUTS
    vec beta, e, stdBeta;
    double AIC, BIC, AICc;
    mat covBeta;
};
/***************************************************
* Function declarations
****************************************************/
// Autocorrelation function
void acf(vec&, int, vec&);
// Regression
void regression(vec, mat, REGmodel&);
void regress(vec, mat, vec&, vec&, vec&, double&, double&, double&);
// Augmented Dickey-Fuller test
double adfTest(vec&, vec, double&, double&, double&);
// Augmented Dickey-Fuller test using criterion identification
int adfTests(vec, double, string);
// Harmonic regression
void harmonicRegress(vec&, mat&, vec, uword, vec&, vec&, vec&, vec&);
// Select harmonics of seasonal component
void selectHarmonics(vec&, mat&, vec, uvec&, vec&, string&);
// Identify AR model via information criterion
void selectAR(vec&, double, string, vec&, vec&, double&, vec&);
// Identify stationary ARMA model
void selectARMA(vec, double, int, string, vec&);
// Estimation of ARMA model by Least Squares (Hannan-Rissanen)
void linearARMA(vec&, vec, vec&, vec&);
// Information criteria
void infoCriteria(double, int, int, double&, double&, double&);
// Gaussianity Bera-Jarque
void beraj(vec&, double&, double&);
// Heteroskedasticity ratio of variances
void heterosk(vec&, double&, double&, int&);
// Binomial CDF calculation
double binoCdf(double, double, double);
vec binoCdf(double, double, vec);
// F CDF calculation (only valid for d1 and d2 even)
template <class T>
T fCdf(T, int, int);
// Student t CDF calculation
double tCdf(double, double);
vec tCdf(vec, double);
// Mean of vector or matrix (in cols) with nan or inf values
double nanMean(vec y);
rowvec nanMean(mat y);
double nanStddev(vec y);
rowvec nanStddev(mat y);
double nanMin(vec y);
// Incomplete beta function 
double betaInc(double, double, double);
// Estimation table
void outputTable(vec, vector<string>&);
// Print vec with format
// void printVec(rowvec, const string, const string, char&);

// Standard deviation of vector or matrix (in cols) with nan or inf values
/***************************************************
 * Function implementations
 ****************************************************/
// Autocorrelation function
void acf(vec& y, int ncoef, vec& acfCoef){
  int nNan;
  vec yn = removeNans(y, nNan);
  yn = (y - mean(yn)) / stddev(yn);
  int nNo, n = yn.n_elem;
  vec prod;
  acfCoef.zeros(ncoef);
  for (int i = 0; i < ncoef; i++){
    prod = yn(span(i + 1, n - 1)) % yn(span(0, n - i - 2));
    nNo = n - 2 * nNan - i - 1;
    acfCoef(i) = sum(prod(find_finite(prod))) / nNo;
  }
}
// Regression
void regression(vec y, mat X, REGmodel& m){
  vec beta, stdBeta, e;
  double AIC, BIC, AICc;
//  mat covBeta;
  regress(y, X, beta, stdBeta, e, BIC, AIC, AICc);
  m.beta = beta;
//  m.covBeta = covBeta;
//  m.stdBeta = sqrt(covBeta.diag());
  m.stdBeta = stdBeta;
  m.X = X;
  m.e = e;
  m.BIC = BIC;
  m.AIC = AIC;
  m.AICc = AICc;
}
void regress(vec y, mat X, vec& beta, vec& stdBeta, vec& eOut, double& BIC, double& AIC, double& AICc){
  eOut = y;
  uvec ind = find_finite(mean(join_rows(y, X), 1));
  X = X.rows(ind);
  y = y.rows(ind);
  int k = X.n_cols;
  mat iX = pinv(X.t() * X), covBeta;
  beta = iX * (X.t() * y);
  vec e = y - X * beta;
  int n = e.n_elem;
  vec varE = (e.t() * e) / (n - k);
  covBeta = varE(0) * iX;
  stdBeta = sqrt(covBeta.diag());
  vec nlv = log(varE * (n - k) / n);
  eOut(ind) = e;
  AIC = nlv(0) + 2 * k / n;
  BIC = nlv(0) + k * log(n) / n;
  if (n - k - 1 != 0)
        AICc = (AIC * n + (2 * k * (1 + k)) / (n - k - 1)) / n;
}
// Augmented Dickey-Fuller test
double adfTest(vec& y, vec lags, double& BIC, double& AIC, double& AICc){
  int n = y.n_elem, maxLag = max(lags);
  // Removing mean
  int aux;
  y = y - mean(removeNans(y, aux));
  // Creating lagged matrix
  mat X = join_rows(join_rows(y(span(maxLag + 1, n - 1)), y(span(maxLag, n - 2))), 
                    lag(y(span(1, n - 1)) - y(span(0, n - 2)), lags));
  // Regression with constant and no trend
  vec beta, stdBeta, e;
  mat covBeta;
  regress(X.col(0), X.cols(span(1, maxLag + 1)), beta, stdBeta, e, BIC, AIC, AICc);
//  stdBeta = sqrt(covBeta.diag());
  return (beta(0) - 1) / stdBeta(0);
}
// Augmented Dickey-Fuller test using criterion identification
int adfTests(vec y, double nModels, string criterion){
  double BIC, AIC, AICc;
  if (nModels > y.n_elem / 3){
      nModels = floor(y.n_elem / 3);
  }
  vec vCrit(nModels);
  vec unitRoot(nModels);
  for (unsigned int i = 0; i < nModels; i++){
    unitRoot(i) = adfTest(y, regspace(1, i + 1), BIC, AIC, AICc);
    if (criterion[0] == 'b'){
      vCrit(i) = BIC;
    } else if (criterion == "aic"){
      vCrit(i) = AIC;
    } else {
      vCrit(i) = AICc;
    }
  }
  if (unitRoot(vCrit.index_min()) > -2){    // with trend
    return(1);
  } else if (unitRoot(vCrit.index_min()) < -5) {  // no trend
    return(0);
  } else {                                  // Not sure
    return(-1);
  }
}
// Harmonic regression
void harmonicRegress(vec& y, mat& u, vec period, uword trendPow, vec& beta, vec& stdBeta, vec& e, vec& yFit){
    // It includes a cubic trend, a constant, harmonics and inputs
    int n = y.n_elem, k = u.n_rows, pos;
    // if (season < 2) season = 4;
    // vec line = regspace(1, 1, season / 2);
    // vec period = season / line;
    vec t = regspace(1, n);
    uvec harm;
    rowvec w;
    bool minus = false;
    int nHarm = 0;
    bool seasonal = true;
    if (any(period == 1)){
        seasonal = false;
    } 
    if (seasonal){
        harm = regspace<uvec>(1, 1, period.n_elem);
        w = 2 * datum::pi / period.t();
        // Checking for Nyquist frequency
        if (any(w.row(0) == datum::pi)){
            minus = true;
        }
        nHarm = harm.n_elem;
    }
    // Regressors (harmonics + constant + trend + regressors)
    mat X(n, nHarm * 2 - minus + 1 + trendPow + k);
    if (seasonal){
        // Setting cos/sin regressors
        X.cols(span(0, nHarm - 1)) = kron(w, t);
        X.cols(span(nHarm, 2 * nHarm - 1 - minus)) = sin(X.cols(span(0, nHarm - 1 - minus)));
        X.cols(span(0, nHarm - 1)) = cos(X.cols(span(0, nHarm - 1)));
    }
    pos = 2 * nHarm - minus;
    // Adding trend
    t = t / n;
    for (uword i = 0; i <= trendPow; i++)
        X.col(pos + i) = pow(t, i);
    pos += trendPow + 1;
    // Exogenous inputs
    if (k > 0){
        X.cols(span(pos, pos + k - 1)) = u.submat(0, 0, k - 1, n - 1).t();
    }
    // Regression
    double BIC, AIC, AICc;
    mat covBeta;
    regress(y, X, beta, stdBeta, e, BIC, AIC, AICc);
//    stdBeta = sqrt(covBeta.diag());
    AIC *= 1;
    BIC *= 1;
    AICc *= 1;
    yFit = X * beta;
}
// Forecasting with harmonic regression
//void forHarmonicReg(vec y, mat u, int h, vec period, uword trendPow, vec& beta,
//                    vec& stdBeta, vec& e, vec& yFor){
//    if (u.n_cols > y.n_elem){
//        uf = u.cols(y.n_elem, u.n_cols - 1);
//        u = u.cols(0, y.n_elem - 1);
//    }
//    harmonicRegress(y, u, period, trendPow, beta, stdBeta, e);
//}
// Select harmonics of seasonal component
void selectHarmonics(vec& y, mat& u, vec period, uvec& harmonics, vec& beta, string& isSeasonal){
  vec stdBeta, e, yFit;
  harmonicRegress(y, u, period, 3, beta, stdBeta, e, yFit);
  vec t = abs(beta) / stdBeta;
  int nHarm = ceil((beta.n_rows - u.n_rows - 4.0) / 2.0);
  uvec aux1 = regspace<uvec>(2 * nHarm - 1, t.n_elem - 1);
  t(aux1).fill(0.0);
  uvec ind = (t > 1.6449);   // 10% confidence interval
  uvec aux = (ind(span(0, nHarm - 1)) + ind(span(nHarm, 2 * nHarm - 1)) > 0);
  uvec harm = regspace<uvec>(1, 1, nHarm);
  harmonics = harm(find(aux)) - 1;
  isSeasonal = "dubious";
  if (any(t > 3)){
    isSeasonal = "true";
  } else if (all(t <= 1.6449)){
    isSeasonal = "false";
  }
}
// Identify AR model via information criterion
void selectAR(vec& y, double maxAR, string criterion, vec& arOrder, vec& eBest, double& minCritAR, vec& arBeta){
  // Removing mean
  int aux;
  y = y - mean(removeNans(y, aux));
  // AR identification
  minCritAR = 1e12;
  vec crit(maxAR + 1);
  mat X = lag(y, regspace(0, maxAR));
  vec beta, stdBeta, e; //, arBeta;
  double BIC, AIC, AICc; // minCritAR = 1e12; //, minCrit = 1e12;
  mat covBeta;
  arOrder.zeros(1);
  for (int i = 0; i <= maxAR; i++){
    if (i == 0){   // AR(0)
      eBest = e = X.col(0);
      AIC = BIC = AICc = log(as_scalar(eBest.t() * eBest) / eBest.n_elem);
    } else {      // AR(1)...
      regress(X.col(0), X.cols(span(1, i)), beta, stdBeta, e, BIC, AIC, AICc);
//      stdBeta = sqrt(covBeta.diag());
    }
    if (criterion == "aic"){
      crit(i) = AIC;
    } else if (criterion == "bic"){
      crit(i) = BIC;
    } else {
      crit(i) = AICc;
    }
    if (crit(i) < minCritAR){
      eBest = e;
      minCritAR = crit(i);
      arBeta = -beta;
      arOrder = i;
    }
  }
}
// Identify stationary ARMA model
void selectARMA(vec y, double period, double maxAR, string criterion, vec& orders, vec& betaOpt){
  vec arOrder;
  double minCritAR;
  vec eBest, arBeta, maxOrders(2); 
  maxOrders.fill((int)period - 1);
  orders.zeros(2);
  betaOpt.fill(0);
  selectAR(y, maxAR, criterion, arOrder, eBest, minCritAR, arBeta);
  if (arOrder(0) == 0){
    return ;
  }
  if (arOrder(0) < 3){
    orders(0) = arOrder(0);
    betaOpt = -arBeta;
    return ;
  }  
  // AR identification
  mat X = lag(y, regspace(0, maxAR));
  vec stdBeta, e;
  mat covBeta;
  double BIC, AIC, AICc, minCrit = 1e12, curCrit1 = minCrit, curCrit2 = minCrit;
  vec beta1, beta2;
  vec ind(2); ind(0) = maxAR; ind(1) = y.n_elem - maxAR - 1;
  bool nonStop = true;
  mat eData = lag(eBest, regspace(1, maxAR));
  do {
    if (orders(0) < maxOrders(0)){
      if (orders(1) == 0){     // Pure AR(p+1)
        regress(X(span(ind(0), ind(1)), 0), X(span(ind(0), ind(1)), span(1, orders(0) + 1)), 
                  beta1, stdBeta, e, BIC, AIC, AICc);
      } else {                // ARMA(p+1, q)
        regress(X(span(ind(0), ind(1)), 0), join_rows(X(span(ind(0), ind(1)), span(1, orders(0) + 1)), 
                  eData.cols(span(0, orders(1) - 1))), beta1, stdBeta, e, BIC, AIC, AICc);
      }
      if (criterion == "aic"){
        curCrit1 = AIC;
      } else if (criterion == "bic"){
        curCrit1 = BIC;
      } else {
        curCrit1 = AICc;
      }
//      stdBeta = sqrt(covBeta.diag());
    }
    if (orders(1) < maxOrders(1)){    
      if (orders(0) == 0){           // Pure MA(q+1)
        regress(X(span(ind(0), ind(1)), 0), eData.cols(span(0, orders(1))), beta2, stdBeta, e, BIC, AIC, AICc);
      } else {                      // ARMA(p, q+1)
        regress(X(span(ind(0), ind(1)), 0), join_rows(X(span(ind(0), ind(1)), span(1, orders(0))), 
                  eData.cols(span(0, orders(1)))), beta2, stdBeta, e, BIC, AIC, AICc);
      }
      if (criterion == "aic"){
        curCrit2 = AIC;
      } else if (criterion == "bic"){
        curCrit2 = BIC;
      } else {
        curCrit2 = AICc;
      }
//      stdBeta = sqrt(covBeta.diag());
    }
    // Decide which model is best in iteration
    if (curCrit1 <= curCrit2 && curCrit1 < minCrit){      // Incrementing AR order best
      orders(0) += 1;
      minCrit = curCrit1;
      betaOpt = beta1;
    } else if(curCrit2 < curCrit1 && curCrit2 < minCrit){  // Incrementing MA order best
      orders(1) += 1;
      minCrit = curCrit2;
      betaOpt = beta2;
    } else {
      nonStop = false;
    }
  } while (nonStop);
  // Correction in case pure AR is selected
  // if (minCrit > minCritAR){
  //   orders(0) = arBeta.n_elem;
  //   orders(1) = 0;
  //   betaOpt = -arBeta;
  // }
  if (orders(0))
    betaOpt(span(0, orders(0) - 1)) = -betaOpt(span(0, orders(0) - 1));
}
// Estimation of ARMA model by Least Squares (Hannan-Rissanen)
void linearARMA(vec& y, vec orders, vec& beta, vec& stdBeta){
  vec aux(2), e;
  int aux1;
  y = y - mean(removeNans(y, aux1));
  // Select best AR
  int maxAR;
  vec arOrder, arBeta;
  double minCritAR;
  aux.resize(2);
  aux(0) = max(orders) * 4;
  aux(1) = 12;
  maxAR = max(aux);
  selectAR(y, maxAR, "bic", arOrder, e, minCritAR, arBeta);
  double BIC, AIC, AICc;
  vec ind(2);
  mat X = lag(y, regspace(0, orders(0)));
  mat covBeta, eData = lag(e, regspace(1, orders(1)));
  int dim;
  aux(0) = X.n_rows;
  aux(1) = eData.n_rows;
  dim = min(aux);
  if (orders(0) == 0){           // pure MA
    regress(X(span(aux(0) - dim, aux(0) - 1), span(0)),
            eData(span(aux(1) - dim, aux(1) - 1), span(0, orders(1) - 1)), beta, stdBeta, e, BIC, AIC, AICc);
  } else if (orders(1) == 0){   // pure AR
    regress(X(span(aux(0) - dim, aux(0) - 1), span(0)), X(span(aux(0) - dim, aux(0) - 1), span(1, orders(0))), 
            beta, stdBeta, e, BIC, AIC, AICc);
  } else {                      // mixed ARMA
    regress(X(span(aux(0) - dim, aux(0) - 1), span(0)), join_rows(X(span(aux(0) - dim, aux(0) - 1), span(1, orders(0))), 
              eData(span(aux(1) - dim, aux(1) - 1), span(0, orders(1) - 1))), beta, stdBeta, e, BIC, AIC, AICc);
  }
//  stdBeta = sqrt(covBeta.diag());
  if (orders(0) > 0)
    beta(span(0, orders(0) - 1)) = -beta(span(0, orders(0) - 1));
}
// Information criteria
void infoCriteria(double llik, int k, int n, double& AIC, double& BIC, double&AICc){
  AIC = -2 * (llik - k) / n;
  BIC = (-2 * llik + k * log(n)) / n;
  if (n - k - 1 > 0)
      AICc = (AIC * n + (2 * k * (1 + k)) / (n - k - 1)) / n;
  else
      AICc = datum::nan;
}
// Gaussianity Bera-Jarque
void beraj(vec& y, double& bj, double& pbj){
    int nNan, n = y.n_elem;
    vec yn = removeNans(y, nNan);
    yn = yn - mean(yn);
    vec media3; media3 = mean(pow(yn, 3));
    vec media4; media4 = mean(pow(yn, 4));
    vec stdb; stdb = sqrt(mean(pow(yn, 2)));
    vec skew = media3 / pow(stdb, 3);
    vec kurto = media4 / pow(stdb, 4) - 3;
    bj = (n - nNan) / 6 * (pow(as_scalar(skew), 2) + pow(as_scalar(kurto), 2) / 4);
    pbj = exp(-bj / 2);
}
// Heteroskedasticity ratio of variances
void heterosk(vec& y, double& F, double& pF, int& df){
  int nNan;
  vec yn = removeNans(y, nNan);
  int n = yn.n_elem;
  int i1 = n / 3, i2 = 2 * n / 3;
  if (remainder(i1, 2) == 0){
    i1 += 1;
    i2 -= 1;
  }
  df = i1 + 1;
  double cov1 = var(yn(span(0, i1)));
  double cov2 = var(yn(span(i2, n - 1)));
  if (cov1 > cov2){
    F = cov2 / cov1;
  } else {
    F= cov1 / cov2;
  }
  pF = 2 * fCdf(F, min(df, 170), min(df, 170));
}

// Binomial CDF calculation
double binoCdf(double k, double n, double p){
  vec pV(1);
  pV(0) = p;
  vec out = binoCdf(k, n, pV);
  return out(0);
}
vec binoCdf(double k, double n, vec p){
  vec pCdf(p.n_elem);
  pCdf.fill(0);
  if (k >= n){
    pCdf.fill(1);
  } else {
    double tagammanp1 = tgamma(n + 1);
    for (int i = 0; i <= k; i++){
      pCdf += (tagammanp1) / (tgamma(i + 1) * tgamma(n - i + 1)) * (pow(p, i) % pow(1 - p, n - i));
    }
  }
  return pCdf;
}
// F CDF calculation (only valid for d1 and d2 even)
template <class T>
T fCdf(T x, int d1, int d2){
  T p = d1 * x / (d1 * x + d2);
  double k = d1 / 2 - 1;
  double n = d2 / 2 + k;
  T pCdf = 1 - binoCdf(k, n, p);
  return pCdf;
}
// Student t CDF calculation
double tCdf(double t, double v) {
  double x = (t + sqrt(t * t + v)) / (2.0 * sqrt(t * t + v));
  return betaInc(v / 2.0, v / 2.0, x);
}
vec tCdf(vec t, double v) {
  int n = t.n_elem;
  vec prob(n);
  for (int i = 0; i < n; i++){
    prob(i) = tCdf(t(i), v);
  }
  return prob;
}
// Mean of vector or matrix (in cols) with nan or inf values
double nanMean(vec y){
  int nNan;
  vec x = removeNans(y, nNan);
  if (x.n_elem == 0){
      return datum::nan;
  } else {
    return mean(x);
  }
}
rowvec nanMean(mat y){
  if (y.has_nan() || y.has_inf()){  // with Nans
    int nNan;
    rowvec my = mean(y);
    uvec ind = find_nonfinite(my);
    int m = ind.n_elem;
    for (int i = 0; i < m; i++){
      vec yi = y.col(ind(i));
      my(ind(i)) = mean(removeNans(yi, nNan));
    }
    return my;
  } else {              // No nans
    return mean(y);
  }
}
// Standard deviation of vector or matrix (in cols) with nan or inf values
double nanStddev(vec y){
  int nNan;
  return stddev(removeNans(y, nNan));
}
rowvec nanStddev(mat y){
  if (y.has_nan() || y.has_inf()){  // with Nans
    int nNan;
    rowvec my = stddev(y);
    uvec ind = find_nonfinite(my);
    int m = ind.n_elem;
    for (int i = 0; i < m; i++){
      vec yi = y.col(ind(i));
      my(ind(i)) = stddev(removeNans(yi, nNan));
    }
    return my;
  } else {              // No nans
    return stddev(y);
  }
}
// Min of vector with nan or inf values
double nanMin(vec y){
  int nNan;
  return min(removeNans(y, nNan));
}
// Incomplete beta function
double betaInc(double a, double b, double x){
  // Abramowitz and Stegun, Handbook of Mathematical Functions
  // Press, WH and Teukolsky, SA (1988), Evaluating Continued 
  //        Fractions and Computing Exponential Integrals
  //         https://doi.org/10.1063/1.4822777
  // Lentz's algorithm for 0 < x < 1
  if (x < 0 || x > 1){
    return(datum::nan);
  }
  if (x > (a + 1) / (a + b + 2) || x > 1 - (b + 1) / (a + b + 2)){
    return(1 - betaInc(b, a , 1 - x));
  }
  double cte,
  dM = 1, 
  C = 1, 
  D = 0, 
  C_D = 0, 
  f = 1,
  den = 1;
  cte = exp(a * log(x) + b * log(1 - x) - log(a) - 
    lgamma(a) - lgamma(b) + lgamma(a + b));
  int m, i = 0;
  do{
    m = i / 2;
    if (i == 0){
      dM = 1;
    } else if (i % 2 == 0){
      dM = (m * (b - m) * x) /((a + 2 * m - 1) * (a + 2 * m));
    } else{
      dM = -((a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1));
    }
    if (abs(den) < 1e-30){
      D = 1e30;
    } else {
      den = 1 + dM * D;
      D = 1 / den;
    }
    if (abs(C) < 1e-30){
      C = 1e-30;
    } else {
      C = 1 + dM / C;
    }
    C_D = C * D;
    f *= C_D;
    i++;
  } while(i < 300 && abs(1 - C_D) > 1e-5);
  if (i >= 300){
    return(datum::nan);
  }
  return(cte * (f - 1));
}
// Regression with nice output
void regressTable(REGmodel m, vector<string>& table){
  //    vec beta, stdBeta, eOut;
  //    double BIC, AIC, AICc;
  //    regress(y, X, beta, stdBeta, eOut, BIC, AIC, AICc);
  char str[70];
  table.clear();
  table.push_back("-------------------------------------------------------------\n");
  snprintf(str, 70, " Regression with %d exogenous\n", (int)m.beta.n_rows);
  table.push_back(str);
  table.push_back("----------------------------------------------------\n");
  table.push_back("            Param       S.E.        |T|    P-value  \n");
  table.push_back("----------------------------------------------------\n");
  mat tp;
  if (m.beta.n_elem > 0){
      tp = join_rows(m.beta, m.stdBeta, abs(m.beta / m.stdBeta));
      // Table of numbers
      char namePar[10];
      for (unsigned i = 0; i < tp.n_rows; i++){
        if (all(m.X.row(i)) == 1)
          snprintf(namePar, 10, "Cnst");
        //            namePar = "Cnst";
        else
          snprintf(namePar, 10, "%s(%d)", "Beta", i);
        //            namePar = "Beta";
        if (abs(tp(i, 0)) > 1e-4){
          snprintf(str, 70, "%*s: %d %12.4f %12.4f\n", 10, namePar, i, tp(i, 0), tp(i, 1));
          //            snprintf(str, 70, "%*s(%d): %12.4f %12.4f %12.4f %10.6f\n", 10, namePar.c_str(), i, tp(i, 0), tp(i, 1), tp(i, 2), 0);
        } else {
          snprintf(str, 70, "%*s: %d %12.3e %12.3e\n", 10, namePar, i, tp(i, 0), tp(i, 1));
          //            snprintf(str, 70, "%*s(%d): %12.3e %12.3e %12.4f %10.6f\n", 10, namePar.c_str(), i, tp(i, 0), tp(i, 1), tp(i, 2), 0);
        }
        table.push_back(str);
      }
  }
  table.push_back("-------------------------------------------------------------\n");
  snprintf(str, 70, "  AIC: %12.4f   BIC: %12.4f   AICc: %12.4f\n", m.AIC, m.BIC, m.AICc);
  table.push_back(str);
  uvec ind = find_finite(m.e);
  m.e = m.e(ind);
  m.e = m.e - mean(m.e);
  vec sigma2 = m.e.t() * m.e / m.e.n_elem;
  double logLik = -(m.e.n_elem / 2.0) * (log(2.0 * datum::pi * sigma2(0)) + 1.0);
  snprintf(str, 70, "           Log-Likelihood: %12.4f\n", logLik);
  table.push_back(str);
  table.push_back("-------------------------------------------------------------\n");
  table.push_back("   Summary statistics:\n");
  table.push_back("-------------------------------------------------------------\n");
  if (m.e.n_elem < 5){
    table.push_back("  All innovations are NaN!!\n");
  } else {
    outputTable(m.e, table);
  }
  table.push_back("-------------------------------------------------------------\n");
}
// Estimation table
void outputTable(vec v, vector<string>& table){
  int nNan, n = v.n_elem;
  vec vn = v - nanMean(v);
  vn = removeNans(vn, nNan);
  uvec outliers = find(abs(vn) > 2.7 * stddev(vn));
  int nOutliers = outliers.n_elem;
  vec vClean = v;
  vClean(outliers).fill(datum::nan);
  // Autocorrelations
  int n4 = n / 4;
  uvec ind;
  if (n4 < 168){
    if (n < 30){
      //ind << 0 << 1 << 2 << 3 << endr;
      ind = {0, 1, 2, 3};
    } else {
      //ind << 0 << 3 << 7 << 11 << endr;
      ind = {0, 3, 7, 11};
    }
  } else {
    //ind << 0 << 3 << 7 << 11 << 23 << 167 << endr;
    ind = {0, 3, 7, 11, 23, 167};
  }
  int ncoef = ind.max() + 1;
  vec acfCoef1, acfCoef2;
  acf(v, ncoef, acfCoef1);
  vec qq= n * (n + 2) * cumsum((acfCoef1 % acfCoef1) / (n - regspace(1, ncoef)));
  vec Q1= qq.rows(ind), Q2;
  if (nOutliers > 0){
    acf(vClean, ncoef, acfCoef2);
    qq= n * (n + 2) * cumsum((acfCoef2 % acfCoef2) / (n - regspace(1, ncoef)));
    Q2= qq.rows(ind);
  } else{
    acfCoef2.zeros(ncoef);
  }
  // Gaussianity
  double bj1 = datum::nan, pbj1 = datum::nan, bj2 = datum::nan, pbj2 = datum::nan;
  beraj(v, bj1, pbj1);
  if (nOutliers > 0){
    beraj(vClean, bj2, pbj2);
  }
  // Heteroskedasticity
  double F1, pF1, F2 = datum::nan, pF2 = datum::nan;
  int df1, df2;
  heterosk(v, F1, pF1, df1);
  if (nOutliers > 0){
    heterosk(vClean, F2, pF2, df2);
  }
  // Table
  char str[70];
  ind = ind + 1;
  snprintf(str, 70, "        Observations: %5.0i        Missing data: %5.0i\n", (int)v.n_elem, nNan);
  table.push_back(str);
  snprintf(str, 70, "        Q(%2.0i): %12.4f         Q(%2.0i): %12.4f\n", (int)ind(0), Q1(0), (int)ind(1), Q1(1));
  table.push_back(str);
  snprintf(str, 70, "        Q(%2.0i): %12.4f         Q(%2.0i): %12.4f\n", (int)ind(2), Q1(2), (int)ind(3), Q1(3));
  table.push_back(str);
  if (ind.n_rows > 4){
    snprintf(str, 70, "        Q(%2.0i): %12.4f         Q(%2.0i): %12.4f\n", (int)ind(4), Q1(4), (int)ind(5), Q1(5));
    table.push_back(str);
  }
  snprintf(str, 70, "  Bera-Jarque: %12.4f       P-value: %12.4f\n", bj1, pbj1);
  table.push_back(str);
  snprintf(str, 70, "       H(%3.0i): %12.4f       P-value: %12.4f\n", df1, F1, pF1);
  table.push_back(str);
  if (nOutliers > 0){
    snprintf(str, 70, "  Outliers (>2.7 ES): %5.0i\n", nOutliers);
    table.push_back(str);
    snprintf(str, 70, "        Q(%2.0i): %12.4f         Q(%2.0i): %12.4f\n", (int)ind(0), Q2(0), (int)ind(1), Q2(1));
    table.push_back(str);
    snprintf(str, 70, "        Q(%2.0i): %12.4f         Q(%2.0i): %12.4f\n", (int)ind(2), Q2(2), (int)ind(3), Q2(3));
    table.push_back(str);
    if (ind.n_rows > 4){
      snprintf(str, 70, "        Q(%2.0i): %12.4f         Q(%2.0i): %12.4f\n", (int)ind(4), Q2(4), (int)ind(5), Q2(5));
      table.push_back(str);
    }
    snprintf(str, 70, "  Bera-Jarque: %12.4f       P-value: %12.4f\n", bj2, pbj2);
    table.push_back(str);
    snprintf(str, 70, "       H(%3.0i): %12.4f       P-value: %12.4f\n", df2, F2, pF2);
    table.push_back(str);
  }
}
// Estimation table
void outputTable(mat v, vector<string>& table){
    int nNan, n = v.n_cols, n4, ncoef, df1;
    vec vn, vClean, acfCoef1, qq, Q1, VnNan(v.n_rows),
        vBj(v.n_rows), vPBj(v.n_rows), vF1(v.n_rows), pvF1(v.n_rows);
    uvec ind, outliers;
    uword nOutliers = 0; nOutliers++;
    double bj1, pbj1, F1, pF1; //, bj2, pbj2, F2 = datum::nan, pF2 = datum::nan;
    n4 = n / 4;
    if (n4 < 168){
        if (n < 30){
            //ind << 0 << 1 << 2 << 3 << endr;
            ind = {0, 1, 2, 3};
        } else {
            //ind << 0 << 3 << 7 << 11 << endr;
            ind = {0, 3, 7, 11};
        }
    } else {
        //ind << 0 << 3 << 7 << 11 << 23 << 167 << endr;
        ind = {0, 3, 7, 11, 23, 167};
    }
    ncoef = ind.max() + 1;
    mat vQ1(ind.n_rows, v.n_rows);
    nOutliers = 0;
    for (uword i = 0; i < v.n_rows; i++){
        vn = v.row(i).t();
        vn -= nanMean(vn);
        vn = removeNans(vn, nNan);
        VnNan(i) = nNan;
        outliers = find(abs(vn) > 2.7 * stddev(vn));
        nOutliers = outliers.n_elem;
        vClean = vn;
        vClean(outliers).fill(datum::nan);
        // Autocorrelations
        // n4 = n / 4;
        // ind;
        // if (n4 < 168){
        //     if (n < 30){
        //         //ind << 0 << 1 << 2 << 3 << endr;
        //         ind = {0, 1, 2, 3};
        //     } else {
        //         //ind << 0 << 3 << 7 << 11 << endr;
        //         ind = {0, 3, 7, 11};
        //     }
        // } else {
        //     //ind << 0 << 3 << 7 << 11 << 23 << 167 << endr;
        //     ind = {0, 3, 7, 11, 23, 167};
        // }
        // ncoef = ind.max() + 1;
        acf(vn, ncoef, acfCoef1);
        qq= n * (n + 2) * cumsum((acfCoef1 % acfCoef1) / (n - regspace(1, ncoef)));
        Q1= qq.rows(ind);
        vQ1.col(i) = Q1;
        // if (nOutliers > 0){
        //     acf(vClean, ncoef, acfCoef2);
        //     qq= n * (n + 2) * cumsum((acfCoef2 % acfCoef2) / (n - regspace(1, ncoef)));
        //     Q2= qq.rows(ind);
        // } else{
        //     acfCoef2.zeros(ncoef);
        // }
        // Gaussianity
        bj1 = datum::nan, pbj1 = datum::nan; //, bj2 = datum::nan, pbj2 = datum::nan;
        beraj(vn, bj1, pbj1);
        vBj(i) = bj1;
        vPBj(i) = pbj1;
        // uvec outliers, ind;
        // if (nOutliers > 0){
        //     beraj(vClean, bj2, pbj2);
        // }
        // Heteroskedasticity
        heterosk(vn, F1, pF1, df1);
        vF1(i) = F1;
        pvF1(i) = pF1;
        // if (nOutliers > 0){
        //     heterosk(vClean, F2, pF2, df2);
        // }
    }
    // // Table
    // // char str[70];
    // char str; //[6 * v.n_rows + 25];
    // printVec(VnNan.t(), "        Missing data:", " %5.0i", str);

    // // int offset = 0;

    // // offset += snprintf(str + offset, sizeof(str) - offset, "Vector data:");

    // // for (size_t i = 0; i < VnNan.n_elem; ++i) {
    // //     offset += snprintf(str + offset, sizeof(str) - offset, " %.2f", VnNan(i));
    // // }



    // // string f1 = "%s", f2 = "%s";
    // // for (uword i = 0; i < v.n_rows; i++){
    // //     f1 += " %6.0i";
    // //     f2 += " %3.2f";
    // // }





    // // ind = ind + 1;
    // // const int width = 5 * v.n_rows + 25;
    // // // char buffer[width];
    // // char str[width];
    // // // string str = "        Missing data:";
    // // // str += snprintf(buffer, width, "%d", VnNan(0));

    // // // snprintf(str, allW, "        Missing data:");
    // // // snprintf(str, allW, "%5d", VnNan(0));
    // // string format = "Missing data: " + f1 + "\n";
    // // snprintf(str, width, format, VnNan);

    // // for (uword i = 0; i < v.n_rows; i++){snprintf(str, sizeof(str) + width, "%*d", width, VnNan(i));}
    // table.push_back(str);

    // snprintf(str, 70, "        Observations: %5.0i        Missing data: %5.0i\n", (int)v.n_elem, nNan);
    // table.push_back(str);
    // snprintf(str, 70, "        Q(%2.0i): %12.4f         Q(%2.0i): %12.4f\n", (int)ind(0), Q1(0), (int)ind(1), Q1(1));
    // table.push_back(str);
    // snprintf(str, 70, "        Q(%2.0i): %12.4f         Q(%2.0i): %12.4f\n", (int)ind(2), Q1(2), (int)ind(3), Q1(3));
    // table.push_back(str);
    // if (ind.n_rows > 4){
    //     snprintf(str, 70, "        Q(%2.0i): %12.4f         Q(%2.0i): %12.4f\n", (int)ind(4), Q1(4), (int)ind(5), Q1(5));
    //     table.push_back(str);
    // }
    // snprintf(str, 70, "  Bera-Jarque: %12.4f       P-value: %12.4f\n", bj1, pbj1);
    // table.push_back(str);
    // snprintf(str, 70, "       H(%3.0i): %12.4f       P-value: %12.4f\n", df1, F1, pF1);
    // table.push_back(str);
    // if (nOutliers > 0){
    //     snprintf(str, 70, "  Outliers (>2.7 ES): %5.0i\n", nOutliers);
    //     table.push_back(str);
    //     snprintf(str, 70, "        Q(%2.0i): %12.4f         Q(%2.0i): %12.4f\n", (int)ind(0), Q2(0), (int)ind(1), Q2(1));
    //     table.push_back(str);
    //     snprintf(str, 70, "        Q(%2.0i): %12.4f         Q(%2.0i): %12.4f\n", (int)ind(2), Q2(2), (int)ind(3), Q2(3));
    //     table.push_back(str);
    //     if (ind.n_rows > 4){
    //         snprintf(str, 70, "        Q(%2.0i): %12.4f         Q(%2.0i): %12.4f\n", (int)ind(4), Q2(4), (int)ind(5), Q2(5));
    //         table.push_back(str);
    //     }
    //     snprintf(str, 70, "  Bera-Jarque: %12.4f       P-value: %12.4f\n", bj2, pbj2);
    //     table.push_back(str);
    //     snprintf(str, 70, "       H(%3.0i): %12.4f       P-value: %12.4f\n", df2, F2, pF2);
    //     table.push_back(str);
    // }
}
// Print vec with format
// void printVec(rowvec data, const string intro, const string format, char& str){
//     int offset = 0;
//     offset += snprintf(str + offset, sizeof(str) - offset, intro);
//     for (size_t i = 0; i < data.n_elem; ++i) {
//         offset += snprintf(str + offset, sizeof(str) - offset, format, data(i));
//     }
//     str += "\n";
// }

