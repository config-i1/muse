/**************************************
 State Space systems class
 Needs Armadillo
***************************************/
// #include <armadillo>
// using namespace arma;
/***************************************************
  * Data structures
****************************************************/
// Model matrices
struct SSmatrix{
    // system matrices of a general State Space model
    // y(t)   = Z a(t) + D   u(t) + C eps(t)
    // a(t+1) = T a(t) + Gam u(t) + R eta(t)
    // Var(eps(t)) = H;  Var(eta(t)) = Q; Cov(eta(t), eps(t)) = S
    mat T, R, Q, Z, C, H, Gam, D, S;   // Gam and S not implemented
    uword ns,               // number of states
          ny;               // number of outputs
    int   n;                // length of time series
    uvec tvMatrices;        // tvp matrices {T, Gam, R, Q, Z, D, C, H, S}
    bool tvp = false;       // system time varying
    // mat delta;              // distance between observations for CT models
    bool identityR = false,  // RQR = Q
         identityC = false; // CHC = H
         // constantRQR = false, // constant RQR
         // constantCHC = false;  // constant CHC
};
// Model structure
struct SSinputs{
   // Inputs
   vec y,                 // output data
       p,                 // vector of parameter values
       pTransform,        // un-concentrated transformed parameters
       p0,                // vector of initial values for parameters
       stdP;              // standard errors of parameters
   mat u;                 // input data
   int h = 24;            // forecast horizon
   bool cLlik = true;    // concentrated log-likelihood on / off
   // user function implementing the model
   std::function <void (vec, SSmatrix*, void*)> userModel;
   void* userInputs;      // inputs needed by the user model
   // Outputs
   vec v,                 // innovations
       yFit,              // fitted values
       F,                 // Variance of fitted values
       yFor,              // output forecasts
       FFor,              // Variance of forecasts
       betaAug,           // betas of augmented KF (including initial states)
       betaAugVar,        // variances of betaAug (idag(iSn))
       criteria,          // identification criterion
       coef;              // Coefs for coef function
   mat a,                 // estimated states
       P,                 // variances of states
       eta,               // estimates of transition perturbations
       covp;              // Covariance of parameters
   SSmatrix system;       // system matrices
   double objFunValue,    // value of objective function at optimum
          outlier;        // critical value for outlier detection
   string estimOk;        // type of estimation convergence
   vector<string> table;  // output table from evaluate()
   // Needed for other purposes
   vec Finf,              // innovation variance before colapsing
       aEnd,              // final state vector estimated
       iF,                // inverse of F
       grad,              // gradient at optimum
       rNrOut;            // Needed for outlier detection
   mat K,                 // Kalman gain for smoothing
       Kinf,              // Kalman gain before colapsing
       PEnd,              // final P estimated
       rOut;              // Needed for outlier detection
   cube NOut;             // Needed for outlier detection
   int d_t = 0,           // colapsing observation
       nonStationaryTerms, // number of non stationary terms in state vector
       flag,              // output of optimization algorithm
       Iter;              // Number of iterations in estimation
   double innVariance;    // innovations variance
   bool exact = true,     // exact or numerical gradient
        verbose,          // intermediate output verbose on / off
        augmented = false; // Augmented KF estimation on / off
   std::function <double (vec&, void*)> llikFUN; // LogLik to select llik or llikAug
};
/****************************************************
 * Defining SSmodel class
 ****************************************************/
// SS system class
class SSmodel{
  protected:
    SSinputs inputs;
  public:
    // Constructor declarations
    SSmodel(){}
    SSmodel(SSinputs);
    SSmodel(SSinputs, SSmatrix);
    // Destructor
    ~SSmodel();
    // Evaluate llik
    double evalLlik(vec);
    // Estimate by Maximum-Likelihood
    void estim();
    void estim(vec);
    // Forecasting system
    void forecast();
    // Kalman filter pass
    void filter();
    void filter(unsigned int);
    // Smoothing pass
    void smooth(bool);
    // Disturbance pass
    void disturb();
    // Validation
    void validate(bool, double);
    // // Getters and setters
    // Get inputs
    SSinputs getInputs(){
      return inputs;
    }
    // Set inputs
    void setInputs(SSinputs inputs){
      this->inputs = inputs;
    }
    // Set system matrices
    void setSystemMatrices(){
      inputs.userModel(inputs.p, &inputs.system, inputs.userInputs);
    }
    // Get Objective function (after estimation)
    double getObjFunValue(){
      return inputs.objFunValue;
    }
    // Print inputs on screen
     void print();
};
/***************************************************
 * Auxiliar function declarations
 ****************************************************/
// Check stationarity of transition matrix (KFinit)
void isStationary(mat&, uvec&);
// Initialize Kalman Filter (llik)
void KFinit(mat&, mat&, uword, vec&, mat&, mat&);
// Matrix operations in KF (llik)
void MFK(mat&, vec&, mat, vec&, mat&, vec&);
void aP(vec& at, mat& Pt, vec& Kt, vec& vt, vec& Mt);
// Correction step of KF (llik)
void KFcorrection(bool, bool, bool, bool, SSinputs*, mat, mat&, vec&,
                  mat&, vec&, double, mat&, mat&, mat&, vec&, uword,
                  vec&, mat&, mat);
// Auxiliar function for computing llik in KF (llik)
void llikCompute(bool, mat, mat, mat, mat, mat&, mat&, mat&);
// Prediction stage in KF (llik)
void KFprediction(bool steadyState, bool colapsed, mat& T, mat& RQRt, vec& at, mat& Pt, mat& Pinft);
// Compute log-likelihood
double llik(vec&, void*);
// Compute log-likelihood with eXogenous inputs
double llikAug(vec&, void*);
// Select differentials (increments)
vec differential(vec p);
// Analytic and numeric gradient of log-likelihood
vec gradLlik(vec&, void*, double, int&);
// Llik hessian (for parameter covariances)
mat hessLlik(void*);
// True filter/smooth/disturb function
void auxFilter(unsigned int, SSinputs&);
// solution to lyapunto equation P = Phi * P * Phi' + Q
mat dlyap(mat Phi, mat Q);
// Check consistency of SS system matrices in SSmatrix
bool checkSystem(SSmatrix&);
// Set time varying matrices to time t
void setT(uword, SSmatrix&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&, mat&);
/****************************************************
 // SS implementations for univariate SS systems
 ****************************************************/
// Constructors with inputs
SSmodel::SSmodel(SSinputs inputs){
  this->inputs = inputs;
}
SSmodel::SSmodel(SSinputs inputs, SSmatrix system){
  this->inputs = inputs;
  this->inputs.system = system;
}
// Destructor
SSmodel::~SSmodel(){}
// Print inputs on screen
 // void SSmodel::print(){
 //   cout << "**************************" << endl;
 //   cout << "Start of SS system:" << endl;
 //   inputs.y.t().print("y:");
 //   inputs.p.t().print("p:");
 //   inputs.p0.t().print("p0:");
 //   inputs.stdP.t().print("stdP:");
 //   inputs.u.print("u:");
 //   cout << "h: " << inputs.h << endl;
 //   cout << "cLlik: " << inputs.cLlik << endl;
 //   inputs.v.t().print("v:");
 //   inputs.betaAug.t().print("betaAug:");
 //   inputs.criteria.print("criteria:");
 //   inputs.eta.t().print("eta:");
 //   cout << "objFunValue: " << inputs.objFunValue << endl;
 //   cout << "outlier: " << inputs.outlier << endl;
 //   cout << "estimOk: " << inputs.estimOk << endl;
 //   inputs.grad.t().print("grad:");
 //   cout << "d_t: " << inputs.d_t << endl;
 //   cout << "nonStationaryTerms: " << inputs.nonStationaryTerms << endl;
 //   cout << "flag: " << inputs.flag << endl;
 //   cout << "innVariance: " << inputs.innVariance << endl;
 //   cout << "exact: " << inputs.exact << endl;
 //   cout << "verbose: " << inputs.verbose << endl;
 //   cout << "augmented: " << inputs.augmented << endl;
 //   cout << "End of SS system:" << endl;
 //   cout << "**************************" << endl;
 // }
// Evaluateing llik
double SSmodel::evalLlik(vec p){
    return inputs.llikFUN(p, &inputs);
}
// Estimation by Maximum-Likelihood
void SSmodel::estim(){
  SSmodel::estim(inputs.p0);
}
void SSmodel::estim(vec p){
  double objFunValue;
  vec grad;
  mat iHess;
  this->inputs.p0 = p;
  wall_clock timer;
  timer.tic();
  int flag = quasiNewton(inputs.llikFUN, gradLlik, p, &inputs, objFunValue, grad, iHess, inputs.verbose);
  // Information criteria
  uvec indNan = find_nonfinite(inputs.y);
  int nNan = inputs.y.n_elem - indNan.n_elem;
  double LLIK, AIC, BIC, AICc;
  LLIK = -0.5 * nNan * (log(2*datum::pi) + objFunValue);
  infoCriteria(LLIK, p.n_elem + inputs.nonStationaryTerms, nNan,
               AIC, BIC, AICc);
  vec criteria(4);
  criteria(0) = LLIK;
  criteria(1) = AIC;
  criteria(2) = BIC;
  criteria(3) = AICc;
  this->inputs.criteria = criteria;
  if (!isfinite(objFunValue))
      flag = 0;
  // Printing results
  if (flag == 1) {
    this->inputs.estimOk = "Q-Newton: Gradient convergence.\n";
  } else if (flag == 2){
    this->inputs.estimOk = "Q-Newton: Function convergence.\n";
  } else if (flag == 3){
      this->inputs.estimOk = "Q-Newton: Parameter convergence.\n";
  } else if (flag == 4){
      this->inputs.estimOk = "Q-Newton: Maximum number of iterations reached.\n";
  } else if (flag == 5){
      this->inputs.estimOk = "Q-Newton: Maximum number of Function evaluations.\n";
  } else if (flag == 6){
      this->inputs.estimOk = "Q-Newton: Unable to decrease objective function.\n";
  } else if (flag == 7){
      this->inputs.estimOk = "Q-Newton: Objective function returns nan.\n";
  } else {
      this->inputs.estimOk = "Q-Newton: No convergence!!\n";
  }
  if (inputs.verbose){
    double nSeconds = timer.toc();
    printf("%s", this->inputs.estimOk.c_str());
    printf("Elapsed time: %10.5f seconds\n", nSeconds);
  }
  this->inputs.p = p;
  this->inputs.objFunValue = objFunValue;
  this->inputs.grad = grad;
  this->inputs.flag = flag;
  this->inputs.v.reset();
}
// Forecasting system
void SSmodel::forecast(){
  // mat RQRt = inputs.system.R * inputs.system.Q * inputs.system.R.t(),
  //     CHCt = inputs.system.C * inputs.system.H * inputs.system.C.t();
  int n = SSmodel::inputs.y.n_elem, k = SSmodel::inputs.u.n_rows;
  inputs.yFor.zeros(inputs.h);
  inputs.FFor.zeros(inputs.h);
  vec at = inputs.aEnd;
  mat Pt;
  if (at.has_nan()){
    filter();
    inputs.yFor = SSmodel::inputs.yFit.tail_rows(SSmodel::inputs.h);
    inputs.FFor = SSmodel::inputs.F.tail_rows(SSmodel::inputs.h);
  } else {
    mat T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt;
    // if (abs(inputs.innVariance - 1) > 1e-4){
      Pt = inputs.PEnd * inputs.innVariance;
    // } else {
    //   Pt = inputs.PEnd; // * inputs.innVariance;
    // }
    mat P0 = Pt;
    // mat Z = inputs.system.Z.row(0);
    uword t = inputs.y.n_elem;
    bool TVP = (inputs.system.Z.n_rows > 1);
    if (k > 0 && inputs.system.Z.n_rows == 1){
        int npar = inputs.betaAug.n_elem;
        inputs.system.D = inputs.betaAug.rows(npar - k, npar - 1).t();
    }
    ////////////////
    /// \brief setT
    ///// Initiialising constant system matrices for t = 0
    setT(0, inputs.system, T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt);
    /////////////////
    for (uword i = 0; i < (uword)inputs.h; i++){
      ////////////////
      setT(t + i, inputs.system, T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt);
      /////////////////
      // if (TVP)
      //     Z = inputs.system.Z.row(t + i);
      inputs.yFor(span(i)) = Z * at;
      if (k == 0){
          inputs.yFor(span(i)) += inputs.system.D;
      } else if (!TVP) {
          inputs.yFor(span(i)) += inputs.system.D * SSmodel::inputs.u.col(n + i);
      }
      inputs.FFor(span(i)) = Z * Pt * Z.t() + CHCt;
      KFprediction(false, true, T, RQRt, at, Pt, P0);
    }
    // if (abs(inputs.innVariance - 1) < 1e-4){
    //   inputs.FFor *= inputs.innVariance;
    // }
  }
}
// Kalman filter pass
void SSmodel::filter(){
  SSmodel::filter(0);
}
void SSmodel::filter(unsigned int smooth){
    auxFilter(smooth, inputs);
}
// Smoothing pass
void SSmodel::smooth(bool outlier){
  if (outlier){
    SSmodel::filter(3);
  } else {
    SSmodel::filter(1);
  }
}
// Disturbance pass
void SSmodel::disturb(){
  SSmodel::filter(2);
}
// Validation
void SSmodel::validate(bool estimateHess, double nPar){
  // Input is inverse of Hessian. Calculated if empty
  uvec auxx;
  // Inverse of hessian and covariance of parameters
  int k = inputs.p.n_elem;
  uvec nn = find_finite(inputs.y);
  mat hess = eye(k, k);
  mat iHess = hess, table0;
  vec t, pValue(k);
//  if (k > 0){
      if (estimateHess){
          hess = hessLlik(&inputs) * 0.5 * nn.n_elem;
          iHess.fill(datum::nan);
          if (hess.is_finite()){
              iHess = pinv(hess);
              iHess.diag() = abs(iHess.diag());
          }
      }
      inputs.stdP = sqrt(iHess.diag());
      t = abs(inputs.p / inputs.stdP);
      pValue = 2 * (1- tCdf(t, nn.n_elem - k));
      uvec aux = find(t > 1000);
      if (aux.n_elem > 0){
        t(aux).fill(datum::inf);
        pValue(aux).fill(0);
      }
      table0 = join_horiz(join_horiz(join_horiz(inputs.p, inputs.stdP), t), pValue);
//  }
  // First part of table
  char str[70];
  inputs.table.clear();
  inputs.table.push_back("-------------------------------------------------------------\n");
  snprintf(str, 70, " %s", inputs.estimOk.c_str());
  inputs.table.push_back(str);
  inputs.table.push_back("-------------------------------------------------------------\n");
  inputs.table.push_back("            Param       S.E.        |T|    P-value     |Grad| \n");
  inputs.table.push_back("-------------------------------------------------------------\n");
  for (unsigned int i = 0; i < nPar; i++){
    snprintf(str, 70, "       %10.4f %10.4f %10.4f %10.4f %10.6f\n", table0(i, 0), table0(i, 1), table0(i, 2), table0(i, 3), abs(inputs.grad(i)));
    inputs.table.push_back(str);
  }
  // Adding inputs betas
  int nu = inputs.u.n_rows;
  if (nu > 0){
    vec betas, stdBetas;
    if (inputs.system.Z.n_rows > 1){
        filter();
        betas = inputs.aEnd.rows(inputs.aEnd.n_rows - nu, inputs.aEnd.n_rows - 1);
        vec dPEnd = inputs.PEnd.diag();
        stdBetas = sqrt(dPEnd.rows(inputs.aEnd.n_rows - nu, inputs.aEnd.n_rows - 1));
    } else {
        int ind = inputs.betaAug.n_elem - nu;
        betas = inputs.betaAug.rows(ind, ind + nu - 1);
        stdBetas = sqrt(inputs.betaAugVar.rows(ind, ind + nu - 1));
    }
    vec tBetas = betas / stdBetas;
    vec pValueBetas = 2 * (1- tCdf(tBetas, nn.n_elem - k));
    for (int i = 0; i < nu; i++){
      snprintf(str, 70, "       %10.4f %10.4f %10.4f %10.4f %10.6f\n", betas(i),
              stdBetas(i), tBetas(i), pValueBetas(i), datum::nan);
      inputs.table.push_back(str);
    }
  }
  uvec ind = find_finite(inputs.y);
  inputs.table.push_back("-------------------------------------------------------------\n");
  snprintf(str, 70, "  AIC: %12.4f   BIC: %12.4f   AICc: %12.4f\n", inputs.criteria(1), inputs.criteria(2), inputs.criteria(3));
  inputs.table.push_back(str);
  snprintf(str, 70, "           Log-Likelihood: %12.4f\n", inputs.criteria(0));
  inputs.table.push_back(str);
  inputs.table.push_back("-------------------------------------------------------------\n");
  // Recovering innovations for tests
  if (inputs.augmented)
    llikAug(inputs.p, &inputs);
  else
    llik(inputs.p, &inputs);
  filter();
  //Second part of table
  inputs.table.push_back("   Summary statistics:\n");
  inputs.table.push_back("-------------------------------------------------------------\n");
  auxx = find_finite(inputs.v);
  if (auxx.n_elem < 5){
    inputs.table.push_back("  All innovations are NaN!!\n");
  } else {
    outputTable(inputs.v, inputs.table);
  }
  inputs.table.push_back("-------------------------------------------------------------\n");
  // Show Table
  // if (show){
      // // for (auto i = inputs.table.begin(); i != inputs.table.end(); i++){
      // //   cout << *i << " ";
      // // }
      // for (unsigned int i = 0; i < inputs.table.size(); i++){
      //   printf("%s ", inputs.table[i].c_str());
      // }
  // }
}
/*************************************************************
 * Implementation of auxiliar functions
 ************************************************************/
// Initializing Kalman Filter
void KFinit(mat& T, mat& RQRt, uword ns, vec& at, mat& Pt, mat& Pinft){
  at.zeros(ns);
  Pt.zeros(ns, ns);
  vec Pinfdiag; Pinfdiag.ones(ns);
  uvec stat;
  isStationary(T, stat);
  if (!stat.is_empty()){
    // int Ns = stat.n_elem;
    Pinfdiag.elem(stat).zeros();
    // Lyapunov for stationary elements
    mat q = RQRt(stat, stat);
    mat t = T(stat, stat);
    mat P2 = dlyap(t, q);
    // int Ns = stat.n_elem, Ns2 = Ns * Ns; mat P2bis = reshape(pinv(eye(Ns2, Ns2) - kron(t, t)) * vectorise(q), Ns, Ns);
    Pt(stat, stat) = P2;
  }
  Pinft = diagmat(Pinfdiag);
}
// Check stationarity of transition matrix (Kfinit)
void isStationary(mat& T, uvec& stat){
  int n = T.n_rows;
  cx_vec eigval(n);
  vec nons, nonstat;
  cx_mat V(n, n);
  double tol = 0.98;
  nons.zeros(n);
  nonstat = nons;
  eig_gen(eigval, V, T);
  nons.elem(find(abs(eigval) >= tol)).ones();
  nonstat.elem(find(abs(V) * nons > 0)).ones();
  stat = find(1 - nonstat);
}
// Update of Mt, Ft and Kt
void MFK(mat& Pt, mat& Z, mat CHCt, vec& Mt, mat& Ft, vec& Kt){
  Mt = Pt * Z.t();
  Ft = Z * Mt + CHCt;
  Kt = Mt / Ft(0, 0);
}
// Update of at and Pt
void aP(vec& at, mat& Pt, vec& Kt, vec& vt, vec& Mt){
  at = at + Kt * vt;
  Pt = Pt - Kt * Mt.t();
}
// Correction step in Kalman Filtering for every t
void KFcorrection(bool miss, bool colapsed, bool steadyState, bool smooth,
                  SSinputs* data, mat CHCt,
                  mat& Finft, vec& vt, double Dt, mat& Ft, mat& iFt, vec& at, mat& Pt,
                  mat& Pinft, vec& Kt, uword t, vec& auxFinf, mat& auxKinf, mat Z){
  vec Mt, Minft, Kinft;
  mat KK;
  rowvec yt = data->y.row(t);
  mat iFinft;
  if (miss){
        vt.fill(datum::nan);
        Kt.fill(0);
        Ft = Z * Pt * Z.t() + CHCt;
    } else {
      vt = yt - Z * at - Dt;
      if (steadyState){
        at = at + Kt * vt;
      } else if(colapsed) {
        MFK(Pt, Z, CHCt, Mt, Ft, Kt);
        aP(at, Pt, Kt, vt, Mt);
        iFt = 1 / Ft(0, 0);
      } else {
        Minft = Pinft * Z.t();
        Finft = Z * Minft;
        if (data->exact || smooth) auxFinf.row(t) = Finft;
        if (Finft(0, 0) > 1e-8){
          Mt = Pt * Z.t();
          Ft = Z * Mt + CHCt;
          iFinft = 1 / Finft(0, 0);
          iFt = iFinft;
          Kinft = Minft * iFinft;
          if (data->exact || smooth) auxKinf.col(t) = Kinft;
          Kt = (Mt - Kinft * Ft) * iFinft;
          aP(at, Pinft, Kinft, vt, Minft);
          KK = Mt * Kinft.t();
          Pt = Pt + Kinft * Ft * Kinft.t() - (KK + KK.t());
        } else {
          MFK(Pt, Z, CHCt, Mt, Ft, Kt);
          aP(at, Pt, Kt, vt, Mt);
          iFt = 1 / Ft(0, 0);
        }
      }
    }
}
// Llik computation inside llik
void llikCompute(bool colapsed, mat Finft, mat vt, mat Ft, mat iFt,
                 mat& v2F, mat& logF, mat& llikValue){
  if (colapsed || Finft(0, 0) < 1e-8){
    v2F  += vt * iFt * vt;
    logF += log(Ft);
  } else {
    llikValue += log(Finft);
  }
}
// KF Prediction step in Kalman filtering for every t
void KFprediction(bool steadyState, bool colapsed, mat& T, mat& RQRt, vec& at, mat& Pt, mat& Pinft){
  at = T * at;
  if (!steadyState){
    Pt = T * Pt * T.t() + RQRt;
  }
  if (!colapsed){
    Pinft = T * Pinft * T.t();
  }
}
// Compute log-likelihood
double llik(vec& p, void* opt_data){
  // Converting void* to SSinputs*
  SSinputs* data = (SSinputs*)opt_data;
  // Running user function model
  data->userModel(p, &data->system, data->userInputs);
  double tolsta = 1e-19;
  uword n,
        ns = data->system.T.n_cols,
        nMiss = 0;
  ///////////
  mat T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt;
  setT(0, data->system, T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt);
  ///////////
  mat Pt,
      Pinft,
      Ft(1, 1),
      Finft(1, 1),
      iFt(1, 1),
      oldPt,
      llikValue(1, 1),
      logF(1, 1),
      v2F(1, 1),
      auxKinf;
  vec at,
      Kt(ns),
      vt(1),
      auxFinf;
  bool colapsed = false,
       steadyState = false,
       miss = false;
  data->innVariance = 1;
  // Initializing variables
  llikValue.fill(0);
  logF.fill(0);
  v2F.fill(0);
  n = data->y.n_rows;
  data->d_t = n;
  // Kfinit
  // KFinit(data->system.T, RQRt, ns, at, Pt, Pinft);
  KFinit(T, RQRt, ns, at, Pt, Pinft);
  oldPt.zeros(ns, ns);
  data->nonStationaryTerms = sum(Pinft.diag());
  data->v = zeros(n);
  data->F = data->v;
  data->iF = data->v;
  if (data->exact){
    data->K = zeros(ns, n);
    auxFinf = data->v;
    auxKinf = data->K;
  }
  // mat Z = data->system.Z.row(0);
  // bool TVP = false;
  // if (data->system.Z.n_rows > 1)
  //     TVP = true;
  // KF loop
  for (uword t = 0; t < n; t++){
    // Data missing
    if (!is_finite(data->y.row(t))){
      steadyState = false;
      miss = true;
      nMiss += 1;
    } else {
      miss = false;
    }
    ////////////////
    if (t > 0 && data->system.tvp){
        setT(t, data->system, T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt);
    }
    /////////////////
    // if (TVP)
    //     Z = Z = data->system.Z.row(t);
    // Correction
    KFcorrection(miss, colapsed, steadyState, data->exact, data, CHCt,
                 Finft, vt, D(0, 0), Ft, iFt, at, Pt, Pinft,
                 Kt, t, auxFinf, auxKinf, Z);
    // llik calculation
    if (!miss && t < n)
      llikCompute(colapsed, Finft, vt, Ft, iFt, v2F, logF, llikValue);
    // Prediction
    KFprediction(steadyState, colapsed, T, RQRt, at, Pt, Pinft);
    // Storing final state and covariance for forecasting
   if (t == n - 1){
     data->PEnd = Pt;
     data->aEnd = at;
   }
    // Checking colapsed
    if (!colapsed){
      if (all(all(abs(Pinft) < 1e-6))){
        colapsed = true;
        data->d_t = t;
        if (data->exact){
          data->Finf = auxFinf.rows(0, t);
          data->Kinf = auxKinf.cols(0, t);
        }
      }
    }
    // Checking steady state
    if (!steadyState){
      if (colapsed && all(all(abs(Pt - oldPt) < tolsta))){
        steadyState = true;
      } else {
        oldPt = Pt;
      }
    }
    // Storing for analytical derivatives
    data->v.row(t)= vt;
    data->F.row(t)= Ft;
    data->iF.row(t)= iFt;
    if (data->exact){
      data->K.col(t)= Kt;
    }
  }
  // System did not colapsed
  if (data->exact && !colapsed){
    data->Finf = auxFinf;
    data->Kinf = auxKinf;
  }
  // Computing llik value
  int nTrue;
  if (data->d_t < (int)(T.n_rows + 10)){
    // Colapsed KF
    nTrue = n - nMiss - 1 - data->d_t;
  } else {
    // KF did not colapsed
    nTrue = n - nMiss - 1 - T.n_rows;
  }
  if (data->cLlik){         // Concentrated Likelihood
      data->innVariance = v2F(0, 0) / nTrue;
      llikValue = log(data->innVariance) + 1 + (llikValue + logF) / nTrue;
  } else {                  // Crude Likelihood
      llikValue = (llikValue + v2F + logF) / nTrue;
  }
  data->objFunValue = llikValue(0, 0);
  // System colapsed
  if ((uword)data->d_t < n){
      data->v(span(0, data->d_t)).fill(datum::nan);
  }
  return llikValue(0, 0);
}
// Compute log-likelihood for model with eXogenous inputs
double llikAug(vec& p, void* opt_data){
  // Augmented Kalman Filter
  SSinputs* data = (SSinputs*)opt_data;
  // Running user function model (setting system matrices)
  data->userModel(p, &data->system, data->userInputs);
  uword ns = data->system.T.n_rows,
        nMiss = 0,
        n = data->y.n_rows,
        nu = data->u.n_rows,
        k = nu + ns;
  double tolsta = 1e-19;
  ///////////
  mat T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt;
  setT(0, data->system, T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt);
  ///////////
  mat Pt(ns, ns),
      oldPt(ns, ns),
      Ft(1, 1),
      FEnd(1, 1),
      llikValue(1, 1),
      At(ns, k),
      Sn(k, k),
      iSn(k, k),
      AtiSn(ns, k),
      VtiSn(1, k),
      PEndZ(ns ,1); //, W(ns, nu);
  vec at(ns),
      vt(1),
      vEnd(1),
      sn(k),
      beta(k),
      Kt(ns),
      iFt(1),
      viFt(1),
      logF(1),
      v2F(1),
      snBeta(1);
  rowvec Vt(k),
         Xt(k);
  bool miss = false,
       steadyState = false;
  at.fill(0);
  Pt = RQRt;
  oldPt.fill(-10);
  At.fill(0);  // At = -Wt;
  At.submat(0, 0, ns - 1, ns - 1) = -eye(ns, ns);
  sn.fill(0);
  Sn.fill(0);
  iFt.fill(1);
  viFt.fill(0);
  logF.fill(0);
  v2F.fill(0);
  Xt.fill(0);
  uvec aux = regspace<uvec>(ns, k - 1);
  // KF loop
  for (uword t = 0; t < n; t++){
    ////////////////
    if (t > 0 && data->system.tvp){
        setT(t, data->system, T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt);
    }
    /////////////////
    if (nu > 0){
      Xt(aux) = data->u.col(t).t();
    }
    miss = !is_finite(data->y.row(t));
    // Checking for missing
    if (miss){
      nMiss++;
      vt.fill(0);
      Kt.fill(0);
      steadyState = false;
    }
    // Main calculations
    if (steadyState){
      if (!miss){
        vt = data->y.row(t) - Z * at;
      }
    } else {
      Ft = Z * Pt * Z.t() + CHCt;
      iFt = 1 / Ft;
      if (!miss){
        vt = data->y.row(t) - Z * at;
        Kt = T * Pt * Z.t() * iFt;
      }
      Pt = T * Pt * T.t() + RQRt - Kt * Ft * Kt.t();
    }
    at = T * at + Kt * vt;
    // Augmented part
    Vt = Xt - Z * At;
    At = T * At + Kt * Vt; // + Wt;
    if (!miss){
      viFt = vt * iFt;
      sn += Vt.t() * viFt;
      Sn += Vt.t() * iFt * Vt;
      v2F += vt * viFt;
      logF += log(Ft);
    }
    // Checking steady state
    if (!steadyState && t > ns){
      if (all(all(abs(Pt - oldPt) < tolsta))){
        steadyState = true;
      } else {
        oldPt = Pt;
      }
    }
  }
  if (Sn.has_nan() || Sn.has_inf()){
    // Algorithm blew up
    llikValue(0, 0) = datum::nan;
  } else {
    iSn = pinv(Sn);
    // Storing final state and covariance for forecasting
    AtiSn = At * iSn;
    data->aEnd = at - AtiSn * sn;
    data->PEnd = Pt + AtiSn * At.t();
    beta = iSn * sn;
    // Computing llik value
    int nTrue = n - nMiss - k;
    snBeta = sn.t() * beta;
    data->innVariance = (v2F(0, 0) - snBeta(0)) / nTrue;
    llikValue = log(data->innVariance) + 1 + (log(det(Sn)) + logF) / nTrue;
    data->objFunValue = llikValue(0, 0);
    data->betaAug = beta;
    data->betaAugVar = data->innVariance * iSn.diag();
  }
  return llikValue(0, 0);
}
// Select differentials (increments)
vec differential(vec p){
  vec signP = sign(p);
  signP(find(signP == 0)).fill(1);
  return max(join_horiz(abs(p), ones(p.n_elem, 1)), 1) % signP * 1e-8;
}
// Analytic and numeric gradient of log-likelihood
vec gradLlik(vec& p, void* opt_data, double llikValue, int& nFuns){
  int nPar = p.n_elem;
  vec grad(nPar),
      p0 = p,
      inc;
  SSinputs* data = (SSinputs*)opt_data;
  nFuns = 0;
  inc = differential(p);
  if (p.has_nan()){
    grad.fill(datum::nan);
    return grad;
  }
  if (data->exact){  // Analytical derivative
    int ns = data->system.T.n_rows,
        n = data->y.n_elem,
        cQ,
        nMiss = 0;
    mat GammaQ(ns, ns),
        Nt(ns, ns),
        RR(ns, ns),
        sysmatQ,
        sysmatR,
        Z = data->system.Z,
        Gamma(ns + 1, ns + 1),
        Qt,
        dQt,
        dRQRt(ns, ns),
        Inew(ns, ns),
        Lt(ns, ns);
    vec rt(ns),
        vt(1),
        Kt(ns),
        GammaD(1),
        iFt(1),
        e(1),
        D(1),
        Kinft(ns),
        Z_Ft(ns);
    double Finft = 0.0;
    bool colapsed = true;
    // Initialising variables
    GammaQ.fill(0);
    GammaD.fill(0);
    Gamma.fill(0);
    cQ = data->system.Q.n_cols;
    sysmatQ = zeros(cQ + 1, cQ + 1);
    sysmatR = zeros(ns + 1, cQ + 1);
    Qt = zeros(cQ + 1, cQ + 1);
    dQt = sysmatQ;
    e.fill(0);
    D.fill(0);
    Nt = GammaQ;
    rt.fill(0);
    Inew.eye();
    // Main Loop
    for (int t = n - 1; t >= 0; t--){
      if (t <= data->d_t){
        colapsed = false;
      }
      vt = data->v.row(t);
      Kt = data->K.col(t);
      iFt = data->iF.row(t) / data->innVariance;
      if (!colapsed){
        Finft = data->Finf(t); // * data->innVariance;
        Kinft = data->Kinf.col(t);
      }
      if (!is_finite(data->y.row(t))){
        e.fill(0);
        D.fill(0);
        nMiss += 1;
      } else if (colapsed || Finft< 1e-8) {
        e = vt * iFt - Kt.t() * rt;
        D = 1 * iFt + Kt.t() * Nt * Kt;
        Lt = Inew - Kt * Z;
        Z_Ft = Z.t() * iFt(0, 0);
        rt = Z_Ft * vt + Lt.t() * rt;
        Nt = Z_Ft * Z + Lt.t() * Nt * Lt;
      } else {
        e = -Kinft.t() * rt;
        D = Kinft.t() * Nt * Kinft;
        if (Finft >= 1e-8){   // Finf not singular
          Lt = Inew - Kinft * Z;
          rt = Lt.t() * rt;
          Nt = Lt.t() * Nt * Lt;
        }
      }
      GammaD += e * e - D;
      RR = rt * rt.t() - Nt;
      GammaQ += RR;
      rt = data->system.T.t() * rt;
      Nt = data->system.T.t() * Nt * data->system.T;
    }
    // Derivatives of RQRt and CHCt
    sysmatQ(span(0, cQ - 1), span(0, cQ - 1)) = data->system.Q;
    sysmatQ(span(cQ), span(cQ)) = data->system.H;
    sysmatR(span(0, ns - 1), span(0, cQ - 1)) = data->system.R;
    sysmatR(span(ns), span(cQ)) = data->system.C;
    Gamma(span(0, ns - 1), span(0, ns - 1)) = GammaQ;
    Gamma(span(ns), span(ns)) = GammaD;
    int nn = n - nMiss - data->d_t - 1;
    for (int i = 0; i < nPar; i++){
      p0 = p;
      p0.row(i) += inc(i);
      data->userModel(p0, &data->system, data->userInputs);
      Qt(span(0, cQ - 1), span(0, cQ - 1)) = data->system.Q;
      Qt(span(cQ), span(cQ)) = data->system.H;
      dQt= (Qt - sysmatQ) / inc(i);
      dRQRt = sysmatR * dQt * sysmatR.t();
      grad.row(i) = -trace(Gamma * dRQRt) / nn;
    }
    nFuns += 1;
  } else {          // Numerical derivative
    vec F1 = p;
    for (int i = 0; i < nPar; i++){
        p0 = p;
        p0(i) += inc(i);
        F1(i) = data->llikFUN(p0, opt_data);
    }
    grad = (F1 - llikValue) / inc;
    nFuns += nPar;
  }
  return grad;
}
// Llik hessian (for parameter covariances)
mat hessLlik(void* optData){
  SSinputs* inputs = (SSinputs*)optData;
  uword nPar = inputs->p.n_elem;
  vec grad(nPar), p0 = inputs->p, inc(nPar);
  mat Hess(nPar, nPar);
  vec grad0 = inputs->grad;
  inc.fill(1e-5);
  double llikValue2 = 0, llikValue0;  // = inputs->objFunValue;
  if (inputs->augmented){
    llikValue0 = llikAug(p0, inputs);
  } else {
    llikValue0 = llik(p0, inputs);
  }
  Hess.fill(0);
  for (uword i = 0; i < nPar; i++){
    p0 = inputs->p;
    p0.row(i) += inc(i);
    // grad0(i) = inputs->llikFUN(p0, inputs);
    if (inputs->augmented){
      grad0(i) = llikAug(p0, inputs);
    } else {
      grad0(i) = llik(p0, inputs);
    }
  }
  for (uword i = 0; i < nPar; i++){
    for (uword j = i; j < nPar; j++){
      p0 = inputs->p;
      p0.row(i) += inc(i);
      p0.row(j) += inc(j);
      if (inputs->augmented){
        llikValue2 = llikAug(p0, inputs);
      } else {
        llikValue2 = llik(p0, inputs);
      }
      Hess(i, j) = as_scalar((llikValue2 - grad0.row(i) - grad0.row(j) + llikValue0)
                               / inc(i) / inc(j));
    }
  }
  if (nPar > 1){
    Hess = Hess + trimatu(Hess, 1).t();
  }
  return Hess;
}
// True filter/smooth/disturb function
void auxFilter(unsigned int smooth, SSinputs& data){
  // smooth (0: filter, 1: smooth, 2: disturb)
  // double tolsta = 0; //1e-7;
  uword n,
        ns,
        nMiss = 0;
  ///////////
  mat T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt;
  setT(0, data.system, T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt);
  ///////////
  mat Pt,
      Pinft,
      Ft(1, 1),
      Finft(1, 1),
      v2F(1, 1),
      iFt(1, 1); //, oldPt;  //, auxKinf;
  vec at,
      Kt,
      vt(1),
      data_F;  //, auxFinf;
  bool colapsed = false,
       steadyState = false,
       miss = false;
  cube cP,
       Pinf;
  // Initialising variables
  uword ny = data.y.n_elem;
  int k = data.u.n_rows;
  vec Nans(data.h); Nans.fill(datum::nan);
  data.y = join_vert(data.y, Nans);
  n = data.y.n_elem;
  data.d_t = n;
  ns = T.n_cols;
  // RQRt = data.system.R * data.system.Q * data.system.R.t();
  // CHCt = data.system.C * data.system.H * data.system.C.t();
  // Inputs part
  rowvec Dt(n, fill::zeros);
  if (k > 0 && Z.n_rows == 1){
    int nn = data.betaAug.n_rows;
    D = data.betaAug.rows(nn - k, nn - 1).t();
    Dt = D * data.u;
  }
  KFinit(T, RQRt, ns, at, Pt, Pinft);
  data.v = zeros(n);
  data.a = zeros(ns, n);
  if  (smooth > 0){
    cP = zeros(ns, ns, n);
    Pinf = zeros(ns, ns, data.d_t + 1);
  }
  data.P = zeros(ns, n);
  data.yFit = data.v;
  data.K = zeros(ns, n);
  data_F = data.v;
  data.iF = data.v;
  data.Finf = zeros(data.d_t + 1);
  data.Kinf = zeros(ns, data.d_t + 1);
  v2F.fill(0);
  Kt.resize(ns);
  // mat Z = data.system.Z.row(0);
  // bool TVP = false;
  // if (data.system.Z.n_rows > 1)
  //     TVP = true;
  // KF loop
  for (uword t = 0; t < n; t++){
    ////////////////
    if (t > 0 && data.system.tvp){
        setT(t, data.system, T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt);
    }
    /////////////////
    // if (TVP)
    //     Z = data.system.Z.row(t);
    if (!colapsed && smooth > 0){
      Pinf.slice(t) = Pinft;
    }
    data.yFit.row(t) = Z * at + Dt(t);
    // Storing for smoothing/disturbing
    data.v.row(t) = data.y.row(t) - data.yFit.row(t);
    if (smooth > 0){
      data.a.col(t) = at;
      cP.slice(t) = Pt;
      data.P.col(t) = Pt.diag();
      if (steadyState){
        data_F.row(t) = data_F.row(t - 1);
      } else {
        data_F.row(t) = Z * Pt * Z.t() + CHCt;
      }
    }
    // Data missing
    if (!is_finite(data.y.row(t))){
      steadyState = false;
      miss = true;
      nMiss++;
    } else {
      miss = false;
    }
    // Correction
    KFcorrection(miss, colapsed, steadyState, smooth, &data, CHCt,
                 Finft, vt, Dt(t), Ft, iFt, at, Pt, Pinft, Kt, t,
                 data.Finf, data.Kinf, Z);
    if (!miss && t < n){
      if (colapsed || Finft(0, 0) < 1e-8){
        v2F  += vt * iFt * vt;
      }
    }
    // Storing information
    data.v.row(t) = vt;
    data_F.row(t) = Ft;
    data.iF.row(t) = iFt;
    if (smooth > 0){
      data.K.col(t) = Kt;
    } else {
      data.yFit.row(t) = Z * at + Dt(t);
      data.a.col(t) = at;
      data.P.col(t) = Pt.diag();
    }
    // Prediction
    KFprediction(steadyState, colapsed, T, RQRt, at, Pt, Pinft);
    // Checking colapsed
    if (!colapsed && all(all(abs(Pinft) < 1e-6))){
        colapsed = true;
        data.d_t = t;
    }
    // Storing final state and covariance for forecasting
    if (t == ny - 1){
      data.PEnd = Pt;
      data.aEnd = at;
    }
  }
  // Smoothing loop
  data.F = data_F;   // For final normalization of innovations
  if (smooth > 0){
    mat Nt(ns, ns),
        Ninfti(ns, ns),
        N2t(ns, ns),
        PPinf(ns, ns); //, RR(ns, ns), sysmatQ, sysmatR, Z = data->system.Z;
    mat Inew(ns, ns),
        Lt(ns, ns),
        Linft(ns, ns),
        LinftNt(ns, ns),
        Ninft(ns, ns);
    vec rt(ns),
        rinft(ns),
        Kinft(ns),
        Z_Ft(ns),
        Z_Finft(ns); //, eta; //, vt(1), Kt(ns), Ft(1), GammaD(1);
      mat QRt,
          Veta,
          pinvVeta;
    bool colapsed = true;
    ///////////
    // mat T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt;
    setT(n - 1, data.system, T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt);
    ///////////
    if (smooth == 2){   // Disturbance
      int rQ = Q.n_rows; //, cR = R.n_cols;
      QRt = Q * R.t();
      data.eta.zeros(rQ, n - data.h);
      Veta.zeros(rQ, rQ);
    }
    // Storing in case of outlier detection
    if (smooth == 3){
      data.rNrOut = zeros(n);
      data.rOut = zeros(ns, n);
      data.NOut = zeros(ns, ns, n);
    }
    // Initialising variables
    Nt.fill(0);
    Ninft = Nt;
    N2t = Nt;
    rt.fill(0);
    rinft = rt;
    Inew.eye();
    // Main Loop
    int intN = n;
    for (int t = n - 1; t >= 0; t--){
      ////////////////
      if (t < (int)n - 1 && data.system.tvp){
          setT(t, data.system, T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt);
      }
      /////////////////
      // if (TVP)
      //   Z = data.system.Z.row(t);
      if (t <= data.d_t){
        colapsed = false;
      }
      vt = data.v.row(t);
      Kt = data.K.col(t);
      iFt = data.iF.row(t);
      Ft = data.F.row(t);
      if (!colapsed){
        Finft = data.Finf.row(t);
        Kinft = data.Kinf.col(t);
      }
      if (is_finite(data.y.row(t))){
        miss = false;
        if (colapsed || Finft(0, 0) < 1e-8) {
          Lt = T - T * Kt * Z;
          Z_Ft = Z.t() * iFt(0);
          rt = Z_Ft * vt + Lt.t() * rt;
          Nt = Z_Ft * Z + Lt.t() * Nt * Lt;
          if (!colapsed){
              rinft = T.t() * rinft;
              N2t = T.t() * N2t * T;
              Ninft = T.t() * Ninft * Lt;
          }
        } else if (Finft(0, 0) >= 1e-8) {
            Lt = T - T * Kinft * Z;
            Linft = -T * Kt * Z;
            Z_Finft = Z.t() * iFt(0, 0);  //  / Finft(0, 0);
            rinft = Z_Finft * vt + Lt.t() * rinft + Linft.t() * rt;
            rt = Lt.t() * rt;
            LinftNt = Lt.t() * Ninft * Linft;
            N2t = -Z_Finft * Z_Finft.t() * Ft(0, 0) + Linft.t() * Nt * Linft +
                LinftNt + LinftNt.t() + Lt.t() * N2t * Lt;
            Ninft = Z_Finft * Z + Lt.t() * Ninft * Lt + Linft.t() * Nt * Lt;
            Nt = Lt.t() * Nt * Lt;
        }
      } else {
        miss = true;
      }
      Pt = cP.slice(t);
      data.a.col(t) += Pt * rt;
      if (!colapsed){
        Pinft = Pinf.slice(t);
        data.a.col(t) += Pinft * rinft;
      }
      data.yFit.row(t) = Z * data.a.col(t) + Dt(t); // + D;
      Pt -= Pt * Nt * Pt;
      if (!colapsed){
        PPinf = Pinft * Ninft * Pt;
        Pt = Pt - PPinf - PPinf.t() - Pinft * N2t * Pinft;
      }
      data.F.row(t) = Z * Pt * Z.t() + CHCt;
      data.P.col(t) = Pt.diag();
      //Disturbance smoother
      if (smooth == 2 && t < intN - data.h){
        Veta = Q - QRt * Nt * QRt.t();
        data.eta.col(t) = QRt * rt / sqrt(Veta.diag());
      }
      // Storing for outlier detection
      if (smooth == 3){
        data.rNrOut.row(t) = rt.t() * pinv(Nt) * rt;
        data.rOut.col(t) = rt;
        data.NOut.slice(t) = Nt;
      }
      // Passing to rt(t-1) and Nt(t-1)
      if (t > 0 && miss){
        rt = T.t() * rt;
        Nt = T.t() * Nt * T;
        if (!colapsed){
          rinft = T.t() * rinft;
          Ninft = T.t() * Ninft * T;
          N2t = T.t() * N2t * T;
        }
      }
    }
  }
  // Post-processing outputs
  int nTrue;
  if (data.d_t < (int)(data.system.T.n_cols + 10)){
    // Colapsed KF
    nTrue = n - nMiss - 1 - data.d_t;
  } else {
    // KF did not colapse
    nTrue = n - nMiss - 1 - data.system.T.n_cols;
  }
  double innVar = data.innVariance;
  if (data.cLlik){         // Concentrated Likelihood
    innVar = v2F(0, 0) / nTrue;
  }
  data.y = data.y.rows(0, ny - 1);
  data_F *= innVar;
  data.P *= innVar;
  data.FFor *= innVar; // * scale;
  // Cleaning innovations
  if ((uword)data.d_t < n - 10){
    data.v(span(0, data.d_t)).fill(datum::nan);
  } else {
    data.v(span(0, sum(ns) + 1)).fill(datum::nan);
  }
  uvec ind = find_finite(data.v);
  if (ind.n_elem < 5){
      ind = regspace<uvec>(sum(ns), data.v.n_elem - data.h - 1);
  }
  if (smooth == 3){
    data.F = data_F.rows(0, max(ind)); // * data.innVariance;
    data.v = data.v.rows(0, max(ind));
  } else {
    data.F = data_F.rows(min(ind), max(ind));
    data.v = data.v.rows(min(ind), max(ind));
    data.v = data.v / sqrt(data.F);
  }
  // Disturbances
  if (smooth == 2){
    data.eta = data.eta.t();
    mat covEta = diagmat(cov(data.eta, 1));
    uvec ind = find_finite(covEta.diag());
    mat icovEta = pinv(sqrt(covEta(ind, ind)));
    data.eta.cols(ind) = (data.eta.cols(ind) * icovEta);
    data.eta = data.eta.t();
    data.eta.replace(datum::nan, 0);
    data.eta.replace(datum::inf, 0);
  }
}
// solution to lyapunov equation P = Phi * P * Phi' + Q
mat dlyap(mat T, mat Q){
    uword n = T.n_cols;
    mat Ur, Phir;
    schur(Ur, Phir, T);  // U * S * U' = T
    mat ceros = zeros(n, n);
    cx_mat U = cx_mat(Ur, ceros), Phi = cx_mat(Phir, ceros);
    // Pass schur to complex schur
    uvec k(2);
    cx_vec mu(2), r(1), c(1), s(1);
    cx_mat G(2, 2);
    for (uword m = n - 1; m >= 1; --m) {
        if (norm(Phi(m, m - 1)) != 0.0) {
            k(0) = m - 1; k(1) = m;  //regspace<uvec>(m - 1, m);
            mu = eig_gen(Phi(k, k)) - Phi(m, m) * ones<cx_vec>(2);
            // r = std::hypot(mu(0), Phir(m, m - 1));
            r = sqrt(norm(mu(0)) + norm(Phi(m, m - 1)));
            c = mu(0) / r;
            s = Phi(m, m - 1) / r;
            G = join_rows(join_cols(c.t(), -s), join_cols(s, c));
            Phi.submat(m - 1, m - 1, m, n - 1) = G * Phi.submat(m - 1, m - 1, m, n - 1);
            Phi.submat(0, m - 1, m, m) = Phi.submat(0, m - 1, m, m) * G.t();
            U.cols(k) = U.cols(k) * G.t();
            Phi(m, m - 1) = cx_double(0.0, 0.0);
        }
    }
    // finding matrix P
    cx_mat Qc = U.t() * Q * U, P = zeros<cx_mat>(n, n);
    c = cx_double(1.0, 0.0);
    for (int j = n-1; j >= 0; --j) {
        for (int i = n-1; i >= 0; --i) {
            r = c - Phi(i,i) * Phi.submat(j, j, j, j).t();
            s = Phi.submat(i, i, i, n - 1) * P.submat(i, j, n - 1, n - 1) *
                    Phi.submat(j, j, j, n - 1).t() + Qc.submat(i, j, i, j);
            if (norm(s) < pow(2, -52)) {
                P(i, j) = cx_double(0.0, 0.0);
            } else if (norm(r) < pow(2, -52)) {
                throw std::runtime_error("SSpace: Lyapunov equation with no solution!!!");
            } else {
                P(i, j) = s(0) / r(0);
            }
        }
    }
    return real(U * P * U.t());
}
// Check consistency of system matrices in SSmatrix
bool checkSystem(SSmatrix& sys){
    bool errorExit = false;
    sys.Gam = sys.S = 0.0;   // Not implemented
    // Consistency of system matrices
    if (sys.T.n_rows < 1){
        printf("%s", "ERROR: Matrix T should not be empty!!!\n");
        errorExit = true;
        return errorExit;
    } else if (sys.R.n_rows < 1){
        printf("%s", "ERROR: Matrix R should not be empty!!!\n");
        errorExit = true;
        return errorExit;
    } else if (sys.Q.n_rows < 1){
        printf("%s", "ERROR: Matrix Q should not be empty!!!\n");
        errorExit = true;
        return errorExit;
    } else if (sys.Z.n_rows < 1){
        printf("%s", "ERROR: Matrix Z should not be empty!!!\n");
        errorExit = true;
        return errorExit;
    } else if (sys.C.n_rows < 1){
        printf("%s", "ERROR: Matrix C should not be empty!!!\n");
        errorExit = true;
        return errorExit;
    } else if (sys.H.n_rows < 1){
        printf("%s", "ERROR: Matrix H should not be empty!!!\n");
        errorExit = true;
        return errorExit;
    }
    sys.ns = sys.T.n_cols;
    sys.ny = 1;
    sys.n = -1;
    // if (sys.Gam.n_rows < 1){
    //     sys.Gam.resize(sys.ns, 1);
    //     sys.Gam.fill(0.0);
    // }
    if (sys.D.n_rows < 1){
        sys.D.resize(sys.ny, 1);
        sys.D.fill(0.0);
    }
    // if (sys.S.n_rows < 1){
    //     sys.S.resize(sys.Q.n_cols, sys.H.n_cols);
    //     sys.S.fill(0.0);
    // }
    mat T, R, Gam, D, S, Q, Z, C, H, RQRt, CHCt;
    setT(0, sys, T, Gam, R, Q, Z, D, C, H, S, RQRt, CHCt);
    if (T.n_rows != T.n_cols){
        printf("%s", "ERROR: Matrix T should be square!!!\n");
        errorExit = true;
        return errorExit;
    } else if (Q.n_rows != Q.n_cols){
        printf("%s", "ERROR: Matrix Q should be square!!!\n");
        errorExit = true;
        return errorExit;
    } else if (H.n_rows != H.n_cols){
        printf("%s", "ERROR: Matrix H should be square!!!\n");
        errorExit = true;
        return errorExit;
    } else if (T.n_rows != R.n_rows){
        printf("%s", "ERROR: Rows of T and R should be the same!!!\n");
        errorExit = true;
        return errorExit;
    } else if (R.n_cols != Q.n_rows){
        printf("%s", "ERROR: Columns of R and rows of Q should be the same!!!\n");
        errorExit = true;
        return errorExit;
    } else if (Z.n_cols != T.n_rows){
        printf("%s", "ERROR: Columns of Z and rows of T should be the same!!!\n");
        errorExit = true;
        return errorExit;
    } else if (Z.n_rows != C.n_rows){
        printf("%s", "ERROR: Rows of Z and C should be the same!!!\n");
        errorExit = true;
        return errorExit;
    } else if (C.n_cols != H.n_rows){
        printf("%s", "ERROR: Columns of C and rows of H should be the same!!!\n");
        errorExit = true;
        return errorExit;
    }
    // Setting TVP properties of system
    vec tvMatrices(9, fill::zeros);
    int n = -1;
    sys.n = -1;
    if (sys.T.n_rows > sys.ns){
        tvMatrices(0) = 1.0;
        sys.n = sys.T.n_rows / sys.ns;
        if (n == -1)
            n = sys.n;
        else if (n != sys.n){
            printf("%s", "ERROR: Time varying matrices of different time dimension!!!\n");
            errorExit = true;
            return errorExit;
        }
    }
    // if (sys.Gam.n_rows > sys.ns){
    //     tvMatrices(1) = 1.0;
    //     sys.n = sys.Gam.n_rows / sys.ns;
    // }
    if (sys.R.n_rows > sys.ns){
        tvMatrices(2) = 1.0;
        sys.n = sys.R.n_rows / sys.ns;
        if (n == -1)
            n = sys.n;
        else if (n != sys.n){
            printf("%s", "ERROR: Time varying matrices of different time dimension!!!\n");
            errorExit = true;
            return errorExit;
        }
    }
    if (sys.Q.n_rows > sys.ns){
        tvMatrices(3) = 1.0;
        sys.n = sys.Q.n_rows / sys.R.n_cols;
        if (n == -1)
            n = sys.n;
        else if (n != sys.n){
            printf("%s", "ERROR: Time varying matrices of different time dimension!!!\n");
            errorExit = true;
            return errorExit;
        }
    }
    if (sys.Z.n_rows > sys.ny){
        tvMatrices(4) = 1.0;
        sys.n = sys.Z.n_rows;
        if (n == -1)
            n = sys.n;
        else if (n != sys.n){
            printf("%s", "ERROR: Time varying matrices of different time dimension!!!\n");
            errorExit = true;
            return errorExit;
        }
    }
    if (sys.D.n_rows > sys.ny){
        tvMatrices(5) = 1.0;
        sys.n = sys.D.n_rows;
    }
    if (sys.C.n_rows > sys.ny){
        tvMatrices(6) = 1.0;
        sys.n = sys.C.n_rows;
        if (n == -1)
            n = sys.n;
        else if (n != sys.n){
            printf("%s", "ERROR: Time varying matrices of different time dimension!!!\n");
            errorExit = true;
            return errorExit;
        }
    }
    if (sys.H.n_rows > sys.ny){
        tvMatrices(7) = 1.0;
        sys.n = sys.H.n_rows / sys.C.n_cols;
        if (n == -1)
            n = sys.n;
        else if (n != sys.n){
            printf("%s", "ERROR: Time varying matrices of different time dimension!!!\n");
            errorExit = true;
            return errorExit;
        }
    }
    // if (sys.S.n_rows > sys.ny){
    //     tvMatrices(8) = 1.0;
    //     sys.n = sys.S.n_rows;
    // }
    if (sum(tvMatrices > 0.0))
        sys.tvp = true;
    sys.tvMatrices = find(tvMatrices);
    if (tvMatrices(2) == 0 && R.n_rows == R.n_cols && sum(sum((R - eye(R.n_cols, R.n_cols)))) == 0)
        sys.identityR = true;
    if (tvMatrices(6) == 0 && C.n_rows == C.n_cols && sum(sum((C - eye(C.n_cols, C.n_cols)))) == 0)
        sys.identityC = true;
    // if (sys.identityR && tvMatrices(3) == 0)
    //     sys.constantRQR = true;
    // if (sys.identityC && tvMatrices(7) == 0)
    //     sys.constantCHC = true;
    return errorExit;
}
// Set time varying matrices to time t
void setT(uword t, SSmatrix& sys, mat& T, mat& Gam, mat& R, mat& Q, mat& Z, mat& D, mat& C, mat& H, mat& S, mat& RQRt, mat& CHCt){
    uword i0 = t * sys.ns, i1 = (t + 1) * sys.ns - 1;
    if (t == 0){
        T = sys.T.rows(i0, i1);
        // Gam = sys.Gam.rows(i0, i1);
        R = sys.R.rows(i0, i1);
        i1 = (t + 1) * sys.R.n_cols - 1;
        Q = sys.Q.rows(i0, i1);
        Z = sys.Z.row(t);
        D = sys.D.row(t);
        C = sys.C.row(t);
        i1 = (t + 1) * sys.C.n_cols - 1;
        H = sys.H.rows(i0, i1);
        // S = sys.S.row(t);
    } else if (sys.tvp){
        for (uword i = 0; i < sys.tvMatrices.n_elem; i++){
            if (sys.tvMatrices(i) == 0)
                T = sys.T.rows(i0, i1);
            // else if (sys.tvMatrices(i) == 1)
            //    Gam = sys.Gam.rows(ind);
            else if (sys.tvMatrices(i) == 2)
                R = sys.R.rows(i0, i1);
            else if (sys.tvMatrices(i) == 3){
                i0 = t * sys.R.n_cols;
                i1 = (t + 1) * sys.R.n_cols - 1;
                Q = sys.Q.rows(i0, i1);
            }
            else if (sys.tvMatrices(i) == 4)
                Z = sys.Z.row(t);
            else if (sys.tvMatrices(i) == 5)
                D = sys.D.row(t);
            else if (sys.tvMatrices(i) == 6)
                C = sys.C.row(t);
            else if (sys.tvMatrices(i) == 7){
                i0 = t * sys.C.n_cols;
                i1 = (t + 1) * sys.C.n_cols - 1;
                H = sys.H.rows(i0, i1);
            }
            // else if (sys.tvMatrices(i) == 8)
            //     S = sys.S.rows(ind);
        }
    }
    if (!sys.identityR)
        RQRt = R * Q * R.t();
    else
        RQRt = Q;
    if (!sys.identityC)
        CHCt = C * H * C.t();
    else
        CHCt = H;
}
