#include <TMB.hpp>

extern "C" void Rf_rmultinom(int, double*, int, int*);

// Multinomial simulator. Move to TMB.
template<class Type>
vector<Type> rmultinom(Type size, vector<Type> prob) {
  int K = prob.size();
  vector<int> rN(K);
  int n = (int) asDouble(size);
  vector<double> p(K);
  for (int i=0; i<K; i++) p[i] = asDouble(prob(i));
  p = p / p.sum();
  Rf_rmultinom(n, p.data(), K, rN.data());
  vector<Type> ans(K);
  for (int i=0; i<K; i++) ans[i] = rN[i];
  return ans;
}

template<class Type>
vector<int> sample_replacement(int size, vector<Type> prob) {
  int K = prob.size();
  vector<Type> N = rmultinom(Type(size), prob);
  vector<int> ans(size);
  int k = 0;
  for (int i=0; i<K; i++) {
    for (int j=0; j<asDouble(N(i)); j++) {
      ans(k) = i; k++;
    }
  }
  return ans;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_SPARSE_MATRIX(Q0);
  DATA_SPARSE_MATRIX(I);
  
  
  /* Fishery dependent data */
  DATA_FACTOR(time);      // Length = haulid length
  DATA_FACTOR(gf);        // Length = haulid length
  DATA_FACTOR(rowid);     // Length = haulid length; haulid matching the rows of the dataframe (rowid)
  DATA_VECTOR(response);  // "shorter"
  DATA_MATRIX(X);         // nrow = response length (Both fixed and random)
  DATA_VECTOR(offset);
  //DATA_IVECTOR(fd_indicator); //Length=response length
  //DATA_IVECTOR(support_area); //Length=nlevels gf
    
  /* Random fields */
  PARAMETER_ARRAY(eta_density);
  PARAMETER_VECTOR(eta_nugget); // Not used ! perhaps later.
  
  /* 5 common fixed effects x number of times */
  PARAMETER(logdelta);        // Length=2
  PARAMETER(logscale);         // Dim = c(ntime,2)
  PARAMETER(logsd_nugget);    // Length = ntime. not used
  PARAMETER(time_corr);
  
  
  /* Parameters */
  PARAMETER_VECTOR(beta);         // Fixed effects
  PARAMETER_VECTOR(beta_r);       // Random effects
  PARAMETER_VECTOR(beta_r_logsd);
  DATA_FACTOR(beta_r_fac);
  PARAMETER_VECTOR(logphi); // NB overdispersion
  PARAMETER_VECTOR(alpha);
  
  /* Stuff for prediction */
  DATA_INTEGER(doPredict);  // Flag
  DATA_MATRIX(Xpredict);
  DATA_FACTOR(Apredict);

  /* logphi */
  DATA_FACTOR(Data);

  /* Flat prior for 'robustification' */
  DATA_SCALAR(huge_sd);

  /* Distance in KM between grid points  */
  DATA_SCALAR(h);

  /* Index calc only: which betas are time effects? */
  DATA_IVECTOR(which_beta_time);

  Type sd_nugget = exp(logsd_nugget);
  Type ans = 0;
  using namespace density;

  /* Add *flat* prior to fixed effects */
  ans -= dnorm(beta, Type(0), huge_sd, true).sum();

  /* Add random effects beta_r */
  for(int i=0; i<beta_r.size(); i++) {
    ans -= dnorm(beta_r(i), Type(0), exp(beta_r_logsd(beta_r_fac(i))), true);
    SIMULATE {
      beta_r(i) = rnorm(Type(0), exp(beta_r_logsd(beta_r_fac(i))));
    }
  }
  vector<Type> beta_full(X.cols());
  beta_full << beta, beta_r;
  if(offset.size() == 0) {
    offset.resize(X.rows());
    offset.setZero();
  }

    /* Optional: Add static field */
  PARAMETER_VECTOR(eta_static);
  PARAMETER_VECTOR(logdelta_static);
  PARAMETER_VECTOR(logscale_static);

  // NOTE: eta_static code must be revised !
  if(eta_static.size() > 0) {
    /* Scale parameters for fields */
    Type scale = exp(logscale_static[0]);
    /* GMRF: Sparse precision matrix */
    Eigen::SparseMatrix<Type> Q = Q0 + exp(logdelta_static[0]) * I;
    GMRF_t<Type> nldens = GMRF(Q);
    ans += SCALE(nldens, scale)(eta_static);
    /* Simulate static field */
    SIMULATE {
      SCALE(nldens, scale).simulate(eta_static);
    }
    int ntimes = NLEVELS(time);
    for(int i=0; i<ntimes; i++) {
       eta_density.col(i) -= eta_static;
    }
    REPORT(ntimes);
  }

  /* Add static covariates */
  DATA_MATRIX(Xs);
  PARAMETER_VECTOR(beta_s);
  vector<Type> mu_static = Xs * beta_s;
  if (beta_s.size()) {
    int ntimes = NLEVELS(time);
    for(int i=0; i<ntimes; i++) {
       eta_density.col(i) -= mu_static;
    }
  }

  /* Time covariance */
  //N01<Type> nldens_time;
  Type phi = time_corr / sqrt(1.0 + time_corr*time_corr);
  AR1_t<N01<Type> > nldens_time = AR1(phi);
  /* Scale parameters for fields */
  Type scale = exp(logscale);
  /* GMRF: Sparse precision matrix */
  Eigen::SparseMatrix<Type> Q = Q0 + exp(logdelta) * I;
  GMRF_t<Type> nldens = GMRF(Q);
  ans += SEPARABLE(SCALE(nldens_time, scale), nldens)(eta_density);
  /* Simulate space-time random field (and remember to add static) */
  SIMULATE {
    SEPARABLE(SCALE(nldens_time, scale), nldens).simulate(eta_density);
  }
  if (beta_s.size()) {
    int ntimes = NLEVELS(time);
    for(int i=0; i<ntimes; i++) {
       eta_density.col(i) += mu_static;
    }
  }
  if(eta_static.size() > 0) {
    int ntimes = NLEVELS(time);
    for(int i=0; i<ntimes; i++) {
      eta_density.col(i) += eta_static;
    }
  }

  /* Nugget */
  if(eta_nugget.size() > 0) {
    ans -= dnorm(eta_nugget, Type(0), sd_nugget, true).sum();
    SIMULATE {
      eta_nugget = rnorm(eta_nugget.size(), Type(0), sd_nugget);
    }
  }

  /* Include preferential sampling likelihood */
  DATA_FACTOR(SupportAreaGroup);
  DATA_IARRAY(SupportAreaMatrix);
  matrix<Type> probAreaMatrix(NLEVELS(gf),
                              NLEVELS(time) ); // For simulation only
  probAreaMatrix.setZero();
  for (int group = 0; group < NLEVELS(SupportAreaGroup); group++) {
    vector<int> support_area = SupportAreaMatrix.col(group);
    vector<Type> logsum (NLEVELS(time));
    logsum.setZero(); logsum = log(logsum); // logsum = -INFINITY
    for(int j=0; j<logsum.size();j++){
      for(int i=0; i<NLEVELS(gf); i++){
        if(support_area(i)) logsum(j) = logspace_add(logsum(j), alpha[group]*eta_density(i,j));
      }
    }
    // Given 'group' calculate space-time probabilities of trawl positions (used for simulation only)
    SIMULATE {
      probAreaMatrix.setZero();
      for(int j=0; j<NLEVELS(time); j++) {
        for (int i=0; i<NLEVELS(gf); i++) {
          if(support_area(i))
            probAreaMatrix(i, j) = exp( alpha[group] * eta_density(i, j) - logsum( j ) );
        }
      }
    }
    for(int i=0; i<rowid.size(); i++){
      if(i==0 || (rowid(i)!=rowid(i-1))){
        int pos = gf(i);
        int tim = time(i);
        // density ~ lambda^alpha
        if (SupportAreaGroup(rowid(i)) == group) {
          ans -= alpha[group]*eta_density(pos,tim) - logsum(tim);
          SIMULATE {
            vector<Type> tmp = probAreaMatrix.col(tim).array();
            gf(i) = sample_replacement(1, tmp)[0];
            for (int k = i; k<rowid.size() && (rowid(k)==rowid(i)); k++) {
              gf(k) = gf(i);
            }
          }
        }
      }
    }
  }
  
  /* Fishery dep data */
  if(response.size() > 0) {
    vector<Type> log_mu(response.size());
    log_mu.setZero(); log_mu = log(log_mu); // log_mu = -INFINITY
    for(int i=0; i<gf.size(); i++) {
      if(eta_static.size() == 0) {
        log_mu(rowid(i)) = logspace_add(log_mu(rowid(i)),
                                        eta_density(gf[i], time[i]) );
      } else {
        log_mu(rowid(i)) = logspace_add(log_mu(rowid(i)),
                                        eta_density(gf[i], time[i])  + eta_static(gf[i]) );
      }
    }
    log_mu = log_mu + X * beta_full + offset;
    /* Use numerical robust negative binomial distribution rather than
       mu = exp(log_mu);
       var = mu + mu*mu / exp(logphi(0));
       ans -= dnbinom2(response, mu, var, true).sum(); */
    vector<Type> log_var_minus_mu = log_mu + log_mu;
    for (int i=0; i < response.size(); i++) log_var_minus_mu(i) -= logphi(Data(i));
    REPORT(log_mu);            // For debugging
    REPORT(log_var_minus_mu);  // For debugging
    ans -= dnbinom_robust(response, log_mu, log_var_minus_mu, true).sum();
    SIMULATE {
      vector<Type> mu = exp(log_mu);
      vector<Type> var = exp(log_var_minus_mu) + mu;
      response = rnbinom2(mu, var);
    }
  }

  SIMULATE {
    // Simulated random effects
    REPORT(eta_density);
    REPORT(eta_static);
    REPORT(beta_r);
    // Simulated data
    REPORT(response);
    REPORT(gf); // FIXME: Positions not yet sampled
    // TMB problem: This woodo is needed to keep the factor levels for the simulated gf.
    Rf_setAttrib( findVar(install("gf"), this -> report),
                  install("levels"),
                  Rf_getAttrib( getListElement(this -> data, "gf" ) , install("levels") ) );
  }

  // Do index calculation ?
  if ( which_beta_time.size() > 0 ) {
    int ntimes = NLEVELS(time);
    vector<Type> logIndex(ntimes);
    logIndex.fill( Type( R_NaReal ) );
    for (int i = 0; i < ntimes; i++) {
      if (which_beta_time[i] != -1) { // -1 == NA !!!
        logIndex(i) =
          beta_full[ which_beta_time[i] ] +
          log( eta_density.col(i).exp().mean() );
      }
    }
    ADREPORT(logIndex);
  }

  if(doPredict){
    array<Type> logindex = eta_density;
    // NOTE: eta_density is the *total* field (including the static if present)
    //         ^^^___ RIGHT !!!!
    // if(eta_static.size() > 0) {
    //  for(int j=0; j<logindex.cols(); j++)
    //    logindex.col(j) += eta_static;
    // }
    // if(mu_static.size() > 0) {
    //  for(int j=0; j<logindex.cols(); j++)
    //    logindex.col(j) += mu_static;
    // }    
    // FIXME: Add common fixed effects to predictions (incude covariates available at every grid point).
    REPORT(logindex);
    ADREPORT(logindex);
    
    /* Note: One can use ADREPORT rather than REPORT to get uncertainties but it may require too much memory */
    if(isDouble<Type>::value) {
      Type avgVariance = nldens.variance().mean() * pow(scale, 2.); // This represents the first term in eq.3 in Kristensen et al. (2013)
      Type delta = exp(logdelta);
      // Spatial decorrelation distance (corr=exp(-1)) - to compute a plot of spatial correlation vs. distance (km)
      Type H = h / log(1 + (delta/2) + sqrt(delta + (delta*delta/4)));
      Type SpatialSD = sqrt(avgVariance); // Convert to standard deviation
      Type ForecastSD = sqrt(1. - phi*phi) * SpatialSD; // High time corr ==> Better ability to forcast
      REPORT(H);
      REPORT(SpatialSD);
      REPORT(ForecastSD);
    }
  }

  return ans;
}
