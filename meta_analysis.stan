data {
  int<lower=1> K;          // number of studies
  int<lower=1> p;          // number of parameters per study (10)
  int<lower=1> q;          // number of study-level covariates (3)
  array[K] vector[p] y;          // study-level parameter means
  array[K] matrix[p,p] S;  // within-study covariance matrices
  matrix[K,q] X;           // study-level covariates
}

parameters {
  vector[p] alpha;                        // intercepts
  matrix[p,q] B;                          // regression slopes
  vector<lower=0>[p] tau;                 // between-study SDs
  cholesky_factor_corr[p] Lcorr;          // Cholesky factor of correlation
}

transformed parameters {
  matrix[p,p] Ltau = diag_pre_multiply(tau, Lcorr);
  matrix[p,p] Tau = multiply_lower_tri_self_transpose(Ltau);  // Tau = Ltau * Ltau'
}

model {
  // Priors
  alpha ~ normal(0, 1);
  to_vector(B) ~ normal(0, 1);
  tau ~ normal(0, 1);
  Lcorr ~ lkj_corr_cholesky(2);

  for (i in 1:K) {
    vector[p] mu_i = alpha + B * to_vector(X[i,]);

    // compute Cholesky factor
    matrix[p,p] L_i = cholesky_decompose(Tau + S[i]);

    // use Cholesky version of multivariate normal
    target += multi_normal_cholesky_lpdf(y[i] | mu_i, L_i);
}


}
