functions {
  vector linear(
    real t,        // time
    vector y,      // state
    matrix A,
    array[] matrix B,
    matrix C,
    vector tp_vec,
    matrix u){  
    
    int m = num_elements(y);
    int n_u = cols(C);
    
    vector[m] dydt;
    vector[n_u] u_t;
    matrix[m,m] B_all = rep_matrix(0,m,m);
    
    for(i in 1:n_u){
      u_t[i] = constant_interpolation(t,tp_vec,u[,i]);
      B_all += u_t[i]*B[i];
    }

    dydt = (A+B_all)*y + C*u_t;
    return dydt;
  }
  
  vector convolve(vector a, vector b) {
    int na = num_elements(a);       // Vector lengths
    int nb = num_elements(b);
    int n_zero_a = nb - 1;          // Zero padding lengths
    int n_zero_b = na - 1;

    vector[nb] b_rev = reverse(b); // The reversed b vector

    vector[na + n_zero_a] a_pad;    // Instantiate zero padded vectors
    vector[nb + n_zero_b] b_rev_pad;

    // Perform zero padding
    a_pad[1 : n_zero_a] = rep_vector(0, n_zero_a);
    b_rev_pad[(nb + 1) : (nb + n_zero_b)] = rep_vector(0, n_zero_b);
    
    // Fill in padded vector
    a_pad[(n_zero_a + 1) : (na + n_zero_a)] = a;
    b_rev_pad[1 : nb] = b_rev;

     // Stan automatically scales the inverse FFT, so no rescaling needed
    return get_real(inv_fft(fft(a_pad) .* conj(fft(b_rev_pad))));
  }
  
  // Function for evaluating numeric solution
  array[] vector numeric_sol(int T, int m, vector z0, array[] real tp, data real rel_tol,
                             data real abs_tol, int max_num_steps, matrix A, array[] matrix B,
                             matrix C, matrix u, vector tp_vec){
    array[T] vector[m] ans;
    ans = ode_ckrk_tol(linear, z0, 0, tp, rel_tol, abs_tol, max_num_steps, A,B,C,tp_vec,u);
    return ans;
    
    
  }
    
  vector stan_convolve(int nPoints, vector a,  vector b) {
    vector [nPoints] out;
    vector [nPoints] b_rev;
    
    b_rev = reverse(b);
    
    for(i in 1:nPoints)  
      out[i] = dot_product(head(a, i), tail(b_rev, i));
   
    return(out);
  }
  
  array[] vector HRF_mu(array[] vector mu, vector tp)
  {
    int m = num_elements(mu[1,]); // number of states
    int T = num_elements(mu[,1]); // number of time points
    
    
    vector[T] HRF_tp; 
    array[T] vector[m] HRF_mu_v;
    
    // HRF(t) evaluated at t in ts
    for (t in 1:T){
      if(tp[t] == 0){
        HRF_tp[t] = 0;
      }else{
        HRF_tp[t] = exp(gamma_lpdf(tp[t] | 6, 1)) - 0.166667*exp(gamma_lpdf(tp[t] | 16, 1));
      }
    }

    for (i in 1:m){
      HRF_mu_v[,i] = to_array_1d(stan_convolve(T, to_vector(mu[,i]), HRF_tp));
    }
    return HRF_mu_v;
  }
  
  real linear_interpolation(real x0, vector x, vector y){
    // assume x is sorted
    real ans;
    int N = num_elements(x);
    int i;
    // return boundary values
    if(x0 <= x[1]){ans = y[1];
    }else if(x0 >= x[N]){ans = y[N];
    }else{
      vector[N] deltas = x0 - x; // positive -> negative 
      i=0;
      for(n in 1:N){
        if(deltas[n]>=0){i+=1;}else{break;}
      }
      
      //print(i);
      
      ans = y[i] + (y[i+1]-y[i])/(x[i+1]-x[i])* (x0-x[i]);
    }
    
    
    return ans;
  }
  
  real constant_interpolation(real x0, vector x, vector y){
    // assume x is sorted
    real ans;
    int N = num_elements(x);
    int i;
    // return boundary values
    if(x0 <= x[1]){ans = y[1];
    }else if(x0 >= x[N]){ans = y[N];
    }else{
      vector[N] deltas = x0 - x; // positive -> negative 
      i=0;
      for(n in 1:N){
        if(deltas[n]>=0){i+=1;}else{break;}
      }
      
      //print(i);
      
      ans = y[i];
    }
    
    
    return ans;
  }
  
  // Function for evaluating ODE solution under scenario A
  vector eval_u_a(int i, vector y, matrix A, vector tp_vec, 
                      int prev_index, int m){
    matrix[m,m] A0 = A;
    vector[m] result;  
    result = (matrix_exp(A0*(tp_vec[i]-tp_vec[prev_index]))*y);
    return result;
  }
  
  // Function for evaluating ODE solution under scenario B
  vector eval_u_b(int i, vector y, matrix A, array[] matrix B, matrix C,
                      vector tp_vec, matrix u, int prev_index, int m){
    matrix[m,m] A0;
    vector[m] b0;
    vector[m] ystar;
    int n_u = cols(C);
    matrix[m,m] B_all = rep_matrix(0,m,m);
    vector[m] result;
    for(j in 1:n_u){
      B_all += B[j]*u[i-1,j];
    }
    A0 = A + B_all;
    b0 = C*(u[i-1,]');
    ystar = -(inverse(A0)*b0);
    result = (matrix_exp(A0*(tp_vec[i]-tp_vec[prev_index]))*(y-ystar)+ystar);
    return result;
  }
  
  // Function for evaluating piecewise analytic solution
  array[] vector analytic_sol(vector z0, matrix A, array[] matrix B, matrix C, 
                      vector tp_vec, matrix u, int m, int n_u, int n_changes,
                      array[] int change_pts){
  
    array[num_elements(tp_vec)] vector[m] ans;
    array[num_elements(tp_vec)+1] vector[m] sol;
    sol[1] = z0;
    
    vector[num_elements(tp_vec)+1] tp_vec_pad;
    tp_vec_pad = append_row(0,tp_vec);
    
    matrix[num_elements(tp_vec)+1,n_u] u_pad;
    u_pad = append_row(rep_row_vector(0,n_u),u);
  
    int prev_index = 1;
    int index;
  
    for (i in 1:n_changes){
  
      index = change_pts[i];
  
      if (sum(u_pad[index-1,]) == 0){
        for (j in (prev_index+1):index) {
          sol[j] = eval_u_a(j,sol[prev_index],A,tp_vec_pad,prev_index,m);
        }
          } else {
              for (j in (prev_index+1):index) {
                sol[j] = eval_u_b(j,sol[prev_index],A,B,C,tp_vec_pad,u_pad,prev_index,m);
              }
            }
  
      prev_index = index; 
      
      
  
    }
    ans = tail(sol, num_elements(tp_vec));
    return ans;
  }
  
}


data {
  int T; // total time points
  int m; // total number of states
  int n_u; // total number of inputs
  int n_changes; // total number of change-points
  array[n_changes] int change_pts; // integer array of index change-points
  int d_A; // total number of parameters in A
  int d_B; // total number of parameters in B
  int d_C; // total number of parameters in C
  array[d_A,2] int A_idxs;
  array[d_B,3] int B_idxs;
  array[d_C,2] int C_idxs;
  array[T] real<lower = 0> tp; // sampled time points
  array[T] vector[m] y_obs; // observed vectors
  matrix[T,n_u] u; // input vectors
  int<lower = 0, upper = 1> conv; // whether the input is convolved with HRF
  real<lower = 1e-3> sigma_nu;
  real<lower = 1e-3> sigma_nu_self;
  real<lower = 1e-3> sigma_z0;
  real<lower = 1e-3> rate_sigma;
  real rel_tol;
  real abs_tol;
  int max_num_steps;
  int<lower =0, upper =1> ode_solver_type;
}

transformed data {
  vector[T] tp_vec = to_vector(tp);
}

parameters {
  array[m] real<lower=1e-3> sigma; // noise
  vector[d_A] nu_A; //parameter in A
  vector[d_B] nu_B; //parameter in A
  vector[d_C] nu_C; //parameter in C
  vector<lower=1e-3>[m] z0;    // initial condition to be inferred
  vector[m] beta;
}

transformed parameters { 
  // vector[m2] theta; // parameter
  // https://mc-stan.org/docs/functions-reference/diagonal-matrix-functions.html
  
}

model {
  matrix[m,m] A = rep_matrix(0,m,m);
  array[n_u] matrix[m,m] B;
  matrix[m,n_u] C = rep_matrix(0,m,n_u);
  
  profile("parameter"){
    // A matrix
  for(i in 1:d_A){
    A[A_idxs[i,1], A_idxs[i,2]] = nu_A[i];
  }
  for (k in 1:m) { 
      A[k, k] = -0.5 * exp(A[k, k]);
    }

  // B matrix
  for(k in 1:n_u){
    // initialize each B matrix
    B[k] = rep_matrix(0,m,m);
  }
  for(i in 1:d_B){
    int k = B_idxs[i,1];
    B[k][B_idxs[i,2],B_idxs[i,3]] = nu_B[i];
  }
  
  // C matrix
  for(i in 1:d_C){
    C[C_idxs[i,1], C_idxs[i,2]] = nu_C[i];
  }
  }
  
  
  array[T] vector[m] mu;
  profile("ode"){
    //function linear
    //vector initial state z0
    //real initial_time 0
    //array[] real times tp
    if(ode_solver_type==0){
      // CKRK algorithm, a 4th/5th order Runge-Kutta method
      mu = ode_ckrk_tol(linear, z0, 0, tp, rel_tol, abs_tol, max_num_steps, A,B,C,tp_vec,u);
    }else if(ode_solver_type==1){
      //Analytic solution
      mu = analytic_sol(z0, A, B, C, tp_vec, u, m, n_u, n_changes, change_pts);
    }

  }
  array[T] vector[m] HRF_mu_v;
  
  //print("mu_1:", mu[,1]);
  //print("HRF_mu_1", HRF_mu_v[,1]);
  if(conv==0){
    HRF_mu_v = mu;
  }else if(conv==1){
    profile("HRV_conv"){
      HRF_mu_v = HRF_mu(mu, tp_vec);
    }
    
  }
  
  
  // prior on parameters
  profile("prior"){
    
  sigma ~ exponential(rate_sigma);
  nu_A ~ normal(0,sigma_nu);
  
  // A: Separate priors based on self connections or not
  for(i in 1:d_A){
    int r = A_idxs[i,1];
    int c = A_idxs[i,2];
    if(r == c){
      // Diagonal (self-connection)
      nu_A[i] ~ normal(0, sigma_nu_self);  // tight prior
    } else {
      // Off-diagonal
      nu_A[i] ~ normal(0, sigma_nu); 
    }
  }
  // B: Separate priors based on self connections or not
  for(i in 1:d_B){
  int k = B_idxs[i,1];
  int r = B_idxs[i,2];
  int c = B_idxs[i,3];
  if(r == c){
    nu_B[i] ~ normal(0, sigma_nu_self);  // tight prior
  } else {
    nu_B[i] ~ normal(0, sigma_nu);
  }
  }
  
  nu_C ~ normal(0,sigma_nu);
  z0 ~ normal(0,sigma_z0);
  beta ~ normal(0, 1);
   
  }
  profile("likelihood"){
  // likelihood
  
  for (t in 1:T) {
      y_obs[t] ~ normal(HRF_mu_v[t] + beta, sigma);
  }
  }
}


 





