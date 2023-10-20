data {
    int N; // number of observations
    vector[N] boost; // boost for each observation (response)
    vector[N] day; // boost for each observation (response)
    int N_ind;
    array[N] int ind;
    int N_bt;
    array[N_bt] real bt_vals;
    array[N_ind] int titre_i;

    array[N_ind] int vac_hist; // vac history matrix
    int n_vac_hist; // # of years vac history 
    array[N_ind] int site;
    int n_site;
    array[N_ind] int sex;
    int n_sex;
    array[N_ind] int age_cat;
    int n_age_cat;

    int n_study;
    array[N_ind] int study;

    matrix[n_study, n_site] prop_s;
    matrix[n_study, n_sex] prop_g;
    matrix[n_study, n_vac_hist] prop_v;
    matrix[n_study, n_age_cat] prop_d;
    matrix[n_study, N_bt] prop_t;
    vector[n_study] prop_study;

    vector[n_site] prop_s_all;
    vector[n_sex] prop_g_all;
    vector[n_vac_hist] prop_v_all;
    vector[n_age_cat] prop_d_all;
    vector[N_bt] prop_t_all;

    // include gender
    int<lower=0, upper=1> ind_g;
    // include age
    int<lower=0, upper=1> ind_d;
    // include vaccine history
    int<lower=0, upper=1> ind_vh;
    // include sites
    int<lower=0, upper=1> ind_s;
    // include individual-level effects
    int<lower=0, upper=1> ind_ind;

    int<lower=1, upper=2> likelihood_type; // 1 = normal, 2 = cauchy, 

}
transformed data {
    real delta = 1e-10;
}
parameters {
    // study-speific population-level parameters
    vector[n_study] boost_peak; // boost at peak
    vector<upper = 0, lower = -0.01>[n_study] wane_s; // gradient of wane

    // GP parameters for titre dependent boosting
    real<lower = 0> alpha_b;
    real<lower = 0> rho_b;
    vector[N_bt] eta_b;

    // GP parameters for titre dependent waning
    real<lower = 0, upper = 0.01> alpha_w;
    real<lower = 0> rho_w;
    vector[N_bt] eta_w;

    //--- parameters for 
    // Covariate parameters
    array[ind_vh] vector[n_vac_hist] z_vh_cat;
    array[ind_s] vector[n_site] z_s_cat;
    array[ind_g] vector[n_sex] z_g_cat;
    array[ind_d] vector[n_age_cat] z_d_cat;
    array[ind_vh] real<lower = 0> sigma_vh;
    array[ind_s] real<lower = 0>  sigma_s;
    array[ind_g] real<lower = 0> sigma_g;
    array[ind_d] real<lower = 0>  sigma_d;

    // Covariate parameters interaction terms
    array[ind_vh] vector[n_vac_hist] z_vh_cat_int;
    array[ind_s] vector[n_site] z_s_cat_int;
    array[ind_g] vector[n_sex] z_g_cat_int;
    array[ind_d] vector[n_age_cat] z_d_cat_int;
    array[ind_vh] real<lower = 0> sigma_vh_int;
    array[ind_s] real<lower = 0>  sigma_s_int;
    array[ind_g] real<lower = 0> sigma_g_int;
    array[ind_d] real<lower = 0>  sigma_d_int;

    // Individual-level parameters
    array[ind_ind] vector[N_ind] z_boost_ind;
    array[ind_ind] vector[N_ind] z_wane_ind;
    array[ind_ind] real<lower = 0> sigma_b_ind;
    array[ind_ind] real<lower = 0> sigma_w_ind;

    real<lower = 0> sigma_gp;
}

transformed parameters{

    vector[N_bt] boost_td;
    matrix[N_bt, N_bt] K_b;
    matrix[N_bt, N_bt] L_K_b;
    K_b = gp_exp_quad_cov(bt_vals, alpha_b, rho_b) + diag_matrix(rep_vector(delta, N_bt));
    L_K_b = cholesky_decompose(K_b);
    boost_td = L_K_b * eta_b;

    vector[N_bt] wane_td;
    matrix[N_bt, N_bt] K_w;
    matrix[N_bt, N_bt] L_K_w;
    K_w = gp_exp_quad_cov(bt_vals, alpha_w, rho_w) + diag_matrix(rep_vector(delta, N_bt));
    L_K_w = cholesky_decompose(K_w);
    wane_td = L_K_w * eta_w;

    //vector[N_bt] sigma_gp;
    //matrix[N_bt, N_bt] K_sigma;
   // matrix[N_bt, N_bt] L_K_sigma;
   // K_sigma = gp_exp_quad_cov(bt_vals, alpha_sigma, rho_sigma) + diag_matrix(rep_vector(delta, N_bt));
   // L_K_sigma = cholesky_decompose(K_sigma);
   // sigma_gp = L_K_sigma * eta_sigma;
 
    vector[N_ind] mu_t;
    vector[N_ind] mu_wane;

    vector[n_vac_hist] mu_b;
    vector[n_site] mu_s;
    vector[n_sex] mu_g;
    vector[n_age_cat] mu_d;
    mu_b = rep_vector(0, n_vac_hist);
    mu_s = rep_vector(0, n_site);
    mu_g = rep_vector(0, n_sex);
    mu_d = rep_vector(0, n_age_cat);

    vector[n_vac_hist] mu_b_int;
    vector[n_site] mu_s_int;
    vector[n_sex] mu_g_int;
    vector[n_age_cat] mu_d_int;
    mu_b_int = rep_vector(0, n_vac_hist);
    mu_s_int = rep_vector(0, n_site);
    mu_g_int = rep_vector(0, n_sex);
    mu_d_int = rep_vector(0, n_age_cat);

    vector[N] titre_est;
    vector[N_ind] boost_ind;
    vector[N_ind] wane_ind;

    // Get the estimated titre landscape given the parameters
    // 
    if(ind_vh > 0) {
        mu_b = z_vh_cat[ind_vh] * sigma_vh[ind_vh]; 
        mu_b_int = z_vh_cat_int[ind_vh] * sigma_vh_int[ind_vh]; 
    } 
    if(ind_d > 0) {
        mu_d = z_d_cat[ind_vh] * sigma_d[ind_vh]; 
        mu_d_int = z_d_cat_int[ind_vh] * sigma_d_int[ind_vh]; 
    } 
    if(ind_s > 0) {
        mu_s = z_s_cat[ind_s] * sigma_s[ind_s]; 
        mu_s_int = z_s_cat_int[ind_s] * sigma_s_int[ind_s]; 
    } 
    if(ind_g > 0) {
        mu_g = z_g_cat[ind_g] * sigma_g[ind_g]; 
        mu_g_int = z_g_cat_int[ind_g] * sigma_g_int[ind_g]; 
    } 


    for (i in 1:N_ind) {
        wane_ind[i] = wane_s[study[i]] + wane_td[titre_i[i]];
        boost_ind[i] = (boost_peak[study[i]] + boost_td[titre_i[i]]) * (1 + mu_b_int[vac_hist[i]] + 
            mu_s_int[site[i]] + mu_g_int[sex[i]] + mu_d_int[age_cat[i]]) + 
            mu_g[sex[i]] + mu_b[vac_hist[i]] + mu_s[site[i]] + mu_d[age_cat[i]];
    }
      /*  if(ind_vh > 0) {
            mu_b[i] = z_vh_cat[ind_vh, vac_hist[i]] * sigma_vh[ind_vh]; 
            mu_b_int[i] = z_vh_cat_int[ind_vh, vac_hist[i]] * sigma_vh_int[ind_vh]; 
        } 
        if(ind_d > 0) {
            mu_d[i] = z_d_cat[ind_d, age_cat[i]] * sigma_d[ind_d]; 
            mu_d_int[i] = z_d_cat_int[ind_d, age_cat[i]] * sigma_d_int[ind_d]; 
        } 
        if(ind_s > 0) {
            mu_s[i] = z_s_cat[ind_s, site[i]] * sigma_s[ind_s]; 
            mu_s_int[i] = z_s_cat_int[ind_s, site[i]] * sigma_s_int[ind_s]; 
        } 
        if(ind_g > 0) {
            mu_g[i] = z_g_cat[ind_g, sex[i]] * sigma_g[ind_g]; 
            mu_g_int[i] = z_g_cat_int[ind_g, sex[i]] * sigma_g_int[ind_g]; 
        } 
        boost_ind[i] = mu_t[i] * (1 + mu_b_int[i] +  mu_s_int[i] + mu_g_int[i] + mu_d_int[i]) + 
            mu_g[i] + mu_b[i] + mu_s[i] + mu_d[i];
        wane_ind[i] = mu_wane[i];*/

    if (ind_ind > 0) {
        boost_ind = boost_ind + z_boost_ind[1] * sigma_b_ind[1];
        wane_ind = wane_ind + z_wane_ind[1] * sigma_w_ind[1];
    }

    for (i in 1:N) {
        titre_est[i] = boost_ind[ind[i]] + day[i] * wane_ind[ind[i]];
    }
}
model {
    // likelihood
    real mu;

    for (i in 1:N) {
        if (likelihood_type == 1) {
            mu = normal_cdf( boost[i] + 1 , titre_est[i], sigma_gp) - normal_cdf( boost[i] - 1 , titre_est[i], sigma_gp);
        } else if (likelihood_type == 2) {
            mu = cauchy_cdf( boost[i] + 1 , titre_est[i], sigma_gp) - cauchy_cdf( boost[i] - 1 , titre_est[i], sigma_gp);
        }

        if (mu <= 0) {
            target += -100000;
        } else {
            target += log(mu) + log(0.33 / prop_study[study[ind[i]]]);
        }
    }

    // Study-speific population-level priors
    boost_peak ~ normal(3, 2);
    wane_s ~ normal(-0.003, 0.002);
    
    // GP priors for titre dependent boosting
    rho_b ~ inv_gamma(5, 5);
    alpha_b ~ normal(0, 1);
    eta_b ~ std_normal(); 
    
    // GP priors for waning dependent boosting
    rho_w ~ inv_gamma(5, 5);
    alpha_w ~ normal(0, 0.002);
    eta_w ~ std_normal(); 

    if(ind_vh > 0) {
        z_vh_cat_int[1] ~ normal(0, 1);
        sigma_vh_int ~ normal(0, 0.1);
        z_vh_cat[1] ~ normal(0, 1);
        sigma_vh ~ normal(0, 1);
    }
    if(ind_s > 0) {
        z_s_cat_int[1] ~ normal(0, 1);
        sigma_s_int ~ normal(0, 0.1);
        z_s_cat[1] ~ normal(0, 1);
        sigma_s ~ normal(0, 1);
    }
    if(ind_g > 0) {
        (z_g_cat_int[1]) ~ normal(0, 1);
        sigma_g_int ~ normal(0, 0.1);
        (z_g_cat[1]) ~ normal(0, 1); 
        sigma_g ~ normal(0, 1);
    }
    if(ind_d > 0) {
        (z_d_cat_int[1]) ~ normal(0, 1);
        sigma_d_int ~ normal(0, 0.1);
        (z_d_cat[1]) ~ normal(0, 1);
        sigma_d ~ normal(0, 1);
    }

    // Individual-level priors
    if (ind_ind > 0) {
        z_boost_ind[1] ~ std_normal();
        sigma_b_ind ~ normal(0, 1);
        z_wane_ind[1] ~ std_normal();
        sigma_w_ind ~ normal(0, 0.002);
    }
    sigma_gp ~ exponential(1);
}
generated quantities {

    vector[n_study] ps_study = rep_vector(0,  n_study);

    vector[N_bt] ps_boost = rep_vector(0,  N_bt);
    vector[N_bt] ps_wane = rep_vector(0,  N_bt);

    vector[n_vac_hist] ps_boost_vh = rep_vector(0, n_vac_hist);
    vector[n_site] ps_boost_s = rep_vector(0, n_site);
    vector[n_sex] ps_boost_g = rep_vector(0, n_sex);
    vector[n_age_cat] ps_boost_d = rep_vector(0, n_age_cat);

    matrix[N_bt, n_vac_hist] ps_boost_vh_t = rep_matrix(0, N_bt, n_vac_hist);
    matrix[N_bt, n_site] ps_boost_s_t = rep_matrix(0, N_bt, n_site);
    matrix[N_bt, n_sex] ps_boost_g_t = rep_matrix(0, N_bt, n_sex);
    matrix[N_bt, n_age_cat] ps_boost_d_t = rep_matrix(0, N_bt, n_age_cat);

    matrix[N_bt, n_vac_hist] ps_boost_vh_inter = rep_matrix(0, N_bt, n_vac_hist);
    matrix[N_bt, n_site] ps_boost_s_inter = rep_matrix(0, N_bt, n_site);
    matrix[N_bt, n_sex] ps_boost_g_inter = rep_matrix(0, N_bt, n_sex);
    matrix[N_bt, n_age_cat] ps_boost_d_inter = rep_matrix(0, N_bt, n_age_cat);

    matrix[N_bt, n_vac_hist] ps_boost_vh_sin = rep_matrix(0, N_bt, n_vac_hist);
    matrix[N_bt, n_site] ps_boost_s_sin = rep_matrix(0, N_bt, n_site);
    matrix[N_bt, n_sex] ps_boost_g_sin = rep_matrix(0, N_bt, n_sex);
    matrix[N_bt, n_age_cat] ps_boost_d_sin = rep_matrix(0, N_bt, n_age_cat);

    real mu_temp_b;
    real mu_temp_w;

    real covar;
    real covar_int;
    real covar_sum;
    real covar_sum_inter;
    real covar_sum_sin;

    for (s in 1:n_study) {
    for (i in 1:N_bt) {
        mu_temp_b = boost_peak[s] + boost_td[i];
        mu_temp_w = wane_s[s] + wane_td[i];
        for (a in 1:n_vac_hist) {
        for (b in 1:n_site) {
        for (c in 1:n_sex) {
        for (d in 1:n_age_cat) {
            covar = 0;
            covar_int = 0;
            covar_sum = 0;
            if (ind_vh > 0) {
                covar += (z_vh_cat[ind_vh, a] * sigma_vh[ind_vh]);
                covar_int += (z_vh_cat_int[ind_vh, a] * sigma_vh_int[ind_vh]);
            }
            if (ind_g > 0) {
                covar += (z_g_cat[ind_g, c] * sigma_g[ind_g]);
                covar_int += (z_g_cat_int[ind_g, c] * sigma_g_int[ind_g]);
            }
            if (ind_d > 0) {
                covar += (z_d_cat[ind_d, d] * sigma_d[ind_d]);
                covar_int += (z_d_cat_int[ind_d, c] * sigma_d_int[ind_d]);
            }
            if (ind_s > 0) {
                covar += (z_s_cat[ind_s, b] * sigma_s[ind_s]);
                covar_int += (z_s_cat_int[ind_s, b] * sigma_s_int[ind_s]);
            }      
            covar_sum = mu_temp_b * (1 + covar_int) + covar;

            ps_study[s] += covar_sum * prop_v[s, a] * prop_s[s, b] * prop_g[s, c] * prop_d[s, d] * prop_t[s, i];
            ps_boost[i] += covar_sum * prop_v[s, a] * prop_s[s, b] * prop_g[s, c] * prop_d[s, d] * prop_study[s];
            ps_wane[i] += (mu_temp_w)  * prop_v[s, a] * prop_s[s, b] * prop_g[s, c] * prop_d[s, d] * prop_study[s];

            ps_boost_vh[a] += covar_sum * prop_s[s, b] * prop_g[s, c] * prop_d[s, d] * prop_t[s, i] * prop_study[s];
            ps_boost_s[b] += covar_sum * prop_v[s, a] * prop_g[s, c] * prop_d[s, d] * prop_t[s, i] * prop_study[s];
            ps_boost_g[c] += covar_sum * prop_v[s, a] * prop_s[s, b] * prop_d[s, d] * prop_t[s, i] * prop_study[s];
            ps_boost_d[d] += covar_sum * prop_v[s, a] * prop_s[s, b] * prop_g[s, c] * prop_t[s, i] * prop_study[s];

            ps_boost_vh_t[i, a] += covar_sum * prop_s[s, b] * prop_g[s, c] * prop_d[s, d] * prop_study[s];
            ps_boost_s_t[i, b] += covar_sum * prop_v[s, a] * prop_g[s, c] * prop_d[s, d] * prop_study[s];
            ps_boost_g_t[i, c] += covar_sum * prop_v[s, a] * prop_s[s, b] * prop_d[s, d] * prop_study[s];
            ps_boost_d_t[i, d] += covar_sum * prop_v[s, a] * prop_s[s, b] * prop_g[s, c] * prop_study[s];

            covar_sum_inter = mu_temp_b * (1 + covar_int);
            covar_sum_sin = mu_temp_b + covar;

            ps_boost_vh_inter[i, a] += covar_sum_inter * prop_s[s, b] * prop_g[s, c] * prop_d[s, d] * prop_study[s];
            ps_boost_s_inter[i, b] += covar_sum_inter * prop_v[s, a] * prop_g[s, c] * prop_d[s, d] * prop_study[s];
            ps_boost_g_inter[i, c] += covar_sum_inter * prop_v[s, a] * prop_s[s, b] * prop_d[s, d] * prop_study[s];
            ps_boost_d_inter[i, d] += covar_sum_inter * prop_v[s, a] * prop_s[s, b] * prop_g[s, c] * prop_study[s];

            ps_boost_vh_sin[i, a] += covar_sum_sin * prop_s[s, b] * prop_g[s, c] * prop_d[s, d] * prop_study[s];
            ps_boost_s_sin[i, b] += covar_sum_sin * prop_v[s, a] * prop_g[s, c] * prop_d[s, d] * prop_study[s];
            ps_boost_g_sin[i, c] += covar_sum_sin * prop_v[s, a] * prop_s[s, b] * prop_d[s, d] * prop_study[s];
            ps_boost_d_sin[i, d] += covar_sum_sin * prop_v[s, a] * prop_s[s, b] * prop_g[s, c] * prop_study[s];
        }
        }
        }
        }
    }   
    }

    vector[N] log_lik;
    for (i in 1:N) {
        if (likelihood_type == 1) {
            log_lik[i] = log(normal_cdf( boost[i] + 1 , titre_est[i], sigma_gp) - normal_cdf( boost[i] - 1 , titre_est[i], sigma_gp)) + log(0.33 / prop_study[study[ind[i]]]);
        } else if (likelihood_type == 2) {
            log_lik[i] = log(cauchy_cdf( boost[i] + 1 , titre_est[i], sigma_gp) - cauchy_cdf( boost[i] - 1 , titre_est[i], sigma_gp)) + log(0.33 / prop_study[study[ind[i]]]);
        }
    }
}
