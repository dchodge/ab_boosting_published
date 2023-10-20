short_l <- 8
short_u <- 150
long_l <- 150
long_u <- 360

coded_labels <- as.character(c(0:5))
coded_relabel <- c(1:6)
names(coded_labels) <- coded_relabel


titres_labels <- c("<10", "10", "20", "40", "80", "160", "320", "640", "1280", "2560")
titre_relabel <- c(1:10)
names(titres_labels) <- titre_relabel

titres_labels_trim <- c("<10", "10", "20", "40", "80", "160", ">160", ">160", ">160", ">160")
titre_relabel_trim <- c(1:10)
names(titres_labels_trim) <- titre_relabel_trim


sex_labels <- c("Female", "Male")
sex_relabel <- c(1:2)
names(sex_labels) <- sex_relabel

site_labels <- c("Adelaide", "Brisbane", "Perth", "Newcastle", "Melbourne", "Sydney")
site_relabel <- c(1:6)
names(site_labels) <- site_relabel


age_labels <- c("<30", "30-39", "40-49", "50-59", "60+") 
age_relabel <- c(1:5)
names(age_labels) <- age_relabel

study_labels <- c("HCW2020", "HCW2021", "HCW2022")
study_relabel <- c(1:length(study_labels))
names(study_labels) <- study_relabel


stan_code_gp_prior <-
'functions {
    vector gp_prior(int N_bt, vector[] x, real alpha, real rho,  vector eta) {
        vector[N_bt] gp_values;
        matrix[N_bt, N_bt] K_b_gq;
        matrix[N_bt, N_bt] L_K_b_gq;
        K_b_gq = gp_exp_quad_cov(x, alpha, rho) + diag_matrix(rep_vector(1e-10, N_bt));
        L_K_b_gq = cholesky_decompose(K_b_gq);
        return L_K_b_gq * eta;
   }
  }'
expose_stan_functions(stanc(model_code = stan_code_gp_prior))