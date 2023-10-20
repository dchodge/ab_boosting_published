change_age <- function(data_df) {
    data_df %>% mutate(
        age_cat = case_when(age < 30 ~ "<30", (age >= 30 & age < 40) ~ "30-39", 
            (age >= 40 & age < 50) ~ "40-49",  (age >= 50 & age < 60) ~ "50-59", age >= 60 ~ "60+"))
}

add_titre_info <- function(df_data) {
    df_data %>% mutate(titre_vals = recode(t, !!!titres_labels)) %>%
     mutate(titre_vals = factor(titre_vals, levels = titres_labels)) 
} 

add_titre_info_trunc <- function(df_data) {
    df_data %>% mutate(titre_vals = recode(t, !!!titres_labels_trim)) %>%
     mutate(titre_vals = factor(titre_vals, levels = unique(titres_labels_trim)) )
} 

base_titre_convert <- function(data_df) {
    data_df %>% mutate(base_titre = as.character(5 * 2 ^ as.numeric(base_titre)))
}

add_age_factor <- function(df) {
    df %>% mutate(age_cat = recode(a, !!!age_cat_levels)) %>% 
    mutate(age_cat = factor(age_cat, levels = age_cat_levels))
}

add_all_factor <- function(df) {
    df %>% 
        mutate(base_titre = factor(base_titre, levels = titres_labels)) %>%
        mutate(site = factor(site, levels = site_levels)) %>% 
        mutate(sex = factor(sex, levels = sex_levels)) %>% 
        mutate(age_cat = factor(age_cat, levels = age_cat_levels)) %>% 
        mutate(prevac = factor(prevac, levels = prevac_levels)) 
}

data_addfactors <- function(df_2clean) {
    sex_key <- df_2clean$sex %>% unique %>% sort
    age_key <- df_2clean$age_cat %>% unique %>% sort 
    site_key <- df_2clean$site %>% unique %>% sort 
    prevac_key <- as.character(0:5)
    df_2clean <- df_2clean %>% mutate(base_titre = factor(base_titre, 
        levels = c("5", "10", "20", "40", "80", "160", "320", "640", "1280", "2560")) ) %>%
        mutate(prevac = factor(prevac, levels = prevac_key )) %>%
        mutate(sex = factor(sex, levels = sex_key )) %>%
        mutate(age_cat = factor(age_cat, levels = age_key )) %>%
        mutate(site = factor(site, levels = site_key)) 
    df_2clean
}

select_and_omit <- function(df_2clean) {
    df_2clean %>% dplyr::select(pid, days, vac_year, base_titre, site, study, sex, age_cat, prevac, titre_change) %>% na.omit
}


make_data_stan <- function(hcw_reg, S) {
    
    hcw_reg <- hcw_reg %>% arrange(pid)
    sero_boost_all_convert <- compose_data(hcw_reg)
    N <- sero_boost_all_convert$n
    ids <- sero_boost_all_convert$pid
    day <- sero_boost_all_convert$days
    titre_change <- sero_boost_all_convert$titre_change

    hcw_reg_meta <- hcw_reg %>% select(pid, vac_year, base_titre, site, study, sex, age_cat, prevac) %>% unique %>% arrange(pid)

    sero_boost_all_convert_meta <- compose_data(hcw_reg_meta)

    titre_base <- sero_boost_all_convert_meta$base_titre
    n_base_titre <- 10

    site <- sero_boost_all_convert_meta$site
    n_site <- sero_boost_all_convert_meta$n_site
    sex <- sero_boost_all_convert_meta$sex
    n_sex <- sero_boost_all_convert_meta$n_sex
    age_cat <- sero_boost_all_convert_meta$age_cat
    n_age_cat <- sero_boost_all_convert_meta$n_age_cat
    vac_hist <- sero_boost_all_convert_meta$prevac
    n_vac_hist <- sero_boost_all_convert_meta$n_prevac
    study <- sero_boost_all_convert_meta$study
    n_study <- length(S)

    get_props_study <- function(sero_boost_all_convert, var,  N, n_study) {
        prop <- matrix(, nrow = n_study, ncol = N)
        for (i in 1:n_study) {
            if (length(sero_boost_all_convert[[var]][sero_boost_all_convert$study == S[i]]) == 0) {
                prop[i, ] <- rep(0, N)
            } else {
                prop[i, ] <- sero_boost_all_convert[[var]][sero_boost_all_convert$study == S[i]] %>%
                        factor(levels = 1:N) %>% table %>% as.numeric %>% `/`(., sum(.))
            }
        }
        prop
    }

    get_props <- function(sero_boost_all_convert, var, N) {
        prop <- vector(mode = "numeric", length = N)
        prop <- sero_boost_all_convert[[var]] %>%
                factor(levels = 1:N) %>% table %>% as.numeric %>% `/`(., sum(.))
        prop
    }

    prop_s <- get_props_study(sero_boost_all_convert, "site", n_site, n_study)
    prop_g <- get_props_study(sero_boost_all_convert, "sex", n_sex, n_study)
    prop_v <- get_props_study(sero_boost_all_convert, "prevac", n_vac_hist, n_study)
    prop_d <- get_props_study(sero_boost_all_convert, "age_cat", n_age_cat, n_study)
    prop_t <- get_props_study(sero_boost_all_convert, "base_titre", n_base_titre, n_study)
    prop_study <- study %>% factor(levels = 1:n_study) %>% table %>% as.numeric %>% `/`(. , sum(.))


    prop_s_all <- get_props(sero_boost_all_convert, "site", n_site)
    prop_g_all <- get_props(sero_boost_all_convert, "sex", n_sex)
    prop_v_all <- get_props(sero_boost_all_convert, "prevac", n_vac_hist)
    prop_d_all <- get_props(sero_boost_all_convert, "age_cat", n_age_cat)
    prop_t_all <- get_props(sero_boost_all_convert, "base_titre", n_base_titre)

    data_list_i <- list(
        N = N,
        boost = titre_change,
        day = day,
        titre_i = titre_base,
        N_bt = n_base_titre,
        bt_vals = 1:n_base_titre,
        N_ind = ids %>% unique %>% length,
        ind = ids,
        n_vac_hist = n_vac_hist,
        vac_hist = vac_hist,

        site = site,
        n_site = n_site,
        sex = sex,
        n_sex = n_sex,
        age_cat = age_cat,
        n_age_cat = n_age_cat,
        n_study = n_study,
        study = study,

        prop_s_all = prop_s_all,
        prop_g_all = prop_g_all,
        prop_v_all = prop_v_all,
        prop_d_all = prop_d_all,
        prop_t_all = prop_t_all,

        prop_s = prop_s,
        prop_g = prop_g,
        prop_v = prop_v,
        prop_d = prop_d,
        prop_t = prop_t,
        prop_study = prop_study,

        # Indicators
        ind_g = 1,
        ind_d = 1,
        ind_s = 1,
        ind_vh = 1,
        ind_ind = 1,
        likelihood_type = 1,

        site_levels = levels(hcw_reg$site),
        age_cat_levels = levels(hcw_reg$age_cat),
        prevac_levels = levels(hcw_reg$prevac),
        sex_levels = levels(hcw_reg$sex)

    )
    data_list_i
}

## All combinations of covariates for regression model
fit_whole_dataset <- function(data_list_i, stanfile, filename) {
    boost_model <- cmdstan_model(here::here("src", "stan", paste0("ab_model", ".stan")), compile = TRUE)
    boost_model_c <- boost_model$compile(stanc_options = list("O1", "canonicalize"))
    fit_model <- boost_model_c$sample(
        data = data_list_i[1:35],    # named list of data
        iter_warmup = 1000,          # number of warmup iterations per chain
        iter_sampling = 1000,            # total number of iterations per chain
        chains = 4,            # number of cores (could use one per chain)
        parallel_chains = 4,
        adapt_delta = 0.8,
        max_treedepth = 15,
        init = 1,
        refresh = 100
    )
    fit_model$save_object(file = here::here("outputs", stanfile, paste0("fit_", filename, ".RData" ) ))
    #file_list <- list.files(path = here::here("outputs", stanfile, "temp" ), pattern = filename)
    #file.remove(paste0(here::here("outputs", stanfile, "temp", sep = ""), file_list))
    fit_model$save_output_files(dir = here::here("outputs", stanfile, "temp" ), basename = filename)
}

posteriors_4_predict <- function(best_fit_h3multi, data_test) {

    posterior_list <- list()

    posterior_list$boost_peak <- best_fit_h3multi$draws("boost_peak", format = "matrix") 
    posterior_list$wane_s <- best_fit_h3multi$draws("wane_s", format = "matrix") 

    posterior_list$boost_td <- best_fit_h3multi$draws("boost_td", format = "matrix") 
    posterior_list$wane_td <- best_fit_h3multi$draws("wane_td", format = "matrix") 

    posterior_list$alpha_b <- best_fit_h3multi$draws("alpha_b") %>% as.numeric
    posterior_list$rho_b <- best_fit_h3multi$draws("rho_b") %>% as.numeric
    posterior_list$eta_b <- best_fit_h3multi$draws("eta_b", format = "matrix") 
    posterior_list$alpha_w <- best_fit_h3multi$draws("alpha_w") %>% as.numeric
    posterior_list$rho_w <- best_fit_h3multi$draws("rho_w") %>% as.numeric
    posterior_list$eta_w <- best_fit_h3multi$draws("eta_w", format = "matrix") 

    if(data_test$ind_vh > 0) {
        posterior_list$z_vh_cat <- best_fit_h3multi$draws("z_vh_cat", format = "matrix") 
        posterior_list$sigma_vh <- best_fit_h3multi$draws("sigma_vh") %>% as.numeric
        posterior_list$z_vh_cat_int <- best_fit_h3multi$draws("z_vh_cat_int", format = "matrix") 
        posterior_list$sigma_vh_int <- best_fit_h3multi$draws("sigma_vh_int") %>% as.numeric
    }
    if(data_test$ind_s > 0) {
        posterior_list$z_s_cat <- best_fit_h3multi$draws("z_s_cat", format = "matrix") 
        posterior_list$sigma_s <- best_fit_h3multi$draws("sigma_s") %>% as.numeric
        posterior_list$z_s_cat_int <- best_fit_h3multi$draws("z_s_cat_int", format = "matrix") 
        posterior_list$sigma_s_int <- best_fit_h3multi$draws("sigma_s_int") %>% as.numeric
    }
    if(data_test$ind_g > 0) {
        posterior_list$z_g_cat <- best_fit_h3multi$draws("z_g_cat", format = "matrix") 
        posterior_list$sigma_g <- best_fit_h3multi$draws("sigma_g") %>% as.numeric 
        posterior_list$z_g_cat_int <- best_fit_h3multi$draws("z_g_cat_int", format = "matrix") 
        posterior_list$sigma_g_int <- best_fit_h3multi$draws("sigma_g_int") %>% as.numeric
    }
    if(data_test$ind_d > 0) {
        posterior_list$z_d_cat <- best_fit_h3multi$draws("z_d_cat", format = "matrix") 
        posterior_list$sigma_d <- best_fit_h3multi$draws("sigma_d") %>% as.numeric
        posterior_list$z_d_cat_int <- best_fit_h3multi$draws("z_d_cat_int", format = "matrix") 
        posterior_list$sigma_d_int <- best_fit_h3multi$draws("sigma_d_int") %>% as.numeric
    }

    posterior_list$sigma_gp <- best_fit_h3multi$draws("sigma_gp") %>% as.numeric
    if(data_test$ind_ind > 0) {
        posterior_list$z_boost_ind <- best_fit_h3multi$draws("z_boost_ind", format = "matrix") 
        posterior_list$sigma_b_ind <- best_fit_h3multi$draws("sigma_b_ind") %>% as.numeric
        posterior_list$z_wane_ind <- best_fit_h3multi$draws("z_wane_ind", format = "matrix") 
        posterior_list$sigma_w_ind <- best_fit_h3multi$draws("sigma_w_ind") %>% as.numeric
    }
    posterior_list$n_samples <- 4000
    posterior_list$n_samples_run <- 100
    posterior_list
}
