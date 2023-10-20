plot_effect_sizes <- function(data_stan, stan_fit, filename, study_labels_alt) {


    custom_plot <- scale_y_continuous(breaks = -5:7, labels = 2^c(-5:7), limits = c(0, 5)) 
    custom_plot2 <- scale_y_continuous(breaks = seq(-4, 1, 0.2), labels = round(2^seq(-4, 1, 0.2), 2)) 

    boosting_s <- stan_fit %>% as_draws_df %>% spread_draws(boost_peak[s]) %>% mutate(study = recode(s, !!!study_labels))
    waning_s <- stan_fit %>% as_draws_df %>% spread_draws(wane_s[s]) %>% mutate(study = recode(s, !!!study_labels)) %>% 
        mutate(wane_s = wane_s* 100)

    boosting_t <- stan_fit %>% as_draws_df %>% spread_draws(boost_td[t]) %>% add_titre_info 
    waning_t <- stan_fit %>% as_draws_df %>% spread_draws(wane_td[t]) %>% add_titre_info 

    p1_s <- boosting_s %>% 
        ggplot() + 
            stat_pointinterval(aes(x = study, y = boost_peak)) + 
            labs(x = "Study year", y = "Posterior distribution \nstudy-year-specific boosting")
    p2_s <- waning_s %>% 
        ggplot() + 
            stat_pointinterval(aes(x = study, y = wane_s)) + 
            labs(x = "Study year", y = "Posterior distribution \nstudy-year-specific titre\n waning over 100 days")

    p1 <- p1_s + p2_s & theme_bw()


    p1_t <- boosting_t %>% 
        ggplot() + 
            stat_pointinterval(aes(x = titre_vals, y = boost_td)) +
            labs(x = "Pre-vaccination HAI titre", y = "Posterior distribution of \ntitre boosting")
    p2_t <- waning_t %>% 
        ggplot() + 
            stat_pointinterval(aes(x = titre_vals, y = wane_td)) + 
            labs(x = "Pre-vaccination HAI titre", y = "Posterior distribution of \ntitre waning")

    p2 <- p1_t + p2_t & theme_bw()


        
    ## Get intercepts and interactions comparison
    ### Intercepts information
    inter_vh <- stan_fit %>% as_draws_df %>% spread_draws(z_vh_cat[i, v], sigma_vh[i]) %>%
        mutate(value = z_vh_cat * sigma_vh) %>% mutate(group = recode(v, !!!coded_labels)) %>%
        select(i, group, value) %>% mutate(covariate = "Vaccine history")
    inter_a <- stan_fit %>% as_draws_df %>% spread_draws(z_d_cat[i, v], sigma_d[i]) %>%
        mutate(value = z_d_cat * sigma_d) %>% mutate(group = recode(v, !!!age_labels)) %>%
        select(i, group, value) %>% mutate(covariate = "Age group")
    inter_g <- stan_fit %>% as_draws_df %>% spread_draws(z_g_cat[i, v], sigma_g[i]) %>%
        mutate(value = z_g_cat * sigma_g) %>% mutate(group = recode(v, !!!sex_labels)) %>%
        select(i, group, value) %>% mutate(covariate = "Sex")
    inter_s <- stan_fit %>% as_draws_df %>% spread_draws(z_s_cat[i, v], sigma_s[i]) %>%
        mutate(value = z_s_cat * sigma_s) %>% mutate(group = recode(v, !!!site_labels)) %>%
        select(i, group, value) %>% mutate(covariate = "Site")

    inter_all <- bind_rows(inter_vh, inter_a, inter_g, inter_s) %>% mutate(post_type = "Intercept")
    
    inter_vh_int <- stan_fit %>% as_draws_df %>% spread_draws(z_vh_cat_int[i, v], sigma_vh_int[i]) %>%
        mutate(value = z_vh_cat_int * sigma_vh_int) %>% mutate(group = recode(v, !!!coded_labels)) %>%
        select(i, group, value) %>% mutate(covariate = "Vaccine history")
    inter_a_int <- stan_fit %>% as_draws_df %>% spread_draws(z_d_cat_int[i, v], sigma_d_int[i]) %>%
        mutate(value = z_d_cat_int * sigma_d_int) %>% mutate(group = recode(v, !!!age_labels)) %>%
        select(i, group, value) %>% mutate(covariate = "Age group")
    inter_g_int <- stan_fit %>% as_draws_df %>% spread_draws(z_g_cat_int[i, v], sigma_g_int[i]) %>%
        mutate(value = z_g_cat_int * sigma_g_int) %>% mutate(group = recode(v, !!!sex_labels)) %>%
        select(i, group, value) %>% mutate(covariate = "Sex")
    inter_s_int <- stan_fit %>% as_draws_df %>% spread_draws(z_s_cat_int[i, v], sigma_s_int[i]) %>%
        mutate(value = z_s_cat_int * sigma_s_int) %>% mutate(group = recode(v, !!!site_labels)) %>%
        select(i, group, value) %>% mutate(covariate = "Site")

    inter_all_int <- bind_rows(inter_vh_int, inter_a_int, inter_g_int, inter_s_int) %>% mutate(post_type = "Interaction")

    all_covar_par <- bind_rows(
        inter_all,
        inter_all_int
    )

    p4 <- all_covar_par %>% 
        ggplot() + geom_vline(xintercept = 0, color = "gray40") +
        stat_pointinterval(aes(x = value, y = group), .width = 0.95) +
        facet_grid(vars(covariate), vars(post_type), scales = "free_y", space = "free_y") +
        labs(x = "Posterior distribution of parameter", y = "Covariate") + theme_bw() 

    p1 / (p2) / p4  + plot_layout(widths = c(1, 1, 2), heights = c(1, 1, 2), guide = "collect") + 
        plot_annotation(tag_levels = "A")
    ggsave(here::here("outputs", "figs", "fits_hcwonly", paste0("posteriors_", filename, ".pdf")), height = 12, width = 10)

}

plot_marginal_post <- function(data_stan, stan_fit, filename, study_labels) {

    custom_plot <- scale_y_continuous(breaks = -5:7, labels = 2^c(-5:7), limits = c(0, 5)) 
    custom_plot2 <- scale_y_continuous(breaks = seq(-4, 1, 0.2), labels = round(2^seq(-4, 1, 0.2), 2)) 

    boosting_t <- stan_fit %>% as_draws_df %>% spread_draws(ps_boost[t]) %>% add_titre_info 
    waning_t <- stan_fit %>% as_draws_df %>% spread_draws(ps_wane[t]) %>% add_titre_info 

    boost_year <- stan_fit %>% as_draws_df %>% spread_draws(ps_study[s]) %>% mutate(s = recode(s, !!!study_labels) )

    p1_t <- boosting_t %>% 
        ggplot() + 
            stat_pointinterval(aes(x = titre_vals, y = ps_boost)) +
            scale_y_continuous(breaks = c(0:5), labels = 2^c(0:5)) + 
            theme_bw() + theme(axis.text.x = element_text(angle = 90)) + 
            labs(x = "Pre-vaccine HAI titre", y = "Marginal posterior \ndistribution of titre boosting")
    p2_t <- waning_t %>% 
        ggplot() + 
            stat_pointinterval(aes(x = titre_vals, y = ps_wane)) + 
             theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
            labs(x = "Pre-vaccine HAI titre", y = "Marginal posterior \ndistribution of titre waning")

    p3_t <- boost_year %>% 
        ggplot() + 
            stat_pointinterval(aes(x = s, y = ps_study)) +            
            scale_y_continuous(breaks = c(0:5), labels = 2^c(0:5)) + 
            theme_bw() + theme(axis.text.x = element_text(angle = 90)) + 
            labs(x = "Study year", y = "Marginal posterior \ndistribution of boosting")

    p1 <- p1_t + p2_t + p3_t 

    mp_vh <- stan_fit %>% as_draws_df %>% spread_draws(ps_boost_vh[v]) %>%
        mutate(v = recode(v, !!!coded_labels), covar = "Vaccine history") %>%
        rename(covar_vals = v, boost = ps_boost_vh)
    mp_s <- stan_fit %>% as_draws_df %>% spread_draws(ps_boost_s[s]) %>%
        mutate(s = recode(s, !!!site_labels), covar = "Site") %>%
        rename(covar_vals = s, boost = ps_boost_s)
    mp_g <- stan_fit %>% as_draws_df %>% spread_draws(ps_boost_g[g]) %>%
        mutate(g = recode(g, !!!sex_labels), covar = "Sex") %>%
        rename(covar_vals = g, boost = ps_boost_g)
    mp_a <- stan_fit %>% as_draws_df %>% spread_draws(ps_boost_d[a]) %>%
        mutate(a = recode(a, !!!age_labels), covar = "Age group") %>%
        rename(covar_vals = a, boost = ps_boost_d)

    md_covar <- bind_rows(mp_vh, mp_s, mp_g, mp_a)

    p2 <- md_covar %>% 
        ggplot() + stat_pointinterval(aes(covar_vals, boost)) + 
        scale_y_continuous(breaks = c(0:5), labels = 2^c(0:5)) + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90)) + 
        facet_grid(cols = vars(covar), scales = "free_x", space = "free") +
        labs(x = "", y = "Marginal posterior \ndistribution of boosting")


    mp_vh_int <- stan_fit %>% as_draws_df %>% spread_draws(ps_boost_vh_t[t, v]) %>%
        mutate(v = recode(v, !!!coded_labels), covar = "Vaccine history") %>%
        rename(covar_vals = v, boost = ps_boost_vh_t)
    mp_s_int <- stan_fit %>% as_draws_df %>% spread_draws(ps_boost_s_t[t,s]) %>%
        mutate(s = recode(s, !!!site_labels), covar = "Site") %>%
        rename(covar_vals = s, boost = ps_boost_s_t)
    mp_g_int <- stan_fit %>% as_draws_df %>% spread_draws(ps_boost_g_t[t, g]) %>%
        mutate(g = recode(g, !!!sex_labels), covar = "Sex") %>%
        rename(covar_vals = g, boost = ps_boost_g_t)
    mp_a_int <- stan_fit %>% as_draws_df %>% spread_draws(ps_boost_d_t[t, a]) %>%
        mutate(a = recode(a, !!!age_labels), covar = "Age group") %>%
        rename(covar_vals = a, boost = ps_boost_d_t)

    md_covar_int <- bind_rows(mp_vh_int, mp_s_int, mp_g_int, mp_a_int) %>% add_titre_info

    p3a <- md_covar_int %>% 
        filter(covar == "Vaccine history") %>% 
        ggplot() + stat_pointinterval(aes(titre_vals, boost, color = covar_vals), position = position_dodge(0.7)) + 
        scale_y_continuous(breaks = c(0:5), labels = 2^c(0:5)) + 
        labs(x = "", y = "Marginal posterior \ndistribution of boosting", color = "Vaccine history")
    p3b <- md_covar_int %>% 
        filter(covar == "Site") %>% 
        ggplot() + stat_pointinterval(aes(titre_vals, boost, color = covar_vals), position = position_dodge(0.7)) + 
        scale_y_continuous(breaks = c(0:5), labels = 2^c(0:5)) + 
        labs(x = "", y = "Marginal posterior \ndistribution of boosting", color = "Site")
    p3c <- md_covar_int %>% 
            filter(covar == "Sex") %>% 
            ggplot() + stat_pointinterval(aes(titre_vals, boost, color = covar_vals), position = position_dodge(0.7)) + 
            scale_y_continuous(breaks = c(0:5), labels = 2^c(0:5)) + 
            labs(x = "", y = "Marginal posterior \ndistribution of boosting", color = "Sex")
    p3d <- md_covar_int %>% 
            filter(covar == "Age group") %>% 
            ggplot() + stat_pointinterval(aes(titre_vals, boost, color = covar_vals), position = position_dodge(0.7)) + 
            scale_y_continuous(breaks = c(0:5), labels = 2^c(0:5)) + 
            labs(x = "", y = "Marginal posterior \ndistribution of boosting", color = "Age group")      
        
    p3 <- p3a / p3b / p3c / p3d

    p1 / p2 / p3 + plot_layout(heights = c(1, 1, 3)) + 
        plot_annotation(tag_levels = "A")
    ggsave(here::here("outputs", "figs", "fits_hcwonly", paste0("marginal_", filename, ".pdf")), height = 12, width = 10)

}