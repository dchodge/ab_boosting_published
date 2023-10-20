get_boost_wane_sample <- function(best_fit_X, subtype) {
    boosting_t <- best_fit_X %>% as_draws_df %>% spread_draws(ps_boost[t]) %>% add_titre_info 
    waning_t <- best_fit_X %>% as_draws_df %>% spread_draws(ps_wane[t]) %>% add_titre_info 

    comparbw_t <- boosting_t %>%
        left_join(waning_t, by = c("t", ".chain", ".iteration", ".draw", "titre_vals")) %>%
        add_titre_info %>%
        mutate(ps_wane = ps_wane * 100) %>% mutate(subtype = subtype)
}


best_fit_vac_compare <- bind_rows(
    get_boost_wane_sample(best_fit_h1only, "A(H1N1) vaccinating"),
    get_boost_wane_sample(best_fit_h3only, "A(H3N2) vaccinating"),
    get_boost_wane_sample(best_fit_h3cell, "A(H3N2) circulating")
) %>% mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H3N2) circulating")))

best_fit_vac_compare_mean <- best_fit_vac_compare %>% group_by(subtype, titre_vals) %>%
    summarise(boost_peak = mean(ps_boost), wane_s = mean(ps_wane))

dfh1vac <- data.frame(t = h1only_stan$titre_i,subtype = "A(H1N1) vaccinating") %>% add_titre_info %>% 
    group_by(subtype, titre_vals) %>% summarise(n = n())
dfh3vac <- data.frame(t = h3only_stan$titre_i,subtype = "A(H3N2) vaccinating") %>% add_titre_info %>% 
    group_by(subtype, titre_vals) %>% summarise(n = n())
dfh3cell <- data.frame(t = h3cell_stan$titre_i,subtype = "A(H3N2) circulating") %>% add_titre_info %>% 
    group_by(subtype, titre_vals) %>% summarise(n = n())

subtype_numbers <- bind_rows(dfh1vac, dfh3vac, dfh3cell) %>% 
    mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H3N2) circulating")))

best_fit_vac_compare_mean <- best_fit_vac_compare_mean %>%
    left_join(subtype_numbers)

########################################################################
####################
########################################################################

p1 <- best_fit_vac_compare %>% ggplot() +
        geom_hline(yintercept = 0, color = "gray80") +
        geom_vline(xintercept = 0, color = "gray80") +
        geom_point(aes(ps_wane, ps_boost, color = titre_vals, shape = subtype), size = 0.1, alpha = 0.1) +
        geom_path(data = best_fit_vac_compare_mean, aes(wane_s, boost_peak, group = subtype), geom = "line") +
        geom_point(data = best_fit_vac_compare_mean, aes(wane_s, boost_peak, fill = titre_vals, size = n, shape = subtype), alpha = 0.8) + theme_bw() +
        scale_shape_manual(values = c(21, 22, 23)) +
        guides(size = "none", color = "none", 
            fill = guide_legend(override.aes = list(shape=21, size = 4)),
            shape = guide_legend(override.aes = list(size = 4)),
        ) +
        scale_y_continuous(breaks = c(-4:7), labels = 2^c(-4:7), limits = c(-1, 5.2)) +
        scale_x_continuous(breaks = seq(-3, 3, 0.1), labels = round(2^seq(-3, 3, 0.1), 2), limits = c(-1.5, 0.3)) +
        theme(axis.text.x = element_text(angle = 90)) +
        labs(y = "Peak fold-rise in HAI titre", x = "Wane in HAI titre after 100 days (fold-rise)", 
            fill = "Pre-vaccination HAI titre", shape = "Strain type", title = "Marginal posterior distribution to strain types for pre-vaccination titre") + 
        facet_grid(cols = vars(subtype))


h1only_plot <- left_join(data.frame(
        ind = h1only_stan$ind,
        boost = h1only_stan$boost
    ),
    data.frame(
        ind = 1:length( h1only_stan$titre_i),
        t = h1only_stan$titre_i,
        v = h1only_stan$vac_hist
    )
) %>% add_titre_info %>% mutate(v = recode(v, !!!coded_labels))

h3only_plot <- left_join(data.frame(
        ind = h3only_stan$ind,
        boost = h3only_stan$boost
    ),
    data.frame(
        ind = 1:length( h3only_stan$titre_i),
        t = h3only_stan$titre_i,
        v = h3only_stan$vac_hist
    )
) %>% add_titre_info %>% mutate(v = recode(v, !!!coded_labels))

h3cell_plot <- left_join(data.frame(
        ind = h3cell_stan$ind,
        boost = h3cell_stan$boost
    ),
    data.frame(
        ind = 1:length( h3cell_stan$titre_i),
        t = h3cell_stan$titre_i,
        v = h3cell_stan$vac_hist
    )
) %>% add_titre_info %>% mutate(v = recode(v, !!!coded_labels))

best_fit_h3only_plot <- best_fit_h3only %>% as_draws_df %>% spread_draws(ps_boost_vh_t[t, v]) %>% add_titre_info %>% mutate(v = recode(v, !!!coded_labels))

h3only_plot_n <- h3only_plot %>% summarise(n = n(), .by = c("titre_vals", "v", "boost"))

p1Di <- best_fit_h3only_plot %>% 
    ggplot() +  
    geom_point(data = h3only_plot_n, aes(titre_vals, boost, color = v, size = n), position = position_dodge(0.7), shape = 0, alpha = 0.5) +
    stat_pointinterval(aes(x = titre_vals, y = ps_boost_vh_t, color = v), size = 2, shape = 21, position = position_dodge(0.7)) +  theme_bw() + 
    guides(size = "none") + 
    labs(x = "Pre-vaccination HAI titre values", y = "Peak fold-rise in HAI titre", color = "Number of \nprevious vaccines", 
        title = "Marginal posterior distribution to A(H3N2) vaccinating strains for pre-vaccination titre\n and vaccine history") + 
    scale_y_continuous(breaks = -1:9, labels = round(2^c(-1:9), 2), limits = c(-1, 9))

best_fit_h1only_plot <- best_fit_h1only %>% as_draws_df %>% spread_draws(ps_boost_vh_t[t, v]) %>% add_titre_info %>% mutate(v = recode(v, !!!coded_labels))

p1Dii <- best_fit_h1only_plot %>% filter(!titre_vals %in% "2560") %>%
    ggplot() +  
    geom_count(data = h1only_plot, aes(titre_vals, boost, color = v), position = position_dodge(0.7), shape = 0, alpha = 0.5) +
    stat_pointinterval(aes(x = titre_vals, y = ps_boost_vh_t, color = v), size = 2, shape = 21, position = position_dodge(0.7)) +  theme_bw() + 
    guides(size = "none") + 
    labs(x = "Pre-vaccination HAI titre values", y = "Peak fold-rise in HAI titre", color = "Number of \nprevious vaccines", 
        title = "Marginal posterior distribution for A(H1N1) vaccinating strains for pre-vaccination titre\n and vaccine history") + 
    scale_y_continuous(breaks = -1:9, labels = round(2^c(-1:9), 2), limits = c(-1, 9))

best_fit_h3cell_plot <- best_fit_h3cell %>% as_draws_df %>% spread_draws(ps_boost_vh_t[t, v]) %>% add_titre_info %>% mutate(v = recode(v, !!!coded_labels))

p1Diii <- best_fit_h3cell_plot %>% filter(!titre_vals %in% "2560") %>%
    ggplot() +  
    geom_count(data = h3cell_plot, aes(titre_vals, boost, color = v), position = position_dodge(0.7), shape = 0, alpha = 0.5) +
    stat_pointinterval(aes(x = titre_vals, y = ps_boost_vh_t, color = v), size = 2, shape = 21, position = position_dodge(0.7)) +  theme_bw() + 
    guides(size = "none") + 
    labs(x = "Pre-vaccination HAI titre values", y = "Peak fold-rise in HAI titre", color = "Number of \nprevious vaccines", 
        title = "Marginal posterior distribution for A(H3N2) circulating strains for pre-vaccination titre\n and vaccine history") + 
    scale_y_continuous(breaks = -1:9, labels = round(2^c(-1:9), 2), limits = c(-1, 9))


p0A <- p1Dii / p1Di / p1Diii + plot_layout(guides = "collect")
p0 <- p1 / p0A + plot_layout(heights = c(1, 3)) + plot_annotation(tag_levels = "A")
ggsave(file = here::here("outputs", "figs", "main", "fig2.pdf"), height = 15, width = 10)


######## ######## ######## ######## ######## 
######## Metrics for manuscript ######## 
######## ######## ######## ######## ######## 


best_fit_vac_compare %>% filter(subtype == "H1N1 vaccine", titre_vals == "<10") %>% pull(ps_boost) %>% 
    quantile(c(0.05, 0.5, 0.975)) %>% 2^.

best_fit_vac_compare %>% filter(subtype == "H1N1 vaccine", titre_vals == "10") %>% pull(ps_boost) %>% 
    quantile(c(0.05, 0.5, 0.975))  %>% 2^.

best_fit_vac_compare %>% filter(subtype == "H1N1 vaccine", titre_vals == "20") %>% pull(ps_boost) %>% 
    quantile(c(0.05, 0.5, 0.975)) %>% 2^.

best_fit_vac_compare %>% filter(subtype == "H1N1 vaccine", titre_vals == "40") %>% pull(ps_boost) %>% 
    quantile(c(0.05, 0.5, 0.975)) %>% 2^.

mp_vh <- best_fit_h1only %>% as_draws_df %>% spread_draws(ps_boost_vh[v]) %>%
    mutate(v = recode(v, !!!coded_labels), covar = "Vaccine history") %>%
    rename(covar_vals = v, boost = ps_boost_vh)

mp_vh %>% filter(covar_vals == 0) %>% pull(boost) %>% quantile(c(0.05, 0.5, 0.975)) %>% 2^.
mp_vh %>% filter(covar_vals %in% 2:5) %>% pull(boost) %>% quantile(c(0.05, 0.5, 0.975)) %>% 2^.


best_fit_vac_compare %>% filter(subtype == "H3N2 vaccine", titre_vals == "<10") %>% pull(ps_boost) %>% 
    quantile(c(0.05, 0.5, 0.975)) %>% 2^.

best_fit_vac_compare %>% filter(subtype == "H3N2 vaccine", titre_vals == "10") %>% pull(ps_boost) %>% 
    quantile(c(0.05, 0.5, 0.975))  %>% 2^.

best_fit_vac_compare %>% filter(subtype == "H3N2 vaccine", titre_vals == "20") %>% pull(ps_boost) %>% 
    quantile(c(0.05, 0.5, 0.975)) %>% 2^.

df_marg_tv <- best_fit_h3only %>% as_draws_df %>% spread_draws(ps_boost_vh_t[t, v]) %>% add_titre_info %>% mutate(v = recode(v, !!!coded_labels)) 

df_marg_tv %>% filter(titre_vals == "<10", v == "0") %>% pull(ps_boost_vh_t) %>% quantile(c(0.05, 0.5, 0.975)) %>% 2^.
df_marg_tv %>% filter(titre_vals == "<10", v == "5") %>% pull(ps_boost_vh_t) %>% quantile(c(0.05, 0.5, 0.975)) %>% 2^.