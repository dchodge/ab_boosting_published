recode_vhist <- c("1" = "<2 vaccines in last 5 seasons", "2" = "<2 vaccines in last 5 seasons", 
    "3" = "2 or more vaccines in last 5 seasons", "4" = "2 or more vaccines in last 5 seasons",
    "5" = "2 or more vaccines in last 5 seasons", "6" = "2 or more vaccines in last 5 seasons")
    
get_marginal_boosts_uncert <- function(best_fit_hXonly) {
    boosting_t <- best_fit_hXonly %>% as_draws_df %>% spread_draws(ps_boost_vh_t[t, v]) %>% add_titre_info_trunc 
    waning_t <- best_fit_hXonly %>% as_draws_df %>% spread_draws(ps_wane[t]) %>% add_titre_info_trunc 

    comparbw_t <- boosting_t %>% left_join(waning_t, by = c("t", ".chain", ".iteration", ".draw", "titre_vals")) %>% add_titre_info_trunc %>%
        rename(boost_peak = ps_boost_vh_t, wane_s = ps_wane)
    comparbw_t %>% mutate(v = recode(v, !!!recode_vhist) ) %>% group_by(t, v, titre_vals, .draw) %>% 
        summarise(boost_peak = mean(boost_peak), wane_s = mean(wane_s)) 
}

## Fig 1b, comoparison of lines and waning from model
get_marginal_boosts <- function(best_fit_hXonly) {
    boosting_t <- best_fit_hXonly %>% as_draws_df %>% spread_draws(ps_boost_vh_t[t, v]) %>% add_titre_info_trunc 
    waning_t <- best_fit_hXonly %>% as_draws_df %>% spread_draws(ps_wane[t]) %>% add_titre_info_trunc 

    comparbw_t <- boosting_t %>% left_join(waning_t, by = c("t", ".chain", ".iteration", ".draw", "titre_vals")) %>% add_titre_info_trunc %>%
        rename(boost_peak = ps_boost_vh_t, wane_s = ps_wane)
    comparbw_t_mean <- comparbw_t %>% mutate(v = recode(v, !!!recode_vhist) ) %>% group_by(titre_vals, v) %>%
         summarise(boost_peak_sd = sd(boost_peak), wane_s_sd = sd(wane_s) )
    comparbw_t_mean_sd <- comparbw_t %>% mutate(v = recode(v, !!!recode_vhist) ) %>% group_by(titre_vals, v) %>%
        summarise(boost_peak = mean(boost_peak), wane_s = mean(wane_s))
    comparbw_t_mean %>% left_join(comparbw_t_mean_sd) %>% 
        mutate(boost_peak_lb = boost_peak - 2 * boost_peak_sd, boost_peak_ub = boost_peak + 2 * boost_peak_sd,
            wane_s_lb = wane_s - 2 * wane_s_sd, wane_s_ub = wane_s + 2 * wane_s_sd
        )
}

convert_to_segment <- function(df) {
    df %>% mutate(x1 = 0, x2 = 220, y1 = boost_peak + 14 * wane_s, y2 = boost_peak + 220 * wane_s )
}

find_dur_4fold <- function(df) {
    df %>% mutate(dur_4fold = pmax(0, (boost_peak - 2)/-wane_s),
        dur_4fold_lb = pmax(0, (boost_peak_lb - 2)/-wane_s), dur_4fold_up = pmax(0, (boost_peak_ub - 2)/-wane_s)) 
}

find_dur_4fold_alt <- function(df) {
    df %>% mutate(dur_4fold = pmax(0, (boost_peak - 2)/-wane_s) )
}

find_dur_above40 <- function(df) {
    df %>% mutate(dur_40 = pmax(0, (boost_peak + (t - 1) - 3)/ -wane_s)) 
}

comparbw_t_meanh1n1_uncert <- get_marginal_boosts_uncert(best_fit_h1only) %>% convert_to_segment %>% find_dur_4fold_alt %>% 
    find_dur_above40 %>%
    mutate(dur_4fold_trunc = pmin(dur_4fold, 365), dur_40_trunc = pmin(dur_40, 365))
comparbw_t_meanh1n1cell_uncert <- get_marginal_boosts_uncert(best_fit_h1cell) %>% convert_to_segment %>% find_dur_4fold_alt %>% 
    find_dur_above40 %>%
    mutate(dur_4fold_trunc = pmin(dur_4fold, 365), dur_40_trunc = pmin(dur_40, 365))
comparbw_t_meanh3n2_uncert <- get_marginal_boosts_uncert(best_fit_h3only) %>% convert_to_segment %>% find_dur_4fold_alt %>% 
    find_dur_above40 %>%
    mutate(dur_4fold_trunc = pmin(dur_4fold, 365), dur_40_trunc = pmin(dur_40, 365))
comparbw_t_meanh3n2cell_uncert <- get_marginal_boosts_uncert(best_fit_h3cell) %>% convert_to_segment %>% find_dur_4fold_alt %>% 
    find_dur_above40 %>%
    mutate(dur_4fold_trunc = pmin(dur_4fold, 365), dur_40_trunc = pmin(dur_40, 365))

comparbw_t_mean_uncert <- bind_rows(
    comparbw_t_meanh1n1_uncert %>% mutate(subtype = "A(H1N1) vaccinating"),
    comparbw_t_meanh1n1cell_uncert %>% mutate(subtype = "A(H1N1) circulating"),
    comparbw_t_meanh3n2_uncert %>% mutate(subtype = "A(H3N2) vaccinating"),
    comparbw_t_meanh3n2cell_uncert %>% mutate(subtype = "A(H3N2) circulating")
)

p1 <- comparbw_t_mean_uncert %>% mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H1N1) circulating", "A(H3N2) circulating"))) %>%
    ggplot() + 
        geom_boxplot(aes(x = titre_vals, y = dur_4fold_trunc, fill = v), alpha = 0.8, outlier.shape = NA) + 
            theme_bw() + 
            facet_grid(cols = vars(subtype)) + 
            scale_y_continuous(
                limits = c(0, 365), breaks = c(0, 50, 100, 150, 200, 250, 300, 365), labels = c("0", "50", "100", "150", "200", "250", "300", ">365")) +
            labs(y = "Days post-vaccination above 4 fold-rise (individual-level)", x = "Pre-vaccination HAI titre", fill = "Vaccine history", color = "Vaccine history") 

p2 <- comparbw_t_mean_uncert %>% mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H1N1) circulating",  "A(H3N2) circulating"))) %>%
    ggplot() + 
        geom_boxplot(aes(x = titre_vals, y = dur_40_trunc, fill = v), alpha = 0.8, outlier.shape = NA) + 
            theme_bw() + 
            facet_grid(cols = vars(subtype)) + 
            scale_y_continuous(
                limits = c(0, 365), breaks = c(0, 50, 100, 150, 200, 250, 300, 365), labels = c("0", "50", "100", "150", "200", "250", "300", ">365")) +
            labs(y = "Days post-vaccination above 4 fold-rise (individual-level)", x = "Pre-vaccination HAI titre", fill = "Vaccine history", color = "Vaccine history") 
p1 / p2


comparbw_t_mean_uncert_marginal_eff <- comparbw_t_mean_uncert %>% select(titre_vals, v, dur_4fold_trunc, subtype, .draw) %>%
    pivot_wider(names_from = "v", values_from = "dur_4fold_trunc") %>% 
    mutate(diff_4fold = `<2 vaccines in last 5 seasons` - `2 or more vaccines in last 5 seasons`) %>% 
    filter(titre_vals != ">160") %>% 
    mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H3N2) circulating"))) 

comparbw_t_mean_uncert_marginal_eff_40 <- comparbw_t_mean_uncert %>% select(titre_vals, v, dur_40_trunc, subtype, .draw) %>%
    pivot_wider(names_from = "v", values_from = "dur_40_trunc") %>% 
    mutate(diff_40 = `<2 vaccines in last 5 seasons` - `2 or more vaccines in last 5 seasons`) %>% 
    filter(titre_vals != ">160") %>% 
    mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H3N2) circulating"))) 
    

p3 <- comparbw_t_mean_uncert_marginal_eff  %>%
     ggplot() + 
        geom_boxplot(aes(x = titre_vals, y = diff_4fold)) + 
        facet_grid(cols = vars(subtype)) + theme_bw() + 
        labs(x = "Pre-vaccination HAI titre", y = "Marginal effect of vaccine history on duration of 4-fold rise")

p4 <- comparbw_t_mean_uncert_marginal_eff_40  %>%
     ggplot() + 
        geom_boxplot(aes(x = titre_vals, y = diff_40)) + 
        facet_grid(cols = vars(subtype)) + theme_bw() + 
        labs(x = "Pre-vaccination HAI titre", y = "Marginal effect of vaccine history on duration of 4-fold rise")

p1 / p2 / p3 / p4






comparbw_t_meanh1n1 <- get_marginal_boosts(best_fit_h1only) %>% convert_to_segment %>% find_dur_4fold
comparbw_t_meanh3n2 <- get_marginal_boosts(best_fit_h3only) %>% convert_to_segment %>% find_dur_4fold
comparbw_t_meanh1n1cell <- get_marginal_boosts(best_fit_h1cell) %>% convert_to_segment %>% find_dur_4fold
comparbw_t_meanh3n2cell <- get_marginal_boosts(best_fit_h3cell) %>% convert_to_segment %>% find_dur_4fold

comparbw_t_meanXX <- bind_rows(
    comparbw_t_meanh1n1 %>% mutate(subtype = "A(H1N1) vaccinating"),
    comparbw_t_meanh3n2 %>% mutate(subtype = "A(H3N2) vaccinating"),
    comparbw_t_meanh1n1cell %>% mutate(subtype = "A(H1N1) circulating"),
    comparbw_t_meanh3n2cell %>% mutate(subtype = "A(H3N2) circulating")
) %>% mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H1N1) circulating", "A(H3N2) circulating")))



comparbw_t_meanXX %>% convert_to_segment %>%
    ggplot() + 
            geom_hline(yintercept = 0, linetype = "dashed") +
            geom_hline(yintercept = 2, linetype = "dashed") +
            geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = titre_vals), alpha = 1, size = 1) + 
            theme_bw() + 
            facet_grid(cols = vars(v), rows = vars(subtype)) +
            scale_y_continuous(breaks = c(-1:5), labels = 2^c(-1:5)) + 
            labs(x = "Days post-vaccination", y = "Model-predicted HAI titre fold-rise", 
                        color = "Pre-vaccination HAI titre") +
            ggtitle("Vaccine-induced HAI kinetics to influenza strains") + 
            theme(strip.text = element_text(size = 12), text = element_text(size = 12)) 

ggsave(file = here::here("outputs", "figs", "main", "fig2.pdf"), height = 12, width = 10)