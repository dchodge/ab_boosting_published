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
    df %>% mutate(dur_40 = pmax(0, (boost_peak + (t - 1) - 3)/ - wane_s)) 
}

find_dur_above80 <- function(df) {
    df %>% mutate(dur_80 = pmax(0, (boost_peak + (t - 1) - 4)/ - wane_s)) 
}

comparbw_t_meanh1n1_uncert <- get_marginal_boosts_uncert(best_fit_h1only) %>% convert_to_segment %>% find_dur_4fold_alt %>% 
    find_dur_above40 %>% find_dur_above80 %>%
    mutate(dur_4fold_trunc = pmin(dur_4fold, 365), dur_40_trunc = pmin(dur_40, 365), dur_80_trunc = pmin(dur_80, 365))
comparbw_t_meanh3n2_uncert <- get_marginal_boosts_uncert(best_fit_h3only) %>% convert_to_segment %>% find_dur_4fold_alt %>% 
    find_dur_above40 %>% find_dur_above80 %>%
    mutate(dur_4fold_trunc = pmin(dur_4fold, 365), dur_40_trunc = pmin(dur_40, 365), dur_80_trunc = pmin(dur_80, 365))

comparbw_t_meanh1n1cell_uncert <- get_marginal_boosts_uncert(best_fit_h1cell) %>% convert_to_segment %>% find_dur_4fold_alt %>% 
    find_dur_above40 %>% find_dur_above80 %>%
    mutate(dur_4fold_trunc = pmin(dur_4fold, 365), dur_40_trunc = pmin(dur_40, 365), dur_80_trunc = pmin(dur_80, 365))
comparbw_t_meanh3n2cell_uncert <- get_marginal_boosts_uncert(best_fit_h3cell) %>% convert_to_segment %>% find_dur_4fold_alt %>% 
    find_dur_above40 %>% find_dur_above80 %>%
    mutate(dur_4fold_trunc = pmin(dur_4fold, 365), dur_40_trunc = pmin(dur_40, 365), dur_80_trunc = pmin(dur_80, 365))

comparbw_t_mean_uncert <- bind_rows(
    comparbw_t_meanh1n1_uncert %>% mutate(subtype = "A(H1N1) vaccinating"),
    comparbw_t_meanh3n2_uncert %>% mutate(subtype = "A(H3N2) vaccinating"),
    comparbw_t_meanh1n1cell_uncert %>% mutate(subtype = "A(H1N1) circulating"),
    comparbw_t_meanh3n2cell_uncert %>% mutate(subtype = "A(H3N2) circulating")
) %>% mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating",  "A(H1N1) circulating", "A(H3N2) circulating"))) %>%
    filter(titre_vals != ">160") %>% 
    pivot_longer(c(dur_4fold_trunc, dur_40_trunc, dur_80_trunc), names_to = "heuristic", values_to = "values") %>%
    mutate(heuristic = recode(heuristic, dur_4fold_trunc = "4-fold boost", dur_40_trunc = "HAI titre \u2265 1:40" , dur_80_trunc = "HAI titre \u2265 1:80" ))

p1 <- comparbw_t_mean_uncert %>% mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H1N1) circulating", "A(H3N2) circulating"))) %>%
    filter(subtype %in% c("A(H1N1) circulating", "A(H3N2) circulating")) %>% 
    ggplot() + 
        geom_boxplot(aes(x = titre_vals, y = values, fill = v), alpha = 0.8, outlier.shape = NA) + 
            theme_bw() + 
            facet_grid(cols = vars(subtype), rows = vars(heuristic)) + 
            scale_y_continuous(
                limits = c(0, 365), breaks = c(0, 50, 100, 150, 200, 250, 300, 365), labels = c("0", "50", "100", "150", "200", "250", "300", ">365")) +
            labs(y = "PPD of days post-vaccination where heuristic holds", x = "Pre-vaccination HAI titre", fill = "Vaccine history", color = "Vaccine history")  + 
            theme(strip.text = element_text(size = 12), text = element_text(size = 12)) 


comparbw_t_mean_uncert_marginal_eff <- comparbw_t_mean_uncert %>% select(titre_vals, v, heuristic, values, subtype, .draw) %>%
    pivot_wider(names_from = "v", values_from = "values") %>% 
    mutate(diff_heuristic = `<2 vaccines in last 5 seasons` - `2 or more vaccines in last 5 seasons`) %>% 
    filter(titre_vals != ">160") %>% 
    mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H1N1) circulating", "A(H3N2) circulating"))) 


p2 <- comparbw_t_mean_uncert_marginal_eff %>% filter(subtype %in% c("A(H1N1) circulating", "A(H3N2) circulating"))  %>%
     ggplot() + 
        geom_boxplot(aes(x = titre_vals, y = diff_heuristic)) + 
        facet_grid(cols = vars(subtype), rows = vars(heuristic)) + theme_bw() + 
        labs(x = "Pre-vaccination HAI titre", y = "Marginal effect of vaccination history on heuristic in days\n(Infrequently vs. frequently vaccinated)")  + 
    theme(strip.text = element_text(size = 12), text = element_text(size = 12)) 


p1 / p2  + plot_annotation(tag_levels = "A")
ggsave(here::here("outputs", "figs", "main", "fig3.png"), height = 14, width = 12)


comparbw_t_mean_uncert %>% filter(subtype == "A(H1N1) vaccinating") %>% 
    group_by(v, titre_vals) %>% ggdist::mean_qi(dur_4fold, .width = 0.95)


comparbw_t_mean_uncert %>% filter(subtype == "A(H3N2) vaccinating") %>% 
    group_by(v, titre_vals) %>% ggdist::mean_qi(dur_4fold, .width = 0.95)