recode_vhist <- c("1" = "<2 vaccines in last 5 seasons", "2" = "<2 vaccines in last 5 seasons", 
    "3" = "2 or more vaccines in last 5 seasons", "4" = "2 or more vaccines in last 5 seasons",
    "5" = "2 or more vaccines in last 5 seasons", "6" = "2 or more vaccines in last 5 seasons")

get_ind_boosts <- function(best_fit_hXonly, hXonly_stan) {
    boostingwane_t <- best_fit_hXonly %>% as_draws_df %>% spread_draws(boost_ind[i], wane_ind[i])

    meta_join_h1 <- data.frame(
        i = 1:hXonly_stan$N_ind,
        titre_i = hXonly_stan$titre_i,
        vac_hist = hXonly_stan$vac_hist
    )
    df_boost_uncert_meta_h1 <- boostingwane_t %>% left_join(meta_join_h1) %>% ungroup %>%
        summarise(boost_ind = mean(boost_ind), wane_ind = mean(wane_ind), .by = c(titre_i, i, vac_hist)) %>%
        mutate(vac_hist = recode(vac_hist, !!!recode_vhist)) %>% rename(t = titre_i) %>% add_titre_info_trunc
}

ind_traj_h1h1 <- get_ind_boosts(best_fit_h1only, h1only_stan)
ind_traj_h3h2 <- get_ind_boosts(best_fit_h3only, h3only_stan)
ind_traj_h3h2cell <- get_ind_boosts(best_fit_h3cell, h3cell_stan)

convert_to_segment <- function(df) {
    df %>% mutate(x1 = 0, x2 = 220, y1 = boost_ind + 0 * wane_ind, y2 = boost_ind + 220 * wane_ind )
}

find_dur_4fold <- function(df) {
    df %>% mutate(dur_4fold = pmax(0, (boost_ind - 2)/ -wane_ind)) 
}

find_dur_above40 <- function(df) {
    df %>% mutate(dur_40 = pmax(0, (boost_ind + (t - 1) - 3)/ -wane_ind)) 
}

ind_traj_h1h1_ind <- ind_traj_h1h1 %>% convert_to_segment %>% find_dur_4fold %>% find_dur_above40
ind_traj_h3h2_ind <- ind_traj_h3h2 %>% convert_to_segment %>% find_dur_4fold %>% find_dur_above40
ind_traj_h3h2cell_ind <- ind_traj_h3h2cell %>% convert_to_segment %>% find_dur_4fold %>% find_dur_above40

ind_traj_subtype <- bind_rows(
    ind_traj_h1h1_ind %>% mutate(subtype = "A(H1N1) vaccinating"), 
    ind_traj_h3h2_ind %>% mutate(subtype = "A(H3N2) vaccinating"), 
    ind_traj_h3h2cell_ind %>% mutate(subtype = "A(H3N2) circulating")) %>%
    mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H3N2) circulating"))) 


function_add_test_4fold <- function(ind_traj_hXhX_ind, name_string) {
    stat.test <- 
        left_join(
            ind_traj_hXhX_ind %>%
                group_by(titre_vals) %>%
                cohens_d(dur_4fold ~ vac_hist),
            ind_traj_hXhX_ind %>%
            group_by(titre_vals) %>%
            t_test(dur_4fold ~ vac_hist) %>%
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj") 
        ) %>% mutate(subtype = name_string) %>% mutate(heuristic = "Four-fold rise")
    stat.test
}

function_add_test_40 <- function(ind_traj_hXhX_ind, name_string) {
    stat.test <- 
        left_join(
            ind_traj_hXhX_ind %>%
                group_by(titre_vals) %>%
                cohens_d(dur_40 ~ vac_hist),
            ind_traj_hXhX_ind %>%
            group_by(titre_vals) %>%
            t_test(dur_40 ~ vac_hist) %>%
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj") 
        ) %>% mutate(subtype = name_string) %>% mutate(heuristic = "1:40 titre")
    stat.test
}

cal_ks_4fold <- function(ind_traj_h3h2_ind, string_name) { 
        (ind_traj_h3h2_ind$titre_vals %>% levels ) %>% map_df(
        function(x) {
            dist_freq <- ind_traj_h3h2_ind %>% filter(titre_vals == x, vac_hist == "2 or more vaccines in last 5 seasons") %>% pull(dur_4fold)
            dist_infreq <- ind_traj_h3h2_ind %>% filter(titre_vals == x, vac_hist == "<2 vaccines in last 5 seasons") %>% pull(dur_4fold)
            ks_test_out <- ks.test(
                dist_freq,
                dist_infreq
            )
            data.frame(
                subtype = string_name,
                titre_val = x,
                statistic = ks_test_out$statistic,
                diff_mean = (mean(dist_infreq) - mean(dist_freq)) / sd(c(dist_freq, dist_infreq)),
                pval = ks_test_out$p.value
            )
        }
    )
}

cal_ks_40 <- function(ind_traj_h3h2_ind, string_name) { 
        (ind_traj_h3h2_ind$titre_vals %>% levels ) %>% map_df(
        function(x) {
            dist_freq <- ind_traj_h3h2_ind %>% filter(titre_vals == x, vac_hist == "2 or more vaccines in last 5 seasons") %>% pull(dur_40)
            dist_infreq <- ind_traj_h3h2_ind %>% filter(titre_vals == x, vac_hist == "<2 vaccines in last 5 seasons") %>% pull(dur_40)
            ks_test_out <- ks.test(
                dist_freq,
                dist_infreq
            )
            data.frame(
                subtype = string_name,
                titre_val = x,
                statistic = ks_test_out$statistic,
                diff_mean = (mean(dist_infreq) - mean(dist_freq)) / sd(c(dist_freq, dist_infreq)),
                pval = ks_test_out$p.value
            )
        }
    )
}

df_h1n1_ks <- cal_ks_4fold(ind_traj_h1h1_ind, "A(H1N1) vaccinating")
df_h3n2_ks <- cal_ks_4fold(ind_traj_h3h2_ind, "A(H3N2) vaccinating")
df_h3n2cell_ks <- cal_ks_4fold(ind_traj_h3h2cell_ind, "A(H3N2) circulating")

df_ks_4fold <- bind_rows(df_h1n1_ks, df_h3n2_ks, df_h3n2cell_ks) %>% mutate(heuristic = "Four-fold rise")

df_h1n1_ks <- cal_ks_40(ind_traj_h1h1_ind, "A(H1N1) vaccinating")
df_h3n2_ks <- cal_ks_40(ind_traj_h3h2_ind, "A(H3N2) vaccinating")
df_h3n2cell_ks <- cal_ks_40(ind_traj_h3h2cell_ind, "A(H3N2) circulating")

df_ks_40 <- bind_rows(df_h1n1_ks, df_h3n2_ks, df_h3n2cell_ks) %>% mutate(heuristic = "1:40 HAI titre")

df_ks <- bind_rows(df_ks_4fold, df_ks_40) %>% 
    mutate(p_sig = case_when(pval < 0.05~"<0.05", pval >= 0.05~">=0.05")) %>%
    mutate(titre_val = factor(titre_val, levels = c("<10", "10", "20", "40", "80", "160", ">160")))

p3 <- df_ks %>%
    mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H3N2) circulating"))) %>%
    filter(titre_val != ">160") %>% 
    ggplot() + 
    geom_col(aes(x = titre_val, y = diff_mean, fill = p_sig), width = 0.5 ) + 
    facet_grid(cols = vars(subtype), rows = vars(heuristic)) + theme_bw() + 
    scale_fill_manual(values = c("red", "gray")) +
    labs(x = "Pre-vaccination HAI titre", y = 
        "Effect size between infrequently\nvaccinated and frequently vaccinated", fill = "p-value from KS-test") + 
    ggtitle("Measure of effect size of vaccination history for a given heuristic")


stat.test_4fold <- bind_rows(
    function_add_test_4fold(ind_traj_h1h1_ind, "A(H1N1) vaccinating"),
    function_add_test_4fold(ind_traj_h3h2_ind, "A(H3N2) vaccinating"),
    function_add_test_4fold(ind_traj_h3h2cell_ind, "A(H3N2) circulating")
) %>% mutate(p.adj = pmax(round(p.adj, 3), 0.001))

stat.test_40 <- bind_rows(
    function_add_test_40(ind_traj_h1h1_ind, "A(H1N1) vaccinating"),
    function_add_test_40(ind_traj_h3h2_ind, "A(H3N2) vaccinating"),
    function_add_test_40(ind_traj_h3h2cell_ind, "A(H3N2) circulating")
) %>% mutate(p.adj = pmax(round(p.adj, 3), 0.001))


#p3 <- bind_rows(stat.test_4fold , stat.test_40) %>% 
 #   mutate(p.adj.signif = recode(p.adj.signif, 
 ##       "ns" = ">=0.05",
   #     "**" = "<0.05",
 ## #      "*" = "<0.05",
  #      "***" = "<0.05",
 #       )) %>%
 #   mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H3N2) circulating"))) %>%
 #   ggplot() + 
 #   geom_hline(yintercept = 0.0, color = "gray10") +
 #   geom_hline(yintercept = 0.2, color = "gray") +
 #   geom_hline(yintercept = 0.5, color = "gray") +
 #   geom_hline(yintercept = 0.8, color = "gray") +
 #   geom_col(aes(x = titre_vals, y = effsize, fill = p.adj.signif), width = 0.5 ) + 
 #   facet_grid(cols = vars(subtype), rows = vars(heuristic)) + theme_bw() + 
 #   scale_fill_manual(values = c("red", "gray")) +
 #   labs(x = "Pre-vaccination HAI titre", y = "Cohen's d (effect size)", fill = "p-value") + 
 #   ggtitle("Measure of effect size of vaccination history for a given heuristic")


p1 <- ind_traj_subtype %>% mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H3N2) circulating"))) %>%
    filter(titre_vals != ">160") %>% 
    ggplot() + 
        geom_boxplot(aes(x = titre_vals, middle = mean(dur_4fold), y = dur_4fold, fill = vac_hist), alpha = 0.2, outlier.shape = NA) + 
        geom_point(aes(x = titre_vals, y = dur_4fold, color = vac_hist), alpha = 0.2, 
            position = position_jitterdodge(jitter.width = 0.2)) + theme_bw() + 
        #geom_text(data = stat.test_4fold, aes(x = titre_vals, label = 
         #   paste0(p.adj, " (", p.adj.signif, ")")), angle = 90, y = 850) +
            facet_grid(cols = vars(subtype)) + 
            scale_y_continuous(
                limits = c(0, 365), breaks = c(0, 50, 100, 150, 200, 250, 300, 365), labels = c("0", "50", "100", "150", "200", "250", "300", ">365")) +
            labs(y = "Days post-vaccination above 4 fold-rise\n (individual-level)", x = "Pre-vaccination HAI titre", fill = "Vaccine history", color = "Vaccine history") 


p2 <- ind_traj_subtype %>% 
    mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H3N2) circulating"))) %>%
    mutate(dur_40 = pmin(dur_40, 365)) %>%
    filter(titre_vals != ">160") %>% 
    ggplot() + 
        geom_boxplot(aes(x = titre_vals, middle = mean(dur_40) , y = dur_40, fill = vac_hist), alpha = 0.2, outlier.shape = NA) + 
        geom_point(aes(x = titre_vals, y = dur_40, color = vac_hist), alpha = 0.2, 
            position = position_jitterdodge(jitter.width = 0.2)) + theme_bw() + 
       # geom_text(data = stat.test_40, aes(x = titre_vals, label = 
       #     paste0(p.adj, " (", p.adj.signif, ")")), angle = 90, y = 850) +
            facet_grid(cols = vars(subtype)) +
            scale_y_continuous(
                limits = c(0, 365), breaks = c(0, 50, 100, 150, 200, 250, 300, 365), labels = c("0", "50", "100", "150", "200", "250", "300", ">365")) +
            labs(y = "Days post-vaccination with titre value above 1:40\n (individual-level)", x = "Pre-vaccination HAI titre", fill = "Vaccine history", color = "Vaccine history") 

p2 / p1 / p3
ggsave(file = here::here("outputs", "figs", "main", "fig4.pdf"), height = 13, width = 13)