# Load all data for plottuing
## H3 data
swr30 = function(string, nwrap=30) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr30 = Vectorize(swr30)

swr20 = function(string, nwrap=20) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr20 = Vectorize(swr20)

df_boost_uncert_h1 <- best_fit_h1only %>% as_draws_df %>% spread_draws(boost_ind[i], wane_ind[i])
meta_join_h1 <- data.frame(
    i = 1:h1only_stan$N_ind,
    titre_i = h1only_stan$titre_i,
    vac_hist = h1only_stan$vac_hist
)
df_boost_uncert_meta_h1 <- df_boost_uncert_h1 %>% left_join(meta_join_h1) %>% ungroup


df_boost_uncert_h3 <- best_fit_h3only %>% as_draws_df %>% spread_draws(boost_ind[i], wane_ind[i])
meta_join_h3 <- data.frame(
    i = 1:h3only_stan$N_ind,
    titre_i = h3only_stan$titre_i,
    vac_hist = h3only_stan$vac_hist
)
df_boost_uncert_meta_h3 <- df_boost_uncert_h3 %>% left_join(meta_join_h3) %>% ungroup


df_boost_uncert_h3_cell <- best_fit_h3cell %>% as_draws_df %>% spread_draws(boost_ind[i], wane_ind[i])
meta_join_h3_cell <- data.frame(
    i = 1:h3cell_stan$N_ind,
    titre_i = h3cell_stan$titre_i,
    vac_hist = h3cell_stan$vac_hist
)
df_boost_uncert_meta_h3_cell <- df_boost_uncert_h3_cell %>% left_join(meta_join_h3_cell) %>% ungroup



df_boost_uncert_meta_un <- bind_rows(
    df_boost_uncert_meta_h1 %>% mutate(subtype = "H1N1"),
    df_boost_uncert_meta_h3 %>% mutate(subtype = "H3N2")
)

recode_vhist <- c("1" = "<2 vaccines in last 5 seasons", "2" = "<2 vaccines in last 5 seasons", 
    "3" = "2 or more vaccines in last 5 seasons", "4" = "2 or more vaccines in last 5 seasons",
    "5" = "2 or more vaccines in last 5 seasons", "6" = "2 or more vaccines in last 5 seasons")
df_boost_uncert_meta <- df_boost_uncert_meta_un %>% mutate(vac_hist = recode(vac_hist, !!!recode_vhist))


individual_level_bw <- df_boost_uncert_meta %>%
    summarise(boost_ind = mean(boost_ind), wane_ind = mean(wane_ind), .by = c(titre_i, i, subtype, vac_hist))
individual_level_bw_cell <- df_boost_uncert_meta_h3_cell %>%
    mutate(vac_hist = recode(vac_hist, !!!recode_vhist)) %>% 
    summarise(boost_ind = mean(boost_ind), wane_ind = mean(wane_ind), .by = c(titre_i, i, vac_hist))


full_pp_dist <- 1:nrow(individual_level_bw) %>% 
    map_df(
        ~data.frame(
            ind = as.numeric(individual_level_bw[.x, 2]),
            subtype = as.character(individual_level_bw[.x, 3]),
            titre_i = as.numeric(individual_level_bw[.x, 1]),
            vac_hist = as.character(individual_level_bw[.x, 4]),
            time = 1:365,
            titre = as.numeric(individual_level_bw[.x, 5]) + c(1:365) * as.numeric(individual_level_bw[.x, 6] )
        )
    )


full_pp_dist_cell <- 1:nrow(individual_level_bw_cell) %>% 
    map_df(
        ~data.frame(
            ind = as.numeric(individual_level_bw_cell[.x, 2]),
            titre_i = as.numeric(individual_level_bw_cell[.x, 1]),
            vac_hist = as.character(individual_level_bw_cell[.x, 3]),
            time = 1:365,
            titre = as.numeric(individual_level_bw_cell[.x, 4]) + c(1:365) * as.numeric(individual_level_bw_cell[.x, 5] )
        )
    )


######################################################################
## Individual-level kinetics  ##
######################################################################

full_pp_traj <- full_pp_dist %>% rename(t = titre_i) %>% add_titre_info_trunc %>%
    group_by(subtype, titre_vals, vac_hist, time) %>% ggdist::mean_qi() 
   
full_pp_traj$vac_hist <- swr30(full_pp_traj$vac_hist)



## Fig 1b, comoparison of lines and waning from model
get_marginal_boosts <- function(best_fit_h1only) {
    boosting_t <- best_fit_h1only %>% as_draws_df %>% spread_draws(ps_boost_vh_t[t, v]) %>% add_titre_info 
    waning_t <- best_fit_h1only %>% as_draws_df %>% spread_draws(ps_wane[t]) %>% add_titre_info 

    comparbw_t <- boosting_t %>% left_join(waning_t, by = c("t", ".chain", ".iteration", ".draw", "titre_vals")) %>% add_titre_info %>%
        rename(boost_peak = ps_boost_vh_t, wane_s = ps_wane)
    comparbw_t_mean <- comparbw_t %>% mutate(v = recode(v, !!!recode_vhist) ) %>% group_by(titre_vals, v) %>% summarise(boost_peak = mean(boost_peak), wane_s = mean(wane_s))
}

convert_to_segment <- function(df) {
    df %>% mutate(x1 = 0, x2 = 220, y1 = boost_peak + 14 * wane_s, y2 = boost_peak + 220 * wane_s )
}

find_dur_4fold <- function(df) {
    df %>% mutate(dur_4fold = pmax(0, (boost_peak - 2)/-wane_s)) 
}

comparbw_t_meanh1n1 <- get_marginal_boosts(best_fit_h1only) %>% convert_to_segment %>% find_dur_4fold
comparbw_t_meanh3n2 <- get_marginal_boosts(best_fit_h3only) %>% convert_to_segment %>% find_dur_4fold
comparbw_t_meanh3n2cell <- get_marginal_boosts(best_fit_h3cell) %>% convert_to_segment %>% find_dur_4fold

comparbw_t_mean %>% convert_to_segment %>%
    ggplot() + 
            geom_hline(yintercept = 0, linetype = "dashed") +
            geom_hline(yintercept = 2, linetype = "dashed") +
            geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2, colour = titre_vals), alpha = 1, size = 1) + 
            theme_bw() + 
            facet_grid(cols = vars(v)) +
            scale_y_continuous(breaks = c(-1:5), labels = 2^c(-1:5)) + 
            labs(x = "Days post-vaccination", y = "Model-predicted HAI titre fold-rise", 
                        color = "Pre-vaccination HAI titre") +
            ggtitle("H1N1: Vaccine-induced HAI kinetics to vaccinating strains")



p1 <- full_pp_traj %>% filter(subtype == "H1N1") %>% filter(time <= 220) %>%
    ggplot() + 
        geom_hline(yintercept = 2, color = "gray50", linetype = "dashed") + 
        geom_hline(yintercept = 0, color = "gray50") + 
        geom_line(aes(x = time, y = titre, color = titre_vals), size = 2, alpha = 0.8) + 
        facet_grid(cols = vars(vac_hist)) + theme_bw() + 
        scale_y_continuous(breaks = 0:4, labels = 2^c(0:4), limits = c(-1, 4.5)) +
          labs(x = "Days post-vaccination", y = "Model-predicted HAI titre fold-rise", 
            color = "Pre-vaccination HAI titre") + 
     ggtitle("H1N1: Vaccine-induced HAI kinetics to vaccinating strains")

p2 <- full_pp_traj %>% filter(subtype == "H3N2") %>% filter(time <= 220) %>%
    ggplot() + 
        geom_hline(yintercept = 2, color = "gray50", linetype = "dashed") + 
        geom_hline(yintercept = 0, color = "gray50") + 
        geom_line(aes(x = time, y = titre, color = titre_vals), size = 2, alpha = 0.8) + 
        facet_grid(cols = vars(vac_hist)) + theme_bw() + 
        scale_y_continuous(breaks = 0:4, labels = 2^c(0:4), limits = c(-1, 4.5)) +
          labs(x = "Days post-vaccination", y = "Model-predicted HAI titre", 
            color = "Pre-vaccination HAI titre") + 
         ggtitle("H3N2: Vaccine-induced HAI kinetics to vaccinating strains")


full_pp_traj_cell <- full_pp_dist_cell %>% rename(t = titre_i) %>% add_titre_info_trunc %>%
    group_by(titre_vals, vac_hist, time) %>% ggdist::mean_qi() 
   
full_pp_traj_cell$vac_hist <- swr30(full_pp_traj_cell$vac_hist)

p3 <- full_pp_traj_cell %>% filter(time <= 220) %>%
    ggplot() + 
        geom_hline(yintercept = 2, color = "gray50", linetype = "dashed") + 
        geom_hline(yintercept = 0, color = "gray50") + 
        geom_line(aes(x = time, y = titre, color = titre_vals), size = 2, alpha = 0.8) + 
        theme_bw() + 
        scale_y_continuous(breaks = 0:4, labels = 2^c(0:4), limits = c(-1, 4.5)) +
        facet_grid(cols = vars(vac_hist)) + theme_bw() + 
          labs(x = "Days post-vaccination", y = "Model-predicted HAI titre fold-rise", 
            color = "Pre-vaccination HAI titre") + 
        ggtitle("H3N2: Vaccine-induced HAI kinetics to circulating strains")

recode_subtype <- c("H1N1" = "A(H1N1) vaccinating", "H3N2" = "A(H3N2) vaccinating")

p4 <- bind_rows(full_pp_traj %>%  mutate(subtype = recode(subtype, !!!recode_subtype)), full_pp_traj_cell %>% mutate(subtype = "A(H3N2) circulating")) %>%
    filter(titre <= 2) %>% group_by(subtype, vac_hist, titre_vals) %>% mutate(j = row_number()) %>% filter(j == 1) %>% 
    mutate(subtype = factor(subtype, levels = c("A(H1N1) vaccinating", "A(H3N2) vaccinating", "A(H3N2) circulating"))) %>%
    ggplot() + 
        geom_col(aes(x = titre_vals, y = time, group = vac_hist, fill = vac_hist), width = 0.6, position = position_dodge(0.6)) + 
        facet_grid(cols = vars(subtype)) + theme_bw() + 
        labs(y = "Days post-vaccination until\n<4 fold-rise", x = "Pre-vaccination HAI titre", fill = "Vaccine history") 


#p1 / p2 / p3 / p4  + plot_layout(guides = "collect", heights = c(1, 1, 1, 2)) + plot_annotation(tag_levels = "A")
#ggsave(file = here::here("outputs", "figs", "main", "fig3.pdf"), height = 14, width = 10)

### Values for manuscripts
full_pp_traj %>% filter(subtype == "H1N1", vac_hist == "<2 vaccines in last 5 seasons") %>%
    filter(titre <= 2) %>% group_by(titre_vals) %>% mutate(j = row_number()) %>% filter(j == 1)
full_pp_traj %>% filter(subtype == "H1N1", vac_hist == "2 or more vaccines in last 5\nseasons") %>%
    filter(titre <= 2) %>% group_by(titre_vals) %>% mutate(j = row_number()) %>% filter(j == 1)

full_pp_traj %>% filter(subtype == "H3N2", vac_hist == "<2 vaccines in last 5 seasons") %>%
    filter(titre <= 2) %>% group_by(titre_vals) %>% mutate(j = row_number()) %>% filter(j == 1)
full_pp_traj %>% filter(subtype == "H3N2", vac_hist == "2 or more vaccines in last 5\nseasons") %>%
    filter(titre <= 2) %>% group_by(titre_vals) %>% mutate(j = row_number()) %>% filter(j == 1)

######################################################################
## Individual-level kinetics and their relationtion to protection ##
######################################################################

# A(H1N1) vaccinating strains

df_sens_plot_h1n1 <- full_pp_dist %>% mutate(pos = titre >= 2) %>% filter(subtype == "H1N1") %>% 
    summarise(sens = sum(pos) / n(), pos = sum(pos), n = n(), .by = c(titre_i, time, subtype, vac_hist)) %>% 
    rename(t = titre_i) %>% add_titre_info_trunc %>% 
    mutate(sd =  sqrt(sens * (1 - sens) / n), lb = sens - 2 * sd, , ub = sens + 2 * sd) %>% filter(t <= 7) %>% mutate(metric = "Inidividuals with a four-fold rise in HAI titre")

df_prot_plot_h1n1 <- full_pp_dist %>% filter(subtype == "H1N1") %>% mutate(titre_abs = (titre_i - 1) + titre) %>%
    mutate(prot = titre_abs >= 3) %>% 
    summarise(sens = sum(prot) / n(), prot = sum(prot), n = n(), .by = c(titre_i, time, subtype, vac_hist)) %>% 
    rename(t = titre_i) %>% add_titre_info_trunc %>% 
    mutate(sd =  sqrt(sens * (1 - sens) / n), lb = sens - 2 * sd, , ub = sens + 2 * sd) %>% mutate(metric = "Individuals with an HAI titre greater than 1:40")

df_prot_info_h1n1 <- bind_rows(
    df_sens_plot_h1n1,
    df_prot_plot_h1n1
)

df_prot_plot_h1n1_diff <- df_prot_info_h1n1 %>% select(subtype, time, vac_hist, titre_vals, metric, sens) %>%
    filter(titre_vals %in% c("<10", "10", "20", "40", "80")) %>%
    pivot_wider(names_from = "vac_hist", values_from = "sens") %>% 
    mutate(`Ratio of proportion of highly vaccinated and poorly vaccinated` = `2 or more vaccines in last 5 seasons`/`<2 vaccines in last 5 seasons` ) %>% 
    pivot_longer(`2 or more vaccines in last 5 seasons`:`Ratio of proportion of highly vaccinated and poorly vaccinated`, 
        names_to = "vac_hist", values_to = "sens")

df_prot_plot_h1n1_diff$vac_hist <- swr30(df_prot_plot_h1n1_diff$vac_hist)
df_prot_plot_h1n1_diff$metric <- swr20(df_prot_plot_h1n1_diff$metric)

# A(H3N2) vaccinating strains

df_sens_plot_h3n2 <- full_pp_dist %>% mutate(pos = titre >= 2) %>% filter(subtype == "H3N2") %>% 
    summarise(sens = sum(pos) / n(), pos = sum(pos), n = n(), .by = c(titre_i, time, subtype, vac_hist)) %>% 
    rename(t = titre_i) %>% add_titre_info_trunc %>% 
    mutate(sd =  sqrt(sens * (1 - sens) / n), lb = sens - 2 * sd, , ub = sens + 2 * sd) %>% filter(t <= 7) %>% mutate(metric = "Inidividuals with a four-fold rise in HAI titre")

df_prot_plot_h3n2 <- full_pp_dist %>% filter(subtype == "H3N2") %>% mutate(titre_abs = (titre_i - 1) + titre) %>%
    mutate(prot = titre_abs >= 3) %>% 
    summarise(sens = sum(prot) / n(), prot = sum(prot), n = n(), .by = c(titre_i, time, subtype, vac_hist)) %>% 
    rename(t = titre_i) %>% add_titre_info_trunc %>% 
    mutate(sd =  sqrt(sens * (1 - sens) / n), lb = sens - 2 * sd, , ub = sens + 2 * sd) %>% mutate(metric = "Individuals with an HAI titre greater than 1:40")

df_prot_info_h3n2 <- bind_rows(
    df_sens_plot_h3n2,
    df_prot_plot_h3n2
)

df_prot_plot_h3n2_diff <- df_prot_info_h3n2 %>% select(subtype, time, vac_hist, titre_vals, metric, sens) %>%
    filter(titre_vals %in% c("<10", "10", "20", "40", "80")) %>%
    pivot_wider(names_from = "vac_hist", values_from = "sens") %>% 
    mutate(`Ratio of proportion of highly vaccinated and poorly vaccinated` =  `2 or more vaccines in last 5 seasons`/`<2 vaccines in last 5 seasons` ) %>% 
    pivot_longer(`2 or more vaccines in last 5 seasons`:`Ratio of proportion of highly vaccinated and poorly vaccinated`, 
        names_to = "vac_hist", values_to = "sens")

df_prot_plot_h3n2_diff$vac_hist <- swr30(df_prot_plot_h3n2_diff$vac_hist)
df_prot_plot_h3n2_diff$metric <- swr20(df_prot_plot_h3n2_diff$metric)



# A(H3N2) circulating strains

df_sens_plot_h3n2cell <- full_pp_dist_cell %>% mutate(pos = titre >= 2) %>% 
    summarise(sens = sum(pos) / n(), pos = sum(pos), n = n(), .by = c(titre_i, time, vac_hist)) %>% 
    rename(t = titre_i) %>% add_titre_info_trunc %>% 
    mutate(sd =  sqrt(sens * (1 - sens) / n), lb = sens - 2 * sd, , ub = sens + 2 * sd) %>% filter(t <= 7) %>% mutate(metric = "Inidividuals with a four-fold rise in HAI titre")

df_prot_plot_h3n2cell <- full_pp_dist_cell %>%  mutate(titre_abs = (titre_i - 1) + titre) %>%
    mutate(prot = titre_abs >= 3) %>% 
    summarise(sens = sum(prot) / n(), prot = sum(prot), n = n(), .by = c(titre_i, time, vac_hist)) %>% 
    rename(t = titre_i) %>% add_titre_info_trunc %>% 
    mutate(sd =  sqrt(sens * (1 - sens) / n), lb = sens - 2 * sd, , ub = sens + 2 * sd) %>% mutate(metric = "Individuals with an HAI titre greater than 1:40")

df_prot_info_h3n2cell <- bind_rows(
    df_sens_plot_h3n2cell,
    df_prot_plot_h3n2cell
)

df_prot_plot_h3n2cell_diff <- df_prot_info_h3n2cell %>% select(time, titre_vals, vac_hist, metric, sens) %>%
    filter(titre_vals %in% c("<10", "10", "20", "40", "80")) %>% 
    pivot_wider(names_from = "vac_hist", values_from = "sens") %>% 
    mutate(`Ratio of proportion of highly vaccinated and poorly vaccinated` =  `2 or more vaccines in last 5 seasons`/`<2 vaccines in last 5 seasons` ) %>% 
    pivot_longer(`2 or more vaccines in last 5 seasons`:`Ratio of proportion of highly vaccinated and poorly vaccinated`, 
        names_to = "vac_hist", values_to = "sens")
df_prot_plot_h3n2cell_diff$metric <- swr20(df_prot_plot_h3n2cell_diff$metric)
df_prot_plot_h3n2cell_diff$vac_hist <- swr30(df_prot_plot_h3n2cell_diff$vac_hist)


# Get the AUC

df_prot_plot_h1n1_diff_A <- df_prot_plot_h1n1_diff %>% filter(vac_hist == "<2 vaccines in last 5 seasons")
df_prot_plot_h1n1_diff_B <- df_prot_plot_h1n1_diff %>% filter(vac_hist == "2 or more vaccines in last 5\nseasons")

df_prot_plot_h1n1_auc <- bind_rows(
    df_prot_plot_h1n1_diff_A %>% group_by(titre_vals, vac_hist, metric) %>% summarise(auc = mean(sens)),
    df_prot_plot_h1n1_diff_B %>% group_by(titre_vals, vac_hist, metric) %>% summarise(auc = mean(sens))
) %>% mutate(subtype = "A(H1N1) vaccinating strains")

df_prot_plot_h3n2_diff_A <- df_prot_plot_h3n2_diff %>% filter(vac_hist == "<2 vaccines in last 5 seasons")
df_prot_plot_h3n2_diff_B <- df_prot_plot_h3n2_diff %>% filter(vac_hist == "2 or more vaccines in last 5\nseasons")

df_prot_plot_h3n2_auc <- bind_rows(
    df_prot_plot_h3n2_diff_A %>% group_by(titre_vals, vac_hist, metric) %>% summarise(auc = mean(sens)),
    df_prot_plot_h3n2_diff_B %>% group_by(titre_vals, vac_hist, metric) %>% summarise(auc = mean(sens))
) %>% mutate(subtype = "A(H3N2) vaccinating strains")


df_prot_plot_h3n2cell_diff_A <- df_prot_plot_h3n2cell_diff %>% filter(vac_hist == "<2 vaccines in last 5 seasons")
df_prot_plot_h3n2cell_diff_B <- df_prot_plot_h3n2cell_diff %>% filter(vac_hist == "2 or more vaccines in last 5\nseasons")

df_prot_plot_h3n2cell_auc <- bind_rows(
    df_prot_plot_h3n2cell_diff_A %>% group_by(titre_vals, vac_hist, metric) %>% summarise(auc = mean(sens)),
    df_prot_plot_h3n2cell_diff_B %>% group_by(titre_vals, vac_hist, metric) %>% summarise(auc = mean(sens))
) %>% mutate(subtype = "A(H3N2) circulating strains")

df_prot_plot_auc <- bind_rows(
    df_prot_plot_h1n1_auc,
    df_prot_plot_h3n2_auc,
    df_prot_plot_h3n2cell_auc
)

vac_hist_vals <- c( "<2 vaccines in last 5 seasons", "2 or more vaccines in last 5\nseasons")

p1 <- df_prot_plot_h1n1_diff %>% filter(vac_hist %in% vac_hist_vals) %>%
    ggplot() + 
        geom_point(aes(x = time, y = sens, color = titre_vals ), size = 0.5, alpha = 0.1 ) + 
        geom_smooth(aes(x = time, y = sens, color = titre_vals ), size = 1, linejoin = "round" ) + 
        guides(fill = "none") + 
        theme_bw() + 
        labs(x = "Days post-vaccination", y = "Proportion", 
            color = "Pre-vaccination HAI titre") + 
        ggtitle("H1N1: Vaccine-induced HAI kinetics\n to vaccinating strains") + 
             facet_grid(cols = vars(vac_hist), rows = vars(metric)) 

p2 <- df_prot_plot_h3n2_diff %>% filter(vac_hist %in% vac_hist_vals) %>%
    ggplot() + 
       geom_hline(yintercept = 0, color = "gray") + 
        geom_point(aes(x = time, y = sens, color = titre_vals ), size = 0.5, alpha = 0.1 ) + 
        geom_smooth(aes(x = time, y = sens, color = titre_vals ), size = 1, linejoin = "round" ) + 
        guides(fill = "none") + 
        theme_bw() + 
        labs(x = "Days post-vaccination", y = "Proportion", 
            color = "Pre-vaccination HAI titre") + 
        ggtitle("H3N2: Vaccine-induced HAI kinetics\n to vaccinating strains") + 
             facet_grid(cols = vars(vac_hist), rows = vars(metric))

p3 <- df_prot_plot_h3n2cell_diff %>% filter(vac_hist %in% vac_hist_vals) %>%
    ggplot() + 
       geom_hline(yintercept = 0, color = "gray") + 
        geom_point(aes(x = time, y = sens, color = titre_vals ), size = 0.5, alpha = 0.1 ) + 
        geom_smooth(aes(x = time, y = sens, color = titre_vals ), size = 1, linejoin = "round" ) + 
        guides(fill = "none") + 
        theme_bw() + 
        labs(x = "Days post-vaccination", y = "Proportion", 
            color = "Pre-vaccination HAI titre") + 
        ggtitle("H3N2: Vaccine-induced HAI kinetics\n to circulating strains") + 
            facet_grid(cols = vars(vac_hist), rows = vars(metric))

p4 <- df_prot_plot_auc %>% mutate(subtype = 
        factor(subtype, levels = 
            c("A(H1N1) vaccinating strains", "A(H3N2) vaccinating strains", "A(H3N2) circulating strains")) ) %>% 
    ggplot() + geom_col(aes(x = titre_vals, y = auc, group = vac_hist, fill = vac_hist), width = 0.6, position = position_dodge(0.6)) + 
    facet_grid(vars(metric), vars(subtype)) + theme_bw() + 
    labs(y = "AUC", x = "Pre-vaccination HAI titre", fill = "Vaccine history") + 
     ggtitle("Area under the curve from 0 to 220 days")


(p1 | p2 | p3) / (p4) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
ggsave(file = here::here("outputs", "figs", "main", "fig4.pdf"), height = 13, width = 17)



### Values for manuscripts

df_prot_plot_h1n1_diff %>% filter(time == 220) %>% arrange(titre_vals) %>% as.data.frame


df_prot_plot_h3n2_diff %>% filter(time == 30) %>% arrange(titre_vals) %>% as.data.frame
df_prot_plot_h3n2_diff %>% filter(time == 30) %>% arrange(titre_vals) %>% as.data.frame