########################################################
########     MODEL A, HISTORY H3N2 VACCINE STRAIN     ########
########################################################

load(file = here::here("outputs", "data_clean", "hcw_reg_2020_H3N2_vac.RData") ) # hcw_reg_2020
load(file = here::here("outputs", "data_clean", "hcw_reg_2021_H3N2_vac.RData") ) # hcw_reg_2021
load(file = here::here("outputs", "data_clean", "hcw_reg_2022_H3N2_vac.RData") ) # hcw_reg_2022

h3only_df <- bind_rows(
            hcw_reg_2020 %>% mutate(study = 1) %>% mutate(vac_year = "2020"),
            hcw_reg_2021 %>% mutate(study = 2) %>% mutate(vac_year = "2021"),
            hcw_reg_2022 %>% mutate(study = 3) %>% mutate(vac_year = "2022")
            ) %>% unique %>% 
            mutate(prevac = str_count(coded, "1")) %>% data_addfactors %>% select_and_omit

save(h3only_df, file = here::here("outputs", "data_model", "h3only_hcwonly_df.RData"))
h3only_stan <- make_data_stan(h3only_df, 1:3)
save(h3only_stan, file = here::here("outputs", "data_model", "h3only_hcwonly_stan.RData"))



# Stuff for manuscript (methods)
h3only_df %>% nrow
h3only_df %>% mutate(pid_only = substr(pid, 1, 7) ) %>% pull(pid_only) %>% unique %>% length
h3only_df %>% group_by(pid) %>% mutate(rn = row_number() ) %>% filter(rn == 1) %>% pull(days) %>% quantile(c(0.025, 0.5, 0.975) )
h3only_df %>% group_by(pid) %>% mutate(rn = row_number() ) %>% filter(rn == 2) %>% pull(days) %>% quantile(c(0.025, 0.5, 0.975) )
h3only_df %>% pull(prevac) %>% table

h3only_df %>% mutate(pid_only = substr(pid, 1, 7) ) %>% select(!c(pid, days)) %>% unique %>% pull(prevac) %>% table

########################################################
########     MODEL B, HISTORY H3N2 CELL STRAIN     ########
########################################################

hcw_reg_2020_h3n2_cell_get <- load(file = here::here("outputs", "data_clean", "hcw_reg_2020_H3N2_cell.RData") ) # hcw_reg_2020
hcw_reg_2020_h3n2_cell <- get(hcw_reg_2020_h3n2_cell_get)
hcw_reg_2021_h3n2_cell_get <- load(file = here::here("outputs", "data_clean", "hcw_reg_2021_H3N2_cell.RData") ) # hcw_reg_2021
hcw_reg_2021_h3n2_cell <- get(hcw_reg_2021_h3n2_cell_get)
hcw_reg_2022_h3n2_cell_get <- load(file = here::here("outputs", "data_clean", "hcw_reg_2022_H3N2_cell.RData") ) # hcw_reg_2021
hcw_reg_2022_h3n2_cell <- get(hcw_reg_2022_h3n2_cell_get)

h3cell_df <- bind_rows(
            hcw_reg_2020_h3n2_cell %>% mutate(study = 1, subtype = "h3", vac_year = "2020"),
            hcw_reg_2021_h3n2_cell %>% mutate(study = 2, subtype = "h3", vac_year = "2021"),
            hcw_reg_2022_h3n2_cell %>% mutate(study = 3, subtype = "h3", vac_year = "2022")
            ) %>% unique %>% 
            mutate(prevac = str_count(coded, "1")) %>% data_addfactors %>% select_and_omit

save(h3cell_df, file = here::here("outputs", "data_model", "h3cell_hcwonly_df.RData"))
h3cell_stan <- make_data_stan(h3cell_df, 1:3)
save(h3cell_stan, file = here::here("outputs", "data_model", "h3cell_hcwonly_stan.RData"))

########################################################
########     MODEL C, HISTORY H1N1 VACCINE STRAIN     ########
########################################################

hcw_reg_2020_h1n1_get <- load(file = here::here("outputs", "data_clean", "hcw_reg_2020_H1N1_vac.RData") ) # hcw_reg_2020
hcw_reg_2020_h1n1 <- get(hcw_reg_2020_h1n1_get)
hcw_reg_2021_h1n1_get <- load(file = here::here("outputs", "data_clean", "hcw_reg_2021_H1N1_vac.RData") ) # hcw_reg_2021
hcw_reg_2021_h1n1 <- get(hcw_reg_2021_h1n1_get)
hcw_reg_2022_h1n1_get <- load(file = here::here("outputs", "data_clean", "hcw_reg_2022_H1N1_vac.RData") ) # hcw_reg_2021
hcw_reg_2022_h1n1 <- get(hcw_reg_2022_h1n1_get)

h1only_df <- bind_rows(
            hcw_reg_2020_h1n1 %>% mutate(study = 1, subtype = "h1", vac_year = "2020"),
            hcw_reg_2021_h1n1 %>% mutate(study = 2, subtype = "h1", vac_year = "2021"),
            hcw_reg_2022_h1n1 %>% mutate(study = 3, subtype = "h1", vac_year = "2022")
            ) %>% unique %>% 
            mutate(prevac = str_count(coded, "1")) %>% data_addfactors %>% select_and_omit
save(h1only_df, file = here::here("outputs", "data_model", "h1only_hcwonly_df.RData"))
h1only_stan <- make_data_stan(h1only_df, 1:3)
save(h1only_stan, file = here::here("outputs", "data_model", "h1only_hcwonly_stan.RData"))


########################################################
########     MODEL D, HISTORY H1N1 CELL STRAIN     ########
########################################################

hcw_reg_2020_h1n1_cell_get <- load(file = here::here("outputs", "data_clean", "hcw_reg_2020_H1N1_cell.RData") ) # hcw_reg_2020
hcw_reg_2020_h1n1_cell <- get(hcw_reg_2020_h1n1_cell_get)
hcw_reg_2021_h1n1_cell_get <- load(file = here::here("outputs", "data_clean", "hcw_reg_2021_H1N1_cell.RData") ) # hcw_reg_2021
hcw_reg_2021_h1n1_cell <- get(hcw_reg_2021_h1n1_cell_get)
hcw_reg_2022_h1n1_cell_get <- load(file = here::here("outputs", "data_clean", "hcw_reg_2022_H1N1_cell.RData") ) # hcw_reg_2021
hcw_reg_2022_h1n1_cell <- get(hcw_reg_2022_h1n1_cell_get)

h1cell_df <- bind_rows(
            hcw_reg_2020_h1n1_cell %>% mutate(study = 1, subtype = "h1", vac_year = "2020"),
            hcw_reg_2021_h1n1_cell %>% mutate(study = 2, subtype = "h1", vac_year = "2021"),
            hcw_reg_2022_h1n1_cell %>% mutate(study = 3, subtype = "h1", vac_year = "2022")
            ) %>% unique %>% 
            mutate(prevac = str_count(coded, "1")) %>% data_addfactors %>% select_and_omit

save(h1cell_df, file = here::here("outputs", "data_model", "h1cell_hcwonly_df.RData"))
h1cell_stan <- make_data_stan(h1cell_df, 1:3)
save(h1cell_stan, file = here::here("outputs", "data_model", "h1cell_hcwonly_stan.RData"))


# Plot some data
load(file = here::here("outputs", "data_model", "h1cell_hcwonly_stan.RData"))


left_join(
    data.frame(
        pid = h3cell_stan$ind,
        boost = h3cell_stan$boost
    ),
    data.frame(
        pid = 1:h3cell_stan$N_ind,
        titre_i = h3cell_stan$titre_i
    )
) %>% 
    ggplot() + 
        geom_count(aes(titre_i, boost))


h3only_df %>% 
    ggplot() + 
        geom_count(aes(base_titre, titre_change))



h3only_df %>% mutate(study = recode(study, !!!study_labels)) %>%
    ggplot() + 
        geom_line(aes(days, titre_change, group = pid), alpha = 0.2) + 
        geom_point(aes(days, titre_change, color = base_titre), alpha = 0.7) + 
        geom_smooth(aes(days, titre_change), method = "lm", size = 2, color = "darkred") +
        facet_wrap(vars(study)) + theme_bw() + 
        scale_color_viridis(discrete = TRUE, option = "D") +
        scale_y_continuous(breaks = seq(-3, 9, 3), labels = 2^seq(-3, 9, 3)) + 
        labs(x = "Days post-vaccinaton", y = "HAI fold-rise", color = "Pre-vac HAI titre")
ggsave(file = here::here("outputs", "figs",  "dataplots", "modelA_hcwonly_sum.pdf"), width = 13)

h3cell_df %>% mutate(study = recode(study, !!!study_labels)) %>%
    ggplot() + 
        geom_count(aes(base_titre, titre_change), alpha = 0.7) + 
        stat_summary(aes(base_titre, titre_change, group = 1), size = 0.5, geom = "line", color = "darkred") +
        stat_summary(aes(base_titre, titre_change), size = 0.5, color = "darkred") +
        facet_wrap(vars(study)) + theme_bw() + 
        scale_color_viridis() +
        scale_y_continuous(breaks = seq(-3, 9, 3), labels = 2^seq(-3, 9, 3)) + 
        labs(x = "Pre-vaccination HAI titre", y = "HAI fold-rise", color = "Pre-vac HAI titre")

#ggsave(file = here::here("outputs", "figs",  "dataplots", "modelB_hcwonly_sum.pdf"), height = 5, width = 8)