
# HCW H3N2 egg
cmdArgs <- commandArgs()
recode_study  <- cmdArgs[1:3]

# HCW H3N2 egg
raw_sero <- read.csv(here::here("data", "hcw", "serology.csv"), sep = ",", header = TRUE, row.names = NULL)
raw_sero <- raw_sero %>%  mutate(titre = log2(titre / 5))
raw_meta <- read.csv (here::here("data", "hcw", "participants.csv"), sep = ",", header = TRUE, row.names = NULL)

meta <- raw_meta %>% select(pid, site, gender, dob, height, weight, recruitment_year, age_screening, bmi)
meta_clean <- meta %>% mutate(dob = ymd(dob), yob = year(dob)) %>% mutate(decade_cat = 
    case_when(
        (yob <= 1969) & (yob >= 1900 ) ~ 1,
        (yob <= 1979) & (yob >= 1970 ) ~ 2,
        (yob <= 1989) & (yob >= 1980 ) ~ 3,
        (yob <= 2010) & (yob >= 1990 ) ~ 4
    )) %>% select(pid, site, gender, age_screening, decade_cat)
bleed_dates_raw <- read.csv(here::here("data", "hcw", "bleed-dates.csv")) %>% mutate(date = ymd(date))
bleed_dates <- bleed_dates_raw %>%
        left_join(bleed_dates_raw %>% filter(day == 0) %>% select(!day) %>% rename(start_date = date), by = c("pid", "year")
    ) %>% 
    mutate(day_true = as.numeric(date - start_date))

sero_vac_strains <- raw_sero %>% left_join(bleed_dates %>% select(pid, year, day, day_true)) %>% left_join(meta_clean) %>% 
    rename(day_cat = day) %>% rename(day = day_true)

raw_vac_full_23 <- read.csv(here::here("data", "hcw", "priorVacc.csv"), sep = ",", header = TRUE, row.names = NULL)

recode_vac <- c("1" = "1", "0" = "0")
vac_full <- raw_vac_full_23 %>% 
    select(c(pid = PID, X2022, X2021, X2020, X2019, X2018, X2017, X2016, X2015))  %>% 
    pivot_longer(!pid, names_to = "vac_year", values_to = "vac_outcome") %>% 
    mutate(vac_year = substr(vac_year, 2, 5)) %>%
    as.data.frame

impute_missing_vacs <- function(vac_hist_XXX) {
    vac_hist_keep <- vac_hist_XXX %>% group_by(pid) %>% summarise(no_missing = sum(is.na(vac_outcome))) %>%
        filter(no_missing <= 1) %>% pull(pid)
    vac_hist_XXX %>% filter(pid %in% vac_hist_keep) %>% group_by(pid) %>% fill(vac_outcome, .direction = "updown")
}

vac_hist_2020 <- vac_full %>% filter(!vac_year %in% c(2020, 2021, 2022)) %>% impute_missing_vacs
vac_hist_2021 <- vac_full %>% filter(!vac_year %in% c(2015, 2021, 2022)) %>% impute_missing_vacs
vac_hist_2022 <- vac_full %>% filter(!vac_year %in% c(2015, 2016, 2022)) %>% impute_missing_vacs

sero_vac_strains_base <- sero_vac_strains %>% filter(day_cat == 0) %>% select(pid, virus, year, titre) %>% rename(base_titre = titre)
sero_vac_strains_reg <- sero_vac_strains %>% filter(day_cat != 0) %>% left_join(sero_vac_strains_base) %>% select(!c(decade_cat, subtype, virus_egg_cell))

make_hcw_year <- function(sero_vac_strains_reg, year_XX, vac_hist_XXX, strain) {
    lol <- sero_vac_strains_reg %>% filter(year == year_XX, virus == strain) %>% 
        mutate(y_type = case_when((day <= short_u & day >= short_l) ~ "short", (day >= long_l & day <= long_u )~ "long") ) %>% 
        filter(!is.na(y_type)) %>% select(pid, day, y_type, titre, site, gender, age = age_screening, base_titre)
    
    hcw_XXX_full <- sero_vac_strains_reg %>% filter(year == year_XX, virus == strain) %>% 
        mutate(y_type = case_when((day <= short_u & day >= short_l) ~ "short", (day >= long_l & day <= long_u )~ "long") ) %>% 
        filter(!is.na(y_type)) %>% select(pid, day, y_type, titre, site, gender, age = age_screening, base_titre) %>% 
        left_join(vac_hist_XXX, by = "pid", multiple = "all") %>% filter(!is.na(vac_year)) %>% 
        arrange(pid)
        
    hcw_XXX_vac <- hcw_XXX_full %>% select(pid, vac_year, vac_outcome) %>% group_by(pid) %>% unique %>% select(!vac_year) %>%
        mutate(row = row_number()) %>% pivot_wider(names_from = row, values_from = vac_outcome) %>% 
        select(pid, `5`, `4`, `3`, `2`, `1`) %>% unite("coded", `5`:`1`, sep = "")
    hcw_XXX_full %>% select(!c(vac_year, vac_outcome, vac_outcome)) %>% left_join(hcw_XXX_vac) %>% unique %>% 
        mutate(pid = paste(pid, year_XX, sep = "_"))

}

make_reg_hcw <- function(data_df) {
    data_df %>% mutate(gender = str_to_title(gender), site = str_to_title(site), titre_change = titre - base_titre, 
        base_titre = as.character(base_titre) ) %>% change_age %>% rename(sex = gender, days = day)
}

subtype <- cmdArgs[4]

hcw_reg_2020 <- make_hcw_year(sero_vac_strains_reg, 2020, vac_hist_2020, cmdArgs[1] ) %>% make_reg_hcw %>% base_titre_convert %>% drop_na
hcw_reg_2021 <- make_hcw_year(sero_vac_strains_reg, 2021, vac_hist_2021, cmdArgs[2] ) %>% make_reg_hcw %>% base_titre_convert %>% drop_na
hcw_reg_2022 <- make_hcw_year(sero_vac_strains_reg, 2022, vac_hist_2022, cmdArgs[3] ) %>% make_reg_hcw %>% base_titre_convert %>% drop_na

save(hcw_reg_2020, file = here::here("outputs", "data_clean", paste0("hcw_reg_2020_", subtype,".RData") ))
save(hcw_reg_2021, file = here::here("outputs", "data_clean", paste0("hcw_reg_2021_", subtype,".RData") ))
save(hcw_reg_2022, file = here::here("outputs", "data_clean", paste0("hcw_reg_2022_", subtype,".RData") ))