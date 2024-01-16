library(flextable)

load(file = here::here("outputs", "data_model", "h3only_hcwonly_df.RData"))
recode_study_table <- c("1" = "2020", "2" = "2021", "3" = "2022") 
h3only_df <- h3only_df %>% mutate(study = recode(study, !!!recode_study_table))

labels <- list(
    variables = list( base_titre = "Pre-vaccine HAI titre",
                    sex = "Sex",
                    age_cat = "Age group (yrs)",
                    site = "Study site",
                    prevac = "Vaccine history \n(# in last 5 years)"),

    groups=list("Study year"))
strata <- c(split(h3only_df, h3only_df$study))
tb1 <- table1::table1(strata, groupspan=c(5), labels = labels, data = h3only_df, caption= "H3N2 vaccine strains" ) %>% as.data.frame
flextable(tb1) %>% 
  save_as_docx(path = here::here("outputs", "tables", "h3only.docx") )

load(file = here::here("outputs", "data_model", "h1only_hcwonly_df.RData"))
h1only_df <- h1only_df %>% mutate(study = recode(study, !!!recode_study_table)) %>%
    mutate(subtype = "H1N1")

labels <- list(
    variables = list( base_titre = "Pre-vaccine HAI titre",
                    sex = "Sex",
                    age_cat = "Age group (yrs)",
                    site = "Study site",
                    prevac = "Vaccine history \n(# in last 5 years)"),

    groups=list("Study year"))
strata <- split(h1only_df, h1only_df$study)
tb1 <- table1::table1(strata, groupspan=c(5), labels = labels, data = h1only_df, caption= "H1N1 vaccine strains" ) %>% as.data.frame
flextable(tb1) %>% 
  save_as_docx(path = here::here("outputs", "tables", "h1only.docx") )


load(file = here::here("outputs", "data_model", "h3cell_hcwonly_df.RData"))
h3cell_df <- h3cell_df %>% mutate(study = recode(study, !!!recode_study_table))

labels <- list(
    variables = list( base_titre = "Pre-vaccine HAI titre",
                    sex = "Sex",
                    age_cat = "Age group (yrs)",
                    site = "Study site",
                    prevac = "Vaccine history \n(# in last 5 years)"),
    groups=list("Study year", "Study year1"))


strata <- c(split(h3cell_df, ~study))
tbA1 <- table1::table1(strata, groupspan=c(5), labels = labels, data = h3cell_df, caption= "H3N2 circulating strains strains" ) %>% as.data.frame
flextable(tbA1) %>% 
  save_as_docx(path = here::here("outputs", "tables", "h3cell.docx") )


hxonly_df <- bind_rows(
    h1only_df %>% mutate(strain_type = "H1N1 vaccine"), 
    h3only_df %>% mutate(strain_type = "H3N2 vaccine"), 
    h3cell_df  %>% mutate(strain_type = "H3N2 circulating") )

label(hxonly_df$base_titre) <- "Pre-vaccine HAI titre"
label(hxonly_df$sex) <- "Sex"
label(hxonly_df$age_cat) <- "Age group (yrs)"
label(hxonly_df$site) <- "Study site"
label(hxonly_df$prevac) <- "Vaccine history \n(# in last 5 years)"
label(hxonly_df$study) <- "Study year"
label(hxonly_df$strain_type) <- "Strain type"

table1::table1(~ base_titre + sex + age_cat + site | strain_type * study, data = hxonly_df, overall = NULL)

table1::table1(~ base_titre + sex + age_cat + site | study * strain_type, labels = labels, data = hxonly_df, caption= "H3N2 circulating strains strains" )

table1::table1(strata, groupspan=c(3), labels = labels, data = hxonly_df, caption= "H3N2 circulating strains strains" )

tbA1 <- table1::table1(strata, groupspan=c(5), labels = labels, data = hxonly_df, caption= "H3N2 circulating strains strains" )

t1flex(tbA1) %>% 
  save_as_docx(path = here::here("outputs", "tables", "h3cell.docx") )