#Load all the necessary packages via pacman
if(!require(pacman)){install.packages("pacman")}

pacman::p_load(
    # utils + data
    "here", "tidyverse", "tidybayes", "socialmixr", "lubridate",
    # modelling
    "rstan", "brms", "cmdstanr", "loo", "modelr", "truncnorm", 
    "table1",
    # plotting
    "ggh4x", "ggplot2", "ggmcmc", "ggthemes", "ggridges", "patchwork", 'viridis',
    # parallel stuff
    "foreach", "doParallel", "furrr", "rstatix"
)

# Load important things
source(here::here("R/utils_func.R"))
source(here::here("R/global_var.R"))