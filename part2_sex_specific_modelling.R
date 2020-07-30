# IAP modelling

library(rjags)
library(tidyverse)
library(metafor)
load.module("glm")
load.module("dic")

# Data ----
mydfs <- list.files("gender.int/", patt = "rds$", full.names = TRUE)
mydfs <- map(mydfs, readRDS)
names(mydfs) <- list.files("gender.int/", patt = "rds$") %>%  str_sub(1, -5)

mydfs <- bind_rows(mydfs, .id = "data_name")
mydfs <- mydfs %>%
  separate(data_name, into = c("outcome", "qualifier"), sep = "\\.", fill = "right")


mydfs_inter <- mydfs %>%
  filter(qualifier %in% c("men", "women")) %>% 
  arrange(outcome, s, qualifier) %>% 
  select(s, outcome, qualifier, est, se) %>% 
  group_by(outcome, s) %>% 
  dplyr::summarise(est = est[2] - est[1],
            se = sum(se)) %>% 
  ungroup()

cvd_inter <- rma(yi = est,sei = se, data = mydfs_inter[mydfs_inter$outcome=="cvd",] ,method = "ML")
predict(cvd_inter, transf = exp, digit=2)

resp_inter <- rma(yi = est,sei = se, data = mydfs_inter[mydfs_inter$outcome=="resp",] ,method = "ML")
predict(resp_inter, transf = exp, digit=2)





