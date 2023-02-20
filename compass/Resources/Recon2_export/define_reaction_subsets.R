rm(list = ls())
library(tidyverse)
'%nin%' <- function(x,y)!('%in%'(x,y))
printf <- function(...) cat(sprintf(...))

rxn_md = read_csv("Recon2_export/rxn_md.csv")

central_carbon = rxn_md %>%
  # turns out that some of the interesting glycolysis reactions (e.g., HEX1 = hexokinase) have confidence = 2
  # filter(rxn_confidence %in% c("0", "4")) %>%
  filter(str_detect(subsystem, "(Glycolysis|Citric acid cycle)")) %>%
  arrange(rxn_code_nodirection) %>%
  #manually exclude some uninteresting reactions
  filter(!str_detect(rxn_code_nodirection, "^AL(C|D)D2")) %>%
  filter(rxn_code_nodirection %nin% c("ETOHMO")) %>%
  filter(!grepl("\\[(l|r|x|n|g|a)\\]", rxn_formula))

#unique just in case
writeLines(central_carbon$rxn_code_nodirection %>% unique(),
           con = "Recon2_export/central_carbon_rxns.txt")


exchange_rxns = rxn_md %>%
  #no need to include the DM_ reactions that have Exchange/demand as their subsystem,
  #David didn't define them as exchange, it seems
  filter(grepl("^EX_", rxn_code_nodirection)) #| (subsystem == "Exchange/demand reaction"))
writeLines(exchange_rxns$rxn_code_nodirection %>% unique(),
           con = "Recon2_export/exchange_rxns.txt")
