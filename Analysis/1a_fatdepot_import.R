#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(tictoc)
library(furrr)

setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Fatdepots/")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/"), ]
options(future.globals.maxSize= 5e9)
plan(multicore, workers = 20, gc = TRUE)#plan(sequential)

############What you need to change ###########
parameters <- default_param()
parameters$ldref <- "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs" #only maf > 0.01

ID_mrbase_exp <- c(NULL)
exp_mrbase <- NULL#paste0("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/", ID_mrbase_exp, "/", ID_mrbase_exp, ".vcf.gz")
ID_server_exp <- c(paste0("dis-13-", 2:4), "dis-2-1", df_index[grepl("dis-6-",id), id], "dis-14-2", "dis-14-7",
                   df_index[grepl("trait-16", id) & grepl("^HDL|^LDL|^logTG", trait) & population %in% c("European", "Mixed"), id],
                   "trait-14-8", "trait-14-10", df_index[grepl("trait-2-", id), id],
                   "dis-5-1", "dis-7-1", "trait-12-2", "trait-13-1", "trait-13-2",
                   df_index[grepl("trait-25", id) & grepl("adjbmi3", trait), id],
                   paste0("trait-27-", 1:2), paste0("trait-28-", 1:2)) %>% unique(.)
exp_server <- paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", ID_server_exp, "/", ID_server_exp, ".vcf.gz")
arguments <- data.table(id = c(ID_mrbase_exp,ID_server_exp), path = c(exp_mrbase,exp_server), pval = 5e-8, r2 =0.01, kb = 1000)
#################################################

inst <- future_map(split(arguments, 1:nrow(arguments)), function(x) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  GagnonMR::get_inst(x$path, r2 = x$r2, pval = x$pval, kb = x$kb, should_write = FALSE,  parameters = parameters)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)
setDT(inst)
out <- future_map(as.list(c(exp_mrbase,exp_server)), function(x, rsiid = unique(inst$SNP)) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  outr<-GagnonMR::extract_outcome_variant(snps = rsiid, outcomes = x, parameters = parameters)
  return(outr)
}, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)

instmvmr <- GagnonMR::convert_outcome_to_exposure(out) %>% as.data.table(.)

sdtrait <- tribble(
  ~trait, ~sex, ~sd.exposure,
  "VATadj",   "Males and Females", 1.18,
  "ASATadj",   "Males and Females", 1.23,
  "GFATadj",   "Males and Females", 1.44,
  "VATadj",   "Males", 1.41,
  "ASATadj",   "Males", 1.21,
  "GFATadj",   "Males", 1.36,
  "VATadj",   "Females", 0.91,
  "ASATadj",   "Females", 1.25,
  "GFATadj",   "Females", 1.52,
)
setDT(sdtrait)
sdtrait <- merge(sdtrait, df_index[grepl("trait-25", id), .(id,clean_variable_name,sex)], by.x = c("trait", "sex"),
      by.y = c("clean_variable_name", "sex"))

list_res<-list(instmvmr=instmvmr,inst=inst)
list_res<- map(list_res, function(x) {
 x <- merge(x, df_index[,.(id, sex)], by.x = "id.exposure", by.y = "id")
 x <- merge(x, sdtrait[,.(id, sd.exposure)], by.x = "id.exposure", by.y = "id", all.x = TRUE)
 x[grepl("trait-25-", id.exposure) & sex %in% c("Males", "Females"), beta.exposure := beta.exposure*sd.exposure]
 x[grepl("trait-25-", id.exposure) & sex %in% c("Males", "Females"), se.exposure := se.exposure*sd.exposure]
 x[,sd.exposure := NULL]
 x[,sex:=NULL]
 return(x)
})

fwrite((list_res$instmvmr), "Data/Modified/all_inst_mvmr.txt")
fwrite(GagnonMR::convert_exposure_to_outcome(list_res$instmvmr), "Data/Modified/all_outcome_mvmr.txt" )
fwrite(list_res$inst, "Data/Modified/inst_all_sign_clump.txt")

message("This script finished without errors")
