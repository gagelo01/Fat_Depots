#!/usr/bin/env Rscript
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(GagnonMR)
library(writexl)
library(furrr)
library(openxlsx)

setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Fatdepots/")
all_inst_mvmr <- fread( "Data/Modified/all_inst_mvmr.txt")
all_outcome_mvmr <- fread( "Data/Modified/all_outcome_mvmr.txt" )
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
res_univariable <- fread("Data/Modified/res_univariate.txt")
zsexdiff <- fread( "Data/Modified/zsexdiff.txt")
harm_univariate <- fread( "Data/Modified/harm_univariate.txt")
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")
FandQ <- fread("Data/Modified/FandQ_univariate.txt")
FandQ <- FandQ[method != "MR Egger", ]
egger_intercept<- fread("Data/Modified/egger_intercept.txt")
resmvmr<-fread( "Data/Modified/res_mvmr.txt")

##
todump<-c("mr_keep.exposure","pval_origin.exposure")
inst_all_sign_clump[,(todump):=NULL]
inst_all_sign_clump <- merge(inst_all_sign_clump, df_index[,.(id,clean_variable_name)], by.x = "id.exposure", by.y = "id")
setnames(inst_all_sign_clump, "clean_variable_name", "clean_variable_name.exposure")


list_table<-list(egger_intercept = egger_intercept, FandQ = FandQ, 
                  res_univariable = res_univariable)

idtoremove<-c("dis-5-1", "dis-14-7", "trait-14-10", res_univariable[(grepl("trait-16", id.outcome) & population.outcome == "European"), unique(id.outcome) ])
todump<-c("outcome", "exposure", "population.exposure", "population.outcome")


list_table <- map(list_table, function(x) {
  
  if(!is.null(x$id.exposure)) {
    x<- x[!(id.exposure %in% idtoremove), ]
  }
  if(!is.null(x$id.outcome)) {
    x<- x[!(id.outcome %in% idtoremove), ]
  }
  x[,forwardMR:=ifelse(grepl("trait-25-", id.exposure), "exposure", "outcome")]
  x[,(todump):=NULL]
  x<- x[order(forwardMR,clean_variable_name.exposure, clean_variable_name.outcome, id.exposure, id.outcome),]
  x[,forwardMR:=NULL]
  setcolorder(x, c("clean_variable_name.exposure", "clean_variable_name.outcome", "id.exposure", "id.outcome"))
  return(x)
})

res_univariable[,type_of_test := NULL]
#########
# 1
dataset <- df_index[id %in% res_univariable$id.exposure, ]
colnom<-c("group_name", "nsnp", "initial_build", "category", "sd", "note", "unit");dataset[,(colnom):=NULL]

dt_title <- data.table(title = paste0("Supplementary Table ", 1:6),
                       caption = c("Description of the datasets used.",
                                   "Instrument used and relevant statistics for each exposures.",
                                   "Instrument strength and heterogeneity statistics for univariable MR",
                                   "Univariable Mendelian Randomization results",
                                   "Univariable Mendelian Randomization Egger's intercept",
                                   "Sex specific Mendelian randomization results and z test"))
list_supdat <- list("Supplementary Data 1" = dataset,
                    "Supplementary Data 2" = inst_all_sign_clump,
                    "Supplementary Data 3" = FandQ,
                    "Supplementary Data 4" = res_univariable, 
                    "Supplementary Data 5" = egger_intercept, 
                    "Supplementary Data 6" = zsexdiff)
                       
                       
#########
col_description<- vector(mode = "list", length = length(list_supdat))
col_description[[1]] <- tribble(
  ~x, ~y,  
  "id", "a unique id", 
  "trait", "a unique name",
  "year", "the year the GWAS was published", 
  "trait", "The author of the GWAS",
  "consortium", "the name of the consortium of the GWAS",
  "sex", "sex included",
  "population", "ancestry",
  "sample_size", "the sample_size",
  "pmid", "the pubmed ID",
  "ncase", "the number of cases",
  "ncontrol", "the number of controls",
  "adjustments", "what variables were included in the GWAS adjustment model",
  "clean_variable_name", "A publication ready name that can be used to plot figures."
) %>% as.data.table(.)

col_description[[2]] <- tribble(
  ~x, ~y,  
  "id.exposure", "a unique id", 
  "chr.exposure", "chromosome",
  "pos.exposure", "Position build 37", 
  "other_allele.exposure", "The non-effect allele",
  "effect_allele.exposure", "The effect allele",
  "beta.exposure", "beta effect estimate",
  "se.exposure", "standard error of the estimate",
  "pval.exposure", "p-valueor of the estimate",
  "eaf.exposure", "effect allele frequency",
  "samplesize.exposure", "sample size",
  "ncase.exposure", "number of cases",
  "SNP", "rsid",
  "ncontrol.exposure", "number of controls",
  "exposure", "A unique name for the exposure",
  "clean_variable_name.exposure", "A publication ready name that can be used to plot figures."
) %>% as.data.table(.)

col_description[[3]] <- tribble(
  ~x, ~y,  
  "id.outcome", "a unique id for the outcome", 
  "id.exposure", "a unique id for the exposure",
  "outcome", "a unique name for the outcome", 
  "exposure", "a unique name for the exposure",
  "method", "The method used for the cochran's Q",
  "Q_df", "cochran's Q degree of freedom",
  "Q_pval", "cochran's Q pvalue",
  "fstat", "F-statistics",
  "rsq", "R squared variance explained",
  "N", "number of SNPs",
  "sample_size.exposure", "sample size of exposure",
  "sample_size.outcome", "sample size of outcome",
  "biasprop", "proportion of bias of the OLS observational estimate that is due to sample overlap",
  "sex.exposure", "the sex of the exposure",
  "population.exposure", "ancestry",
  "sex.outcome", "the sex of the outcome",
  "population.outcome", "ancestry",
  "clean_variable_name.exposure", "A publication ready name that can be used to plot figures.",
  "clean_variable_name.outcome", "A publication ready name that can be used to plot figures."
) %>% as.data.table(.)

col_description[[4]] <- tribble(
  ~x, ~y,  
  "id.outcome", "a unique id for the outcome", 
  "id.exposure", "a unique id for the exposure",
  "outcome", "a unique name for the outcome", 
  "exposure", "a unique name for twosample MR",
  "method", "The method used for the cochran's Q",
  "b", "beta",
  "se", "standard error",
  "pval", "p-value",
  "lci", "lower confidence interval",
  "uci", "upper confidence interval",
  "sex.exposure", "the sex of the exposure",
  "population.exposure", "ancestry",
  "sex.outcome", "the sex of the outcome",
  "population.outcome", "ancestry",
  "clean_variable_name.exposure", "A publication ready name that can be used to plot figures.",
  "clean_variable_name.outcome", "A publication ready name that can be used to plot figures."
) %>% as.data.table(.)

col_description[[5]] <- tribble(
  ~x, ~y,  
  "id.outcome", "a unique id for the outcome", 
  "id.exposure", "a unique id for the exposure",
  "outcome", "a unique name for the outcome", 
  "exposure", "a unique name for twosample MR",
  "egger_intercept", "Egger's intercept",
  "se", "standard error",
  "pval", "p-value",
  "sex.exposure", "the sex of the exposure",
  "population.exposure", "ancestry",
  "sex.outcome", "the sex of the outcome",
  "population.outcome", "ancestry",
  "clean_variable_name.exposure", "A publication ready name that can be used to plot figures.",
  "clean_variable_name.outcome", "A publication ready name that can be used to plot figures."
) %>% as.data.table(.)

col_description[[6]] <- tribble(
  ~x, ~y,
  "clean_variable_name.exposure", "A publication ready name that can be used to plot figures.",
  "clean_variable_name.outcome", "A publication ready name that can be used to plot figures.",
  "method", "the two sample MR method",
  "population.exposure", "ancestry",
  "population.outcome", "ancestry",
  "b", "beta",
  "se", "standard error",
  "pval", "p-value",
  "sexdiff_z_statlci", "the Z score for sex difference",
  "sexdiff_z_stat_pval", "the p-value for sex difference Z score test"
) %>% as.data.table(.)

bold_st <- createStyle(textDecoration = "Bold")
wb <- createWorkbook()
for(i in 1:length(list_supdat)) {
addWorksheet(wb, sheetName =  dt_title[i, title])

title <- dt_title[i,paste0(title, " : ", caption)]
writeData(wb, sheet = i, x = title, startCol = 1, startRow = 1)
addStyle(wb, sheet = i, style = bold_st, rows = 1, cols = 1:2)
writeData(wb, sheet = i, x = col_description[[i]], startCol = 1, startRow = 2)
addStyle(wb, sheet = i, style = bold_st, rows = 2:col_description[[i]][,.N+2], cols = 1)
deleteData(wb, sheet = i, rows = 2, cols = 1:2, gridExpand = TRUE)
writeData(wb, sheet = i, x = list_supdat[[i]], startCol = 1, startRow = col_description[[i]][,.N+4])
addStyle(wb, sheet = i, style = bold_st, rows = col_description[[i]][,.N+4], cols = 1:length(colnames(list_supdat[[i]])), gridExpand = TRUE, stack = TRUE)
}
saveWorkbook(wb, file = "Results/supplementary_tables.xlsx", overwrite = TRUE)

message("This script finished without errors")

