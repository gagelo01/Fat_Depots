#!/usr/bin/env Rscript
library(TwoSampleMR)
library(tidyverse)
library(data.table)
library(GagnonMR)
library(writexl)
library(furrr)

setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Fatdepots/")
all_inst_mvmr <- fread( "Data/Modified/all_inst_mvmr.txt")
all_outcome_mvmr <- fread( "Data/Modified/all_outcome_mvmr.txt" )
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")

##############
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ao <- fread("/mnt/sda/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao_small <- ao[id %in% list.files("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/"), ]
ao_small[, tomirge := tolower(id)]
options(future.globals.maxSize= 5e9)
plan(multicore, workers = 20, gc = TRUE)#plan(sequential)
#######################
#change only this section 
#univariable
idserver_exposure <-  df_index[grepl("trait-25", id) & grepl("adjbmi3", trait), id]
idserver_outcome <- c(paste0("dis-13-", 2:4), "dis-2-1", df_index[grepl("dis-6-",id), id], "dis-14-2", "dis-14-7",
                      df_index[grepl("trait-16", id) & grepl("^HDL|^LDL|^logTG", trait) & population %in% c("European", "Mixed"), id],
                      "trait-14-8", "trait-14-10", df_index[grepl("trait-2-", id), id],
                      "dis-7-1", "trait-12-2", "trait-13-1", "trait-13-2",
                      paste0("trait-27-", 1:2), paste0("trait-28-", 1:2))
idmrbase_exposure <- c(NULL)
idmrbase_outcome <- NULL
#multivariable
u<-c("logOR", "log odds", "log orr")
exposure <-c("trait-25-19")
exposurefor<-"trait-14-8"#The phenotype you want to correct for, but will be added as "exposure" which means there instruments will be selected
correctfor <- NULL #The phenotype you want to correct for , but will be added as "correctfor" which means there instruments won't selected
if(!is.null(exposurefor)){k <- purrr::cross2(exposure,exposurefor)}
k<-lapply(k, unlist)
mvmr_exposure <- k[sapply(k, function(x) (length(unique(x))==length(x)))]
mvmr_exposure <- unique(mvmr_exposure)
mvmr_outcome <- idserver_outcome
pval <- c(1)
if(is.null(correctfor)){correctfor<-"NULL"}
arguments_mvmr <- purrr::cross(list(mvmr_exposure, mvmr_outcome, correctfor, pval))
test<- sapply(arguments_mvmr, function(x) !any(grepl(paste(x[[3]], collapse = "|"), x[[1]])))
arguments_mvmr <- arguments_mvmr[test]
###############
###############
k <- c(df_index[unit %in% u,id],  ao_small[unit %in% u, id])
# exp_vec<-c(all_inst_mvmr[tolower(exposure) %in% tolower(idmrbase_exposure), unique(exposure)],
#            df_index[id %in% idserver_exposure, trait])
# out_vec<-c(all_inst_mvmr[tolower(exposure) %in% tolower(idmrbase_outcome), unique(exposure)],
#            df_index[id %in% idserver_outcome, trait])
# arguments_uni <- rbind(tidyr::crossing(exposure = exp_vec[!(exp_vec%in%k)], outcome = out_vec),
                       # tidyr::crossing(exposure = out_vec[!(out_vec%in%k)], outcome = exp_vec))
arguments_uni <- rbind(tidyr::crossing(id.exposure = idserver_exposure, id.outcome = idserver_outcome),
                       tidyr::crossing(id.exposure = idserver_outcome, id.outcome = idserver_exposure))
setDT(arguments_uni)
arguments_uni <- distinct(arguments_uni)
setDT(arguments_uni)
arguments_uni <- arguments_uni[!(id.exposure == id.outcome),]
big_index <- rbind(ao_small, df_index, fill = TRUE)
arguments_uni <- merge(arguments_uni, big_index[,.(id, sex)], by.x = "id.exposure", by.y = "id")
arguments_uni <- merge(arguments_uni, big_index[,.(id, sex)], by.x = "id.outcome", by.y = "id")
arguments_uni <- arguments_uni[sex.x == sex.y, ]


harm_univariate <- map(split(arguments_uni, 1:nrow(arguments_uni)), function(x) {
  harm <- TwoSampleMR::harmonise_data(exposure_dat = inst_all_sign_clump[id.exposure == x[,id.exposure], ],
                                      outcome_dat = all_outcome_mvmr[id.outcome == x[,id.outcome], ],
                                      action = 1)
  return(harm)}) %>% rbindlist(.,fill = TRUE)


setDT(harm_univariate)
harm_univariate[, exposure_outcome :=paste0(exposure, "_", outcome)]
egger_intercept <- mr_pleiotropy_test(harm_univariate)
setDT(harm_univariate)
harm_univariate <- harm_univariate[!(grepl("trait-2-", id.exposure) & SNP == "rs1121980"),] #in the FTO region

harm_univariate <- TwoSampleMR::steiger_filtering(harm_univariate) %>% as.data.table(.)
harm_univariate <- harm_univariate[mr_keep==TRUE,]
fwrite(harm_univariate, "Data/Modified/harm_univariate.txt")
harm_univariate <- harm_univariate[steiger_dir == TRUE,]

list_harm_univariate <- split(harm_univariate, harm_univariate$exposure_outcome)

list_res_univariate <- future_map(list_harm_univariate, function(x, idexp=idserver_exposure) {
  GagnonMR::all_mr_methods(x, short = FALSE,#!(unique(x$id.exposure)%in%idexp), 
                           skip_presso = !(unique(x$id.exposure)%in%idexp))}, #FALSE)}, 
  .options = furrr_options(seed = TRUE))

res_univariate <- rbindlist(list_res_univariate, fill = TRUE)

FandQ <- lapply(list_harm_univariate, function(x) {
  res <- TwoSampleMR::mr_heterogeneity(x)
  if(nrow(res)==0){res<-data.frame(exposure = x$exposure, outcome = x$outcome)}
  x <- TwoSampleMR::add_rsq(x)
  res$fstat<-GagnonMR::fstat_fromdat(list(x))
  res$rsq <- sum(x$rsq.exposure)
  setDT(x)
  res$N <- x[,.N]
  return(res)
}) %>% rbindlist(.,fill = TRUE)

######### Include biasprop
FandQ <- merge(FandQ, df_index[,.(id, sample_size)], by.x = "id.exposure", by.y = "id")
setnames(FandQ, "sample_size", "sample_size.exposure")
FandQ <- merge(FandQ, df_index[,.(id, sample_size)], by.x = "id.outcome", by.y = "id")
setnames(FandQ, "sample_size", "sample_size.outcome")
FandQ[,sample_size.exposure := as.numeric(sample_size.exposure)]
FandQ[, sample_size.outcome := as.numeric(sample_size.outcome)]
FandQ[,sample_overlap := as.numeric(sample_size.exposure)/as.numeric(sample_size.outcome)]
FandQ[,sample_overlap := ifelse(sample_overlap>1, 1/sample_overlap, sample_overlap)]
FandQ[, biasprop := sample_overlap / fstat]


#######
setDT(egger_intercept)

varclean <- distinct(df_index[id %in% inst_all_sign_clump$id.exposure,.(id, clean_variable_name)])
vecclean<- c(ALT = "Alanine transaminase", AST = "Aspartate transferase", ASATadjBMI = "ASATadj", 
             GfatadjBMI = "GFATadj", VATadjBMI = "VATadj", "fracture risk" = "Fracture")
varclean[clean_variable_name%in%names(vecclean),clean_variable_name := vecclean[clean_variable_name]]

list_res<-list(FandQ= FandQ, res_univariate=res_univariate, egger_intercept=egger_intercept)
list_res <- map(list_res, function(x) {
  x <- merge(x, df_index[,.(id,sex,population)], by.x = "id.exposure", by.y = "id")
  x <- merge(x, df_index[,.(id,sex,population)], by.x = "id.outcome", by.y = "id", suffixes = c(".exposure", ".outcome"))
  x <- merge(x, distinct(varclean), by.x = "id.exposure", by.y = "id")
  x <- merge(x, distinct(varclean), by.x = "id.outcome", by.y = "id", suffixes = c(".exposure", ".outcome"))
  return(x)
})

fwrite(list_res$FandQ, "Data/Modified/FandQ_univariate.txt")
fwrite(list_res$res_univariate, "Data/Modified/res_univariate.txt")
fwrite(list_res$egger_intercept, "Data/Modified/egger_intercept.txt")

#mvmr
#mvmr
# performmvmr <- function(exposure_vec, outcome_vec, correctfor=NULL,pval_threshold = 1,
#                         clump_exp_arg = "none", pval_inst = 5e-8,
#                         clump_r2 = 0.01, clump_kb = 1000) { #clump_exp_arg either none, first, or second
#   message(paste0("MVMR for the effect of ", paste(exposure_vec, collapse = " + "), " on ", outcome_vec, " while correcting for ",  paste(correctfor, collapse = " + ")))
#   exposure_dat <- inst_all_sign_clump[inst_all_sign_clump$id.exposure %in% exposure_vec & inst_all_sign_clump$pval.exposure < pval_inst,]
#   k<-all_inst_mvmr[(all_inst_mvmr$id.exposure %in% correctfor) & (all_inst_mvmr$SNP %in% unique(exposure_dat$SNP)),]
#   exposure_dat<- rbind(exposure_dat,k, fill = TRUE)
#   d1 <- all_inst_mvmr[(all_inst_mvmr$id.exposure %in% c(exposure_vec,correctfor)) & (all_inst_mvmr$SNP %in% unique(exposure_dat$SNP)),]
#   
#   if(clump_exp_arg == "none") {clump_exp<-NULL} else if(clump_exp_arg == "first"){clump_exp<-exposure_vec[1]} else if(clump_exp_arg == "second"){clump_exp<-exposure_vec[2]}
#   inst_mvmr <- prepare_for_mvmr(exposure_dat = exposure_dat, d1 =d1,clump_r2 = clump_r2, clump_kb = clump_kb, pval_threshold = pval_threshold, clump_exp = clump_exp, harmonise_strictness = 1)
#   
#   exposure_outcome_harmonized <- TwoSampleMR::mv_harmonise_data(exposure_dat = inst_mvmr,
#                                                                 outcome_dat = all_outcome_mvmr[all_outcome_mvmr$id.outcome == outcome_vec,],
#                                                                 harmonise_strictness = 1)
#   mvmr_results <- GagnonMR::mv_multiple_MendelianRandomization(exposure_outcome_harmonized = exposure_outcome_harmonized)
#   mvmr_results$clump_exposure <- mvmr_results$clump_exp %>% ifelse(is.null(.), "none", .)
#   mvmr_results <- mvmr_results[!(mvmr_results$exposure %in% correctfor),]
#   mvmr_results$correctfor <- paste(correctfor, collapse = " + ")
#   mvmr_results <- merge(distinct(all_outcome_mvmr[all_outcome_mvmr$id.outcome == outcome_vec,.(id.outcome, outcome)]), mvmr_results, by = "outcome")
#   mvmr_results <- merge(distinct(inst_mvmr[,.(id.exposure, exposure)]), mvmr_results, by = "exposure")
#   mvmr_results[ , other_exposures:=apply(.SD, 1, function(x) setdiff(exposure_vec,  x)), .SDcols = "id.exposure"]
#   return(mvmr_results)
# }
# 
# resmvmr <- future_map(arguments_mvmr, function(x) {
#   resnone<-performmvmr(exposure_vec = x[[1]],
#                        outcome_vec = x[[2]],
#                        correctfor=if(x[[3]]=="NULL"){NULL}else{x[[3]]},
#                        clump_exp_arg = "none",
#                        pval_inst = x[[4]])
#   return(resnone)
# }, .options = furrr_options(seed = TRUE))
# 
# names(resmvmr) <- sapply(arguments_mvmr, function(x) paste0(paste(x[[1]],collapse = " + "), " ~ ", x[[2]], " correctfor = ", x[[3]]," (pval=", x[[4]],")"))
# resmvmr <- rbindlist(resmvmr)
# fwrite(resmvmr, "Data/Modified/res_mvmr.txt")

message("This script finished without errors")
