#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(ckbplotr)
setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Fatdepots/")

df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
res_univariable <- fread("Data/Modified/res_univariate.txt")
FandQ <- fread("Data/Modified/FandQ_univariate.txt")
egger_intercept <- fread("Data/Modified/egger_intercept.txt")
#########
Z_test_twosided <- function(mean.x, mean.y, se.x, se.y) {
  
  # pooled standard error
  pooled_se <- sqrt((se.x^2 + se.y^2)/2)
  
  # calculate z-score
  z_score <- (mean.x - mean.y)/pooled_se
  
  # two-tailed p-value
  p_value <- 2*pnorm(-abs(z_score))
  
  return(data.table(z_stat = z_score,
                    z_stat_pval = p_value))
}

##########
k <- res_univariable[sex.exposure%in%c("Males", "Females"),]
k<- k[,.(clean_variable_name.exposure, clean_variable_name.outcome, sex.exposure, population.exposure, population.outcome, method, b, se, pval)]
k <- dcast(k, clean_variable_name.exposure + clean_variable_name.outcome + method + population.exposure + population.outcome ~ sex.exposure, value.var = c("b", "se", "pval"))
list_dat<-split(k, 1:k[,.N])

zsexdiff <-map(list_dat, function(x){
res<- x[,Z_test_twosided(mean.x = b_Females, mean.y = b_Males, se.x = se_Females, se.y = se_Males )]
setnames(res, colnames(res), paste0("sexdiff_", colnames(res)))
res  <- cbind(x, res)
return(res)}) %>% rbindlist(., fill = TRUE)

fwrite(zsexdiff, "Data/Modified/zsexdiff.txt")

######

# k<- res_univariable[grepl("trait-25", id.exposure)]
# k<-k[!method %in% c("MR Egger","Weighted mode") & sex.exposure == "Males and Females",  ]
# mmm <- k[ , all(pval<0.05), by=c("id.exposure", "id.outcome", "exposure", "outcome")]
# mmm[exposure=="mri_vatadjbmi3" & V1==TRUE ,]
# 
# 
# k<- res_univariable[grepl("trait-25", id.outcome)]
# k<-k[!method %in% c("MR Egger","Weighted mode") & sex.exposure == "Males and Females",  ]
# mmm <- k[ , all(pval<0.05)&(all(b<0)|all(b>0)), by=c("id.exposure", "id.outcome", "exposure", "outcome")]
# k <- merge(k, mmm, by = c("id.exposure", "id.outcome", "exposure", "outcome"))
# k[V1==TRUE,][,.(exposure, outcome, method, nsnp, b, se, pval)]
# mmm[V1==TRUE,]
# mmm[exposure=="mri_vatadjbmi3" & V1==TRUE,]
