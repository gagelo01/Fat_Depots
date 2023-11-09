#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Fatdepots/")

df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
res_univariable <- fread("Data/Modified/res_univariate.txt")
zsexdiff <- fread( "Data/Modified/zsexdiff.txt")
harm_univariate <- fread( "Data/Modified/harm_univariate.txt")
inst_all_sign_clump <- fread( "Data/Modified/inst_all_sign_clump.txt")
FandQ <- fread("Data/Modified/FandQ_univariate.txt")

#########
res_univariable <- res_univariable[!(grepl("trait-16", id.outcome) & population.outcome == "European"), ]
res_univariable <- res_univariable[!(grepl("trait-16", id.exposure) & population.exposure == "European"), ]

#####prepare data
sexdifftrue<-zsexdiff[!(method %in% c("MR Egger","Weighted mode", "Contamination mixture")),]
sexdifftrue <- sexdifftrue[, all(sexdiff_z_stat_pval<0.05), by = c("clean_variable_name.exposure", "clean_variable_name.outcome", "population.exposure", "population.outcome") ]
setnames(sexdifftrue, "V1", "is_different")
zsexdiff <- merge(zsexdiff, sexdifftrue, by = c("clean_variable_name.exposure", "clean_variable_name.outcome", "population.exposure", "population.outcome"))

k<-res_univariable[!method %in% c("MR Egger","Weighted mode"),  ]
k <- k[ , all(pval<0.05)&(all(b<0)|all(b>0)), by=c("id.exposure", "id.outcome", "exposure", "outcome")]
setnames(k, "V1", "is_causal")
res_univariable <- merge(res_univariable, k, by = c("id.exposure", "id.outcome", "exposure", "outcome"))

######
return_format_data<-function(data) {
  return(data[, paste0(format(round(exp(b), digits = 2), nsmall = 2), " 95% CI=", format(round(exp(lci), digits = 2), nsmall = 2), "-",  format(round(exp(uci), digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))])
}
return_format_data_noexp <-function(data) {
  return(data[, paste0(format(round(b, digits = 2), nsmall = 2) , " 95% CI=", format(round(lci, digits = 2), nsmall = 2), "-",  format(round(uci, digits = 2), nsmall = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))])
}
return_format_fstat <-function(data) {
  k<-data[, paste0(N, " SNPs (r2 = ", round(rsq*100, digits =2), "%; F-statistics = ",  round(fstat, digits = 0), ")")]
  names(k) <- data[, paste0(exposure, " on ", outcome)]
  return(k)
}

#### Abstract
df_index[id %in% res_univariable[grepl("trait-25", id.exposure), id.exposure], max(sample_size)]
k<- res_univariable[grepl("trait-25", id.exposure)]
k<-k[!method %in% c("MR Egger","Weighted mode") & sex.exposure == "Males and Females",  ]
mmm <- k[ , all(pval<0.05), by=c("id.exposure", "id.outcome", "exposure", "outcome")]
mmm[exposure=="mri_gfatadjbmi3" & V1==TRUE ,]
mmm[exposure=="mri_asatadjbmi3" & V1==TRUE ,]
mmm[exposure=="mri_vatadjbmi3" & V1==TRUE ,]

# Results
#para 1
FandQ[grepl("trait-25-", id.exposure), min(sample_overlap)]
FandQ[grepl("trait-25-", id.exposure), max(sample_overlap)]
FandQ[grepl("trait-25-", id.exposure), min(biasprop)]
FandQ[grepl("trait-25-", id.exposure), max(biasprop)]
0.05/54
#para 2
FandQ[grepl("trait-25", id.exposure),min(fstat)]
# para 2
dttrait<- res_univariable[sex.exposure=="Males and Females",][!(id.outcome %in% df_index[id%in%res_univariable$id.outcome & unit == "log odds", id]),][!grepl("dis-5-1", id.outcome),]

myexp<-"mri_gfatadjbmi3"
myid <- res_univariable[exposure==myexp, ][which.max(nsnp),id.outcome]
FandQ[id.outcome == myid & exposure == myexp, ] %>% return_format_fstat
dttrait[is_causal == TRUE & exposure == myexp, unique(clean_variable_name.outcome)]
dttrait[exposure == myexp & method == "Inverse variance weighted", mean(abs(b))]


myexp<-"mri_asatadjbmi3"
myid <- res_univariable[exposure==myexp, ][which.max(nsnp),id.outcome]
FandQ[id.outcome == myid & exposure == myexp, ] %>% return_format_fstat
dttrait[is_causal == TRUE & exposure == myexp, unique(clean_variable_name.outcome)]
dttrait[exposure == myexp & method == "Inverse variance weighted", mean(abs(b))]

myexp<-"mri_vatadjbmi3"
myid <- res_univariable[exposure==myexp, ][which.max(nsnp),id.outcome]
FandQ[id.outcome == myid & exposure == myexp, ] %>% return_format_fstat
dttrait[exposure == myexp & method == "Inverse variance weighted", mean(abs(b))]

dttrait[is_causal == TRUE & exposure == myexp, unique(clean_variable_name.outcome)]
dttrait[is_causal == TRUE & exposure == myexp & method == "Inverse variance weighted" & outcome == "LDL_GLGC2022_Mixed_MalesandFemales", ]  %>% return_format_data_noexp(.)
dttrait[is_causal == TRUE & exposure == myexp & method == "Inverse variance weighted" & outcome == "aspartateaminotransferase_AST", ]  %>% return_format_data_noexp(.)

dttrait[grepl("trait-25", id.exposure), all(is_causal==FALSE),by = "outcome"][V1==TRUE,]
#para 3
dtdis<- res_univariable[sex.exposure=="Males and Females",][id.outcome %in% df_index[id%in%res_univariable$id.outcome & unit == "log odds", id],][!grepl("dis-5-1", id.outcome),]
dtdis[is_causal==TRUE, .(exposure, outcome)] %>% distinct
dtdis[exposure == "mri_vatadjbmi3" & method == "Inverse variance weighted"]
dtdis[exposure == "mri_vatadjbmi3" & method == "Inverse variance weighted" & outcome == "Mahajan_Type2diabetes",] %>% return_format_data(.)
dtdis[exposure == "mri_vatadjbmi3" & method == "Inverse variance weighted" & outcome == "cad_aragam_combined"]%>% return_format_data(.)

# para 4 infering the direction of causality
vecrev <- res_univariable[grepl("trait-25", id.outcome),][is_causal==TRUE,][,.(id.exposure, id.outcome)] %>% distinct(.)

propsteiger <- harm_univariate[, sum(!steiger_dir)/.N, by = c("id.exposure", "id.outcome","exposure","outcome")][order(V1),]
setnames(propsteiger, "V1", "steigerprop")
propsteiger[grepl("trait-25",id.exposure), mean(steigerprop)]
propsteiger[grepl("trait-25",id.outcome), mean(steigerprop)]


dtreverse<-merge(vecrev, propsteiger, by = c("id.exposure", "id.outcome"))
dtreverse[exposure!="Fat_Liver", ][order(steigerprop),]

# para 5 sex specific effects
df_index[grepl("trait-25-", id) & id %in% res_univariable$id.exposure, max(sample_size), by = "sex"]
k<- FandQ[grepl("trait-25-",id.exposure) & sex.exposure %in% c("Males", "Females") & sex.exposure==sex.outcome, .SD[which.max(N)],by = c("exposure")]
k %>% return_format_fstat


zsexdiff[method == "Inverse variance weighted" & grepl("adjBMI", clean_variable_name.exposure) & is_different == TRUE,]


k <- zsexdiff[, c("clean_variable_name.exposure", "clean_variable_name.outcome", "population.exposure", "population.outcome", "is_different")] %>% distinct()
k <- merge(k, res_univariable, by = c("clean_variable_name.exposure", "clean_variable_name.outcome", "population.exposure", "population.outcome"))

k[is_causal == TRUE & is_different == TRUE & sex.exposure %in% c("Males", "Females"),][]
k[is_different==TRUE, .(clean_variable_name.exposure, clean_variable_name.outcome)] %>% distinct(.)

zsexdiff[is_different == TRUE & clean_variable_name.exposure %in% c("GfatadjBMI", "ASATadjBMI", "VATadjBMI"),] 
k[is_different==TRUE & clean_variable_name.exposure == "ASATadjBMI",]
k[is_different==TRUE & clean_variable_name.exposure == "VATadjBMI" & method == "Inverse variance weighted",]
k[is_different==TRUE & is_causal==TRUE & method == "Inverse variance weighted",]

##Multivariable MR
res_univariate_x <- fread(res_univariate, "Data/Modified/reviewer_response_res_univariate.txt")
res_mvmr_x <- fread("Data/Modified/reviewer_response_res_mvmr.txt")
res_mvmr_x[id.exposure=="trait-25-10", ] %>% return_format_data
################

#Methods

################
#Para 1
df_index[grepl("trait-25", id), max(sample_size), by = "sex"]
df_index[grepl("trait-25", id), unique(note)]
dtdis$clean_variable_name.outcome %>% unique %>% length
dttrait[outcome!="Fat_Liver"]$clean_variable_name.outcome %>% unique(.) %>% length(.)
res_univariable[grepl("trait-25", id.exposure), .(id.exposure, id.outcome)] %>% distinct(.) %>% .[,.N]


################

#Discussion

################
#Para 4
k <- res_univariable[sex.exposure == "Males and Females" & clean_variable_name.exposure %in% c("GfatadjBMI","VATadjBMI"),]
k <- k[method=="Inverse variance weighted",.(exposure, outcome, b)]
setnames(k, "b","beta")
k <- dcast(data = k,formula = outcome~exposure, value.var = "beta")
k[,absolute_diff := abs(mri_gfatadjbmi3) - abs(mri_vatadjbmi3)]
k
