#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(ckbplotr)
setwd("/home/gagelo01/workspace/Projects/small_MR_exploration/Fatdepots/")

format_right_var_forest <- function(input) {
  if(is.factor(input)) {
    message("format_right_var_forest does not accept factor; coercing to character")
    input<-as.character(input)
  }
  myvector<-vector(mode="integer", length = length(input))
  for(i in 1:length(input)) {
    input[is.na(input)]<-"dummy"
    if(i==1) {
      value = input[i]
    } else {

      if(input[i-1]!=input[i]){
        value = input[i]
      }  else {
        value = NA
      }

    }
    myvector[i]<-value
  }

  myvector[myvector=="dummy"]<-NA
  return(myvector)
}

df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
res_reverse <- fread("Data/Modified/res_univariate.txt")
zsexdiff <- fread("Data/Modified/zsexdiff.txt")



res_univariable <- res_reverse[grepl("trait-25-", id.exposure)]
dt<- rbind(
  data.table(nom = c(paste0("trait-2-", 1:6), paste0("dis-6-", 1:3)), value = "glucose_homeostasie"),
  data.table(nom = c(paste0("dis-13-", 2:4), "dis-4-1", "dis-5-1", "dis-14-1", paste0("trait-13-", 1:2)), value = "vascular"),
  data.table(nom = c("dis-2-1", "trait-14-8", paste0("trait-27-", 1:2)), value = "liver"),
  data.table(nom = c("dis-7-1", "trait-12-2"), value = "kidney"),
  data.table(nom = unique(res_univariable[grepl("trait-16-", id.outcome), id.outcome]), value = "Lipids"),
  data.table(nom =   paste0("trait-28-", 1:2), value = "Bone")
)
vec_system<-dt$value;names(vec_system)<-dt$nom
res_univariable[, system.outcome := vec_system[id.outcome]]

mycolour<-c(rgb(237,125, 49, maxColorValue = 255),
            rgb(114.5, 18, 84, maxColorValue = 255 ),
            rgb(8,161, 217, maxColorValue = 255))
res_univariable[, clean_variable_name.exposure := factor(clean_variable_name.exposure, levels = c("GFATadj", "ASATadj", "VATadj"))]
names(mycolour) <-  levels(res_univariable$clean_variable_name.exposure)
res_univariable[,colour := mycolour[clean_variable_name.exposure]]
res_univariable <- res_univariable[order(system.outcome,clean_variable_name.outcome, clean_variable_name.exposure, sex.outcome, method),]

#############
my_forest_fancy <- function(
data,
col.right = c("pval"),
col.right.heading = list(c("Effect (95% CI)", "P value")),
exponentiate = TRUE,
xlab="",
col.left,
col.left.heading,
col.toformat,
col.header,
yesleftcol = TRUE) {
stopifnot(all(c("b", "lci", "uci", "panel", "colour") %in% colnames(data)))
data[, estimate := round(b, digits =2) ]
data[, lci := round(lci, digits = 2)]
data[, uci :=  round(uci, digits = 2)]
data[,pval := formatC(pval, format = "e", digits = 1)%>%gsub("NA", "",.)]
# data[, variable := as.character(1:.N), by = "panel"]

colnom<-c("panel",  col.left, "estimate", "lci", "uci", col.right, col.header, "colour") %>% unique
data <- data[, .SD, .SDcols = colnom]

leftnorepeat<-col.left[1]
list_dat <- split(data, data$panel)

if(yesleftcol) {
  list_dat <- map(list_dat, function(doA) {
    k<-doA[!is.na(get(leftnorepeat)),unique(get(leftnorepeat))]
    list_data<-vector(mode = "list", length = length(k))
    for(i in 1:length(k)) {
      insert <- doA[NA,]
      insert[,panel:=doA$panel[1]]
      list_data[[i]]<-rbind(doA[get(leftnorepeat)==k[i],], insert)
    }

    doA<-rbindlist(list_data)
    doA<- doA[1:(.N-1),]

    ###format
    doA[, (col.toformat) := lapply(.SD, format_right_var_forest) , .SDcols = col.toformat]

    return(doA)
  })
}


list_dat <- map(list_dat, function(doA) {
  doA[,variable:=1:.N]
  doA[,dummy:=as.character(NA)]
  })

mylabels <- data.frame(heading1 = as.character(list_dat[[1]][,get(col.header["heading1"])]),
                       heading2 = as.character(list_dat[[1]][,get(col.header["heading2"])]),
                       heading3 = as.character(list_dat[[1]][,get(col.header["heading3"])]),
                       variable = as.character(list_dat[[1]]$variable))

k<-make_forest_plot(panels = list_dat,
                    col.key = "variable",
                    row.labels = mylabels,
                    exponentiate = exponentiate,
                    pointsize = 2,
                    rows = unique(mylabels$heading1),
                    col.stderr = NULL,
                    col.lci = "lci",
                    col.uci = "uci",
                    col.right = col.right,
                    col.right.heading = col.right.heading,
                    col.left         = rev(col.left),
                    col.left.heading = list(rev(col.left.heading)),
                    xlab = xlab,#"NAFLD risk per 1 SD \n higher metabolomic measure",
                    blankrows = c(0,0,0,1),
                    colour = "colour", #"colour"
                    fill = "colour",
                    panel.headings = levels(data$panel),
                    scalepoints = FALSE,
                    envir = environment(),
                    # cicolour = "grey50",
                    col.left.hjust = c(0.5, rep(0, length(col.left)-1)),
                    col.right.hjust = 0.5,
                    shape = 22,
                    ciunder = TRUE,
                    stroke = 0,
                    nullval = ifelse(exponentiate == TRUE, 1, 0))#ifelse(exponentiate == TRUE, 1, 0))

return(k)
}

#disease
k<- res_univariable[sex.exposure=="Males and Females",][method == "Inverse variance weighted",][id.outcome %in% df_index[id%in%res_univariable$id.outcome & unit == "log odds", id],][!grepl("dis-5-1|dis-14-7", id.outcome),]
k[,panel:="a"]
k[,dummy:=as.character("NA")]
my_forest_fancy(
  data=k,
  col.right = c("pval"),
  col.right.heading = list(c("Effect (95% CI)", "P value")),
  exponentiate = TRUE,
  xlab="",
  col.left = c("clean_variable_name.outcome", "clean_variable_name.exposure", "nsnp"),
  col.left.heading=c("Outcomes (Diseases)","Exposures", "n SNP"),
  col.toformat = c("clean_variable_name.outcome"),
  col.header = c(heading1="dummy", heading2="dummy", heading3="dummy"))

ggsave(paste0("Results/", "figure_adipo_disease", ".png"),
       width=630/72,height=340/72, units="in", scale=1, dpi = 300,
       device = "png")
#Traits
k<- res_univariable[sex.exposure=="Males and Females",][method == "Inverse variance weighted",][!(id.outcome %in% df_index[id%in%res_univariable$id.outcome & unit == "log odds", id]),][!grepl("dis-5-1|trait-14-10", id.outcome),]
k <- k[!(grepl("trait-16", id.outcome) & population.outcome == "European"), ]
k[,panel:="a"]
k[,dummy:=as.character("NA")]
# debugonce(my_forest_fancy)
my_forest_fancy(
  data=k,
  col.right = c("pval"),
  col.right.heading = list(c("Effect (95% CI)", "P value")),
  exponentiate = FALSE,
  xlab="",
  col.left = c("clean_variable_name.outcome", "clean_variable_name.exposure", "nsnp"),
  col.left.heading=c("Outcomes (Traits)","Exposures", "n SNP"),
  col.toformat = c("clean_variable_name.outcome"),
  col.header = c(heading1="dummy", heading2="dummy", heading3="dummy"))

ggsave(paste0("Results/", "figure_adipo_traits", ".png"),
       width=609/72,height=623/72, units="in", scale=1, dpi = 300,
       device = "png")
###sex specific
k <- res_univariable[method == "Inverse variance weighted",][grepl("dis-6-|trait-16-|trait-2-|dis-13-", id.outcome),]
k <- k[!(grepl("trait-16", id.outcome) & population.outcome == "European"),]
k[,sex:=c("Females"="Women", "Males" = "Men", "Males and Females"="All")[sex.outcome]]
k <- k[sex %in% c("Women", "Men"),]
k[,sex:=factor(sex, levels = c("Men", "Women"))]
k[,panel:="a"]
k[,dummy:=as.character(NA)]
k <- merge(k, zsexdiff, by = c("clean_variable_name.exposure", "clean_variable_name.outcome","population.exposure", "population.outcome", "method"))
k[,sexdiff_z_stat_pval := formatC(sexdiff_z_stat_pval, format = "e", digits = 1)%>%gsub("NA", "",.)]
k[, clean_variable_name.exposure := factor(clean_variable_name.exposure, levels = c("GFATadj", "ASATadj", "VATadj"))]
k <- k[order(system.outcome,clean_variable_name.outcome, clean_variable_name.exposure, sex.outcome, method),]
mycolour<-c("#5884E5", "#9E131E")
names(mycolour) <-  levels(k$sex)
k[,colour := mycolour[sex]]

##sex specific

my_forest_fancy(
# data=k,
data = k[grepl("dis-", id.outcome)],
col.right = c("pval", "sexdiff_z_stat_pval"),
col.right.heading = list(c("Effect (95% CI)", "P value", "P value \n(sex difference)")),
exponentiate = FALSE,
xlab="",
col.left = c("clean_variable_name.outcome", "clean_variable_name.exposure", "sex", "nsnp"),
col.left.heading=c("Outcomes","Exposures", "Sex", "n SNP"),
col.toformat = c("clean_variable_name.outcome", "clean_variable_name.exposure", "sexdiff_z_stat_pval"),
col.header = c(heading1="dummy", heading2="dummy", heading3="dummy"))

ggsave(paste0("Results/", "figure_adipo_sexspecific", ".png"),
       width=630/72,height=340/72, units="in", scale=1, dpi = 300,
       device = "png")

message("This script finished without errors")
