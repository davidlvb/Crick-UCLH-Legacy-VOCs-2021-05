# 
# Crick/UCLH Legacy Cohort Neutralisation Data Analysis
# Against SARS-CoV-2 Variants of Concern B.1.617.2, B.1.351, and B.1.1.7
#
# 15 May 2021
#
# David LV Bauer
# RNA Virus Replication Laboratory
# The Francis Crick Institute
#
# This work is licensed under a CC BY 4.0 license and is free to use with attribution
# 

### Remove all stored data from environment & reset plot ###
rm(list = ls())
dev.off()

### Load Legacy Study Data table ###
load("Crick_Legacy_2021-24-05_B-1-617-2_PUBLIC.Rda")

### Load required packages ### 
# library(plotly)
library(tidyverse)
library(readxl)
library(khroma)
library(boot)
library(extrafont)
library(svglite)
library(ggpubr)
library(stats)
library(cowplot)
library(broom)
library(furniture)
library(rms)
# source("geom_split_violin.R")

### Subset data for further analysis
dtHashed <- dtHashed %>% 
  filter(COVID_vaccStatus %in% c(1,2)) # Ignore individuals who are in seroconversion window following dose 1 and dose 2 of vaccine
dt <- dtHashed


### Set constants for various functions / plots
runBootstrapStats <- FALSE  # Set to TRUE in order to calculate boostrap stats
runFontImport <- FALSE # Set to TRUE in order to import fonts into R environment
strainOrder <- c("Wuhan1", "D614G", "Kent", "SAfrica", "India2")
referenceIC50 <- 2^9.802672



######################################################
#          Set up functions for later use            #
######################################################

if(runFontImport){
  font_import(path = "fonts/", prompt = FALSE)
}

boot_foldMedian <- function(d, i){
  d2 <- d[i,]
  firstCol <- d2 %>% pull(2)
  secondCol <- d2 %>% pull(3)
  foldMedian <- (median(firstCol, na.rm=TRUE) / median(secondCol, na.rm=TRUE))
  return(foldMedian)
}

CIbootstrap <- function(testTable, strainComparison, rpt=5000) {
  statcheck <- testTable %>% filter(strain %in% strainComparison) %>% pivot_wider(values_from = ic50, names_from = strain)
  results_boot <- boot(statcheck, boot_foldMedian, R=rpt)
  results_ci <- boot.ci(results_boot, type="basic")
  print(results_boot)
  print(results_ci)
  print("######################################################")
  print(strainComparison)
  print(c("Median1", "Median2", "Fold-decrease"))
  print(c(median(pull(statcheck, 2), na.rm=TRUE),
          median(pull(statcheck, 3), na.rm=TRUE),
          median(pull(statcheck, 2), na.rm=TRUE)/ median(pull(statcheck, 3), na.rm=TRUE) )
  )
}

pTtest12 <- function(relevantData, strainSel="Wuhan1", resultTable = 0){
  testTable <- relevantData %>% filter(strain ==strainSel) %>% ungroup() %>% 
    select(sample_barcode, bc_visit, ic50, ic50_log2) %>%
    mutate(bc_visit=replace(bc_visit, bc_visit > 1, 2))
  visit1 <- testTable %>% filter(bc_visit == 1) %>% pull(ic50_log2)
  visit2 <- testTable %>% filter(bc_visit == 2) %>% pull(ic50_log2)
  broomTable <- t.test(visit1, visit2, paired = TRUE) %>% tidy() %>% add_column(strain=strainSel, .before=1)
  if (is.numeric(resultTable)){
    return(broomTable)
  } else {
    return( bind_rows (resultTable, broomTable))
  }
}

pWtest12 <- function(relevantData, strainSel="Wuhan1", resultTable = 0){
  testTable <- relevantData %>% filter(strain ==strainSel) %>% ungroup() %>% 
    select(sample_barcode, bc_visit, ic50, ic50_log2) %>%
    mutate(bc_visit=replace(bc_visit, bc_visit > 1, 2))
  visit1 <- testTable %>% filter(bc_visit == 1) %>% pull(ic50_log2)
  visit2 <- testTable %>% filter(bc_visit == 2) %>% pull(ic50_log2)
  broomTable <- wilcox.test(visit1, visit2, paired = TRUE) %>% tidy() %>% add_column(strain=strainSel, .before=1)
  if (is.numeric(resultTable)){
    return(broomTable)
  } else {
    return( bind_rows (resultTable, broomTable))
  }
}






######################################################
#    Table S1. Generate table of participant data    #
######################################################

studyData <- dtHashed
studyData$shortDoseInterval <- studyData$COVID_daysBetweenJabs < 50
dose1cohort <- studyData %>% filter(COVID_vaccStatus == 1,  sampleOrderInVaccStatus == 1)
dose2cohort <- studyData %>% filter(COVID_vaccStatus == 2,  sampleOrderInVaccStatus == 1)
doseBothCohort <- studyData %>% filter(bc_participant_id %in% dose1cohort$bc_participant_id, 
                                       bc_participant_id %in% dose2cohort$bc_participant_id, 
                                       COVID_vaccStatus == 2,
                                       sampleOrderInVaccStatus == 1)

participants <- studyData %>% group_by(bc_participant_id) %>% arrange(collection_date) %>%  slice_head(n=1)
participants$inDose1cohort <- participants$bc_participant_id %in% dose1cohort$bc_participant_id
participants$inDose2cohort <- participants$bc_participant_id %in% dose2cohort$bc_participant_id

participants$sex <- factor(participants$sex)
participants$site <- factor(participants$site)
participants$COVID_sumJabs_niceLabels <- c("One dose", "Two doses")[participants$inDose1cohort]

participants <- participants %>% ungroup()


tab1 <- table1(ungroup(dose1cohort),
              site,
              age,
              sex,
              BMI,
              ethnicity_3,
              splitby = "shortDoseInterval",
              test = TRUE,
              # assume parametric distribution:
              param = TRUE, 
              na.rm = TRUE,
              output = "text2", export = "supplementary_table_1dose")
tab2 <- table1(ungroup(dose2cohort),
               site,
               age,
               sex,
               BMI,
               ethnicity_3,
               splitby = "shortDoseInterval",
               test = TRUE,
               # assume parametric distribution:
               param = TRUE, 
               na.rm = TRUE,
               output = "text2", export = "supplementary_table_2dose")

tabB <- table1(ungroup(participants),
               site,
               age,
               sex,
               BMI,
               ethnicity_3,
               test = TRUE,
               # assume parametric distribution:
               param = TRUE, 
               na.rm = TRUE,
               output = "text2", export = "supplementary_table_Bdose")

# Test if difference in age or BMI between cohorts
t.test(dose1cohort$age, dose2cohort$age)
t.test(dose1cohort$BMI, dose2cohort$BMI)

# Determine median days since dose 2
quantile(dose2cohort$COVID_daysSinceJab2, probs = 0.95)
median(dose2cohort$COVID_daysSinceJab2)
fivenum(dose2cohort$COVID_daysSinceJab2)
as.integer(dose2cohort$COVID_daysSinceJab2) %>% gghistogram(fill = "#BBCCEE", 
                                                            xlab = "Days Since Dose 2", 
                                                            ylab="Number of Participants",
                                                            binwidth = 7,
                                                            add = "median"
                                                            )

# Narrative in manuscript...
# Using a high-throughput live-virus SARS-CoV-2 neutralisation assay, we determined 
#  NAb titres (NAbTs) in 250 participants (median age = 42 years, IQR = 33-52 years)
nrow(participants)
fivenum(participants$age)
# following either 1 dose (n = 186, median days post-1st dose = 36, IQR = 23-38 days)
nrow(dose1cohort)
fivenum(dose1cohort$COVID_daysSinceJab1)
# and 2 doses (n= 195, median days post-2nd dose = 28, IQR = 21-37 days)
nrow(dose2cohort)
fivenum(dose2cohort$COVID_daysSinceJab2)
# of BNT162b2 (Pfizer-BioNTech) against five SARS-CoV-2 strains.

# Summarise ELISA Status for each cohort
dose1cohort %>% group_by(elisa_result) %>% summarise (n())
dose2cohort %>% group_by(elisa_result) %>% summarise (n())


########################################################################
#   Panel 1. Vaccine responses per strain following 2nd dose Pfizer    #
########################################################################

relevantData <- dose2cohort %>%
  pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50")
relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)

relevantData %>% group_by(strain) %>% summarise(counts = n(), qtile = c(0.25, 0.5, 0.75),  value = quantile(ic50, c(0.25, 0.5, 0.75), na.rm=TRUE))


outplot <- ggplot(relevantData, aes(x=strain, y=ic50, color = strain, label = sample_barcode)) + 
  scale_colour_muted() +
  # scale_color_brewer(palette="Set1") + 
  geom_hline(yintercept = referenceIC50,linetype = 3 ) +
  geom_violin(trim=TRUE) + 
  # stat_summary(fun.data = give.n, geom = "text", fun.y = median,
  #              position = position_dodge(width = 0.75)) +
  # scale_y_log10(limits = c(0.1, 2^12))
  scale_y_continuous(trans='log2', 
                     breaks=c(5, 10, 64, 256, 1024, 5120), 
                     labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
                     minor_breaks=c(32, 128, 512, 2048)) +
  ylab(bquote('Virus Neutralisation, '~IC[50]~ '')) +
  geom_jitter(shape=20, position=position_jitter(0.2), alpha=0.3) + 
  # facet_wrap(~ COVID_vaccStatus, labeller = label_both, nrow = 1) +
  stat_summary(fun=median, geom = "point", color="black",  shape=5, size=1, stroke=1) + 
  theme_bw(base_family = "Helvetica Neue Thin") +
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + 
  theme(axis.text.y = element_text(size=12))  + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y=element_text(size=15),
    axis.title.x=element_blank()
  )

# outplot

Panel1 <- outplot

ggsave("FIGURE-Panel1.svg", outplot, width=50, height=15, units="cm")

# ggplotly(outplot)

### Stats for Panel 1 ###
grouped <- relevantData %>% group_by(strain)
summaryTable <- grouped %>% summarise(count=n(), median=median(ic50, na.rm = T))
summaryTable$dMedWu1 <- summaryTable$median[1] / summaryTable$median
summaryTable$dMedD61 <- summaryTable$median[2] / summaryTable$median
summaryTable$dMedKnt <- summaryTable$median[3] / summaryTable$median
summaryTable$dMedSAf <- summaryTable$median[4] / summaryTable$median
summaryTable$dMedIn2 <- summaryTable$median[5] / summaryTable$median


summaryTable


testTable <- relevantData %>% ungroup %>% select(sample_barcode, strain, ic50) 
pairwise.wilcox.test(testTable$ic50, testTable$strain, paired = TRUE) 

if (runBootstrapStats == TRUE){
  # CIbootstrap(testTable, c("Wuhan1", "D614G"))
  CIbootstrap(testTable, c("Wuhan1", "Kent"))
  CIbootstrap(testTable, c("Wuhan1", "SAfrica"))
  CIbootstrap(testTable, c("Wuhan1", "India2"))
  
  CIbootstrap(testTable, c("D614G", "Kent"))
  CIbootstrap(testTable, c("D614G", "SAfrica"))
  CIbootstrap(testTable, c("D614G", "India2"))
  
  CIbootstrap(testTable, c("SAfrica", "India2"))
}


########################################################################
#   Panel 1firstdoseonly. Vaccine responses per strain following 1st dose Pfizer    #
########################################################################


relevantData <- dose1cohort %>%
  pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50")
relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)

outplot <- ggplot(relevantData, aes(x=strain, y=ic50, color = strain, label = sample_barcode)) + 
  scale_colour_muted() +
  # scale_color_brewer(palette="Set1") + 
  geom_hline(yintercept = referenceIC50,linetype = 3 ) +
  geom_violin(trim=TRUE) + 
  # stat_summary(fun.data = give.n, geom = "text", fun.y = median,
  #              position = position_dodge(width = 0.75)) +
  # scale_y_log10(limits = c(0.1, 2^12))
  scale_y_continuous(trans='log2', 
                     breaks=c(5, 10, 64, 256, 1024, 5120), 
                     labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
                     minor_breaks=c(32, 128, 512, 2048)) +
  ylab(bquote('Virus Neutralisation, '~IC[50]~ '')) +
  geom_jitter(shape=20, position=position_jitter(0.2), alpha=0.3) + 
  # facet_wrap(~ COVID_vaccStatus, labeller = label_both, nrow = 1) +
  stat_summary(fun=median, geom = "point", color="black",  shape=5, size=1, stroke=1) + 
  theme_bw(base_family = "Helvetica Neue Thin") +
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12)) + 
  theme(axis.text.y = element_text(size=12))  + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y=element_text(size=15),
    axis.title.x=element_blank()
  )

# outplot

Panel1firstdoseonly <- outplot

ggsave("FIGURE-Panel1firstdoseonly.svg", outplot, width=50, height=15, units="cm")

# ggplotly(outplot)

### Stats for Panel 1 ###
grouped <- relevantData %>% group_by(strain)
summaryTable <- grouped %>% summarise(count=n(), median=median(ic50, na.rm = T))
summaryTable$dMedWu1 <- summaryTable$median[1] / summaryTable$median
summaryTable$dMedD61 <- summaryTable$median[2] / summaryTable$median
summaryTable$dMedKnt <- summaryTable$median[3] / summaryTable$median
summaryTable$dMedSAf <- summaryTable$median[4] / summaryTable$median
summaryTable$dMedIn2 <- summaryTable$median[5] / summaryTable$median


summaryTable


testTable <- relevantData %>% ungroup %>% select(sample_barcode, strain, ic50) 
pairwise.wilcox.test(testTable$ic50, testTable$strain, paired = TRUE) 

if (runBootstrapStats == TRUE){
  # CIbootstrap(testTable, c("Wuhan1", "D614G"))
  CIbootstrap(testTable, c("Wuhan1", "Kent"))
  CIbootstrap(testTable, c("Wuhan1", "SAfrica"))
  CIbootstrap(testTable, c("Wuhan1", "India2"))
  
  CIbootstrap(testTable, c("D614G", "Kent"))
  CIbootstrap(testTable, c("D614G", "SAfrica"))
  CIbootstrap(testTable, c("D614G", "India2"))
  
  CIbootstrap(testTable, c("SAfrica", "India2"))
}




########################################################################
#      Panel 2a. Correlation of IC50 between virus variants           #
########################################################################

relevantData <- dose2cohort  %>%
  rename_with(~str_replace(., "_ic50", "")) 


dataset <- c("D614G", "Kent", "SAfrica", "India2")
breakList <- c(5, 10, 64, 256, 1024, 5120)

outplot <- ggscatter(relevantData, x = "Wuhan1", y = dataset,
                     combine = TRUE, ylab = "Variant IC50", xlab = "Wuhan1 IC50", size=2, color=alpha("black", alpha=0.3), shape=20,
                     add = "reg.line", conf.int = TRUE, 
                     add.params = list(color = "#CC6677", fill = alpha("#CC6677", 0.2)),
                     nrow=1,
) +
  scale_y_continuous(trans='log2', 
                     breaks=breakList, 
                     labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
                     minor_breaks=c(32, 128, 512, 2048)
                     ) +
  ylab(bquote('Variant '~IC[50]~ '')) +
  scale_x_continuous(trans='log2', 
                     breaks=breakList, 
                     labels=c("", "[<40]", "64", "", "1024", ""),
                     limits = c(breakList[1], tail(breakList, n = 1)),
                     minor_breaks=c(32, 128, 512, 2048)   # , limits = c(min(breakList), max(breakList))
                     )  +
  xlab(bquote('Wuhan '~IC[50]~ '')) +
  stat_cor(method = "spearman", label.sep="\n", label.y.npc = 0.93, size=3) + 
  theme_bw(base_family = "Helvetica Neue Thin") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank()
  )

# outplot
Panel2a <- outplot
ggsave("FIGURE-Panel2a.svg", outplot, width=15, height=5, units="cm")




########################################################################
#      Panel 2b. Correlation of IC50 between virus variants           #
########################################################################

relevantData <- dose2cohort %>%
  rename_with(~str_replace(., "_ic50", "")) 


breakList <- c(5, 10, 64, 256, 1024, 5120)

outplot <- ggscatter(relevantData, x = "India2", y = "SAfrica",
                     combine = TRUE, ylab = "", xlab = "India2 IC50", size=2, color=alpha("black", alpha=0.3), shape=20,
                     add = "reg.line", conf.int = TRUE, 
                     add.params = list(color = "#88CCEE", fill = alpha("#88CCEE", 0.2)),
                     nrow=1,
) +  
  scale_y_continuous(trans='log2', 
                     breaks = NULL,
                     labels = NULL) +
                     # breaks=breakList, 
                     # labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
                     # minor_breaks=c(32, 128, 512, 2048)) +
  scale_x_continuous(trans='log2', 
                     breaks=breakList, 
                     labels=c("", "[<40]", "64", "", "1024", ""),
                     limits = c(breakList[1], tail(breakList, n = 1)),
                     minor_breaks=c(32, 128, 512, 2048)   # , limits = c(min(breakList), max(breakList))
  )  +
  xlab(bquote('India2 '~IC[50]~ '')) +
  theme_bw(base_family = "Helvetica Neue Thin") +
  ggtitle("SAfrica") +
  theme(plot.title = element_text(size=8.8, hjust=0.5))+
  stat_cor(method = "spearman", label.sep="\n", label.y.npc = 0.93, size=3) + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_blank(),
    strip.background = element_blank()
  )

# outplot
Panel2b <- outplot
ggsave("FIGURE-Panel2b.svg", outplot, width=15, height=5, units="cm")


########################################################################
#   Panel 2ab. COMBINE Correlation of IC50 between virus variants      #
########################################################################

outplot <- plot_grid(Panel2a, Panel2b, align = "h", axis = "bt", rel_widths = c(4, .9)) + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
Panel2ab <- outplot
ggsave("FIGURE-Panel2ab.svg", outplot, width=15, height=5, units="cm")



########################################################################
# Panel 3 Correlation between participant age and IC50 across variants #
########################################################################

relevantData <- dose2cohort %>%
  pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50")
relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)



dataset <- c("D614G", "Kent", "SAfrica", "India2")
breakList <- c(5, 10, 64, 256, 1024, 5120)

relevantData$ic50_log2 <- log2(relevantData$ic50)
plotBreaks <- c(5, 10, 64, 256, 1024, 5120)
plotBreaks_log2 <- log2(plotBreaks)
plotBreaks_minor <- c(32, 128, 512, 2048)
plotBreaks_minor_log2 <- log2(plotBreaks_minor)


outplot<-ggscatter(relevantData, x="age", y="ic50_log2", color = "strain", alpha=0.3, shape=20, size=2,
                   add = "reg.line", 
                   add.params = list(color = "black", fill = "grey"), # Customize reg. line
                   conf.int = TRUE, # Add confidence interval
) + 
  scale_y_continuous(
    breaks=plotBreaks_log2, 
    labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
    minor_breaks=plotBreaks_minor_log2
  )+
  scale_x_continuous(
    breaks=c(40, 60) 
  ) +
  xlab("Age")+
  ylab( bquote('Virus Neutralisation, '~IC[50]~ '')) +
  stat_cor(method = "spearman", label.sep="\n", label.y.npc = 0.2, size=3) + 
  theme_bw(base_family = "Helvetica Neue Thin") +
  theme(legend.position="none")+
  scale_colour_muted() +
  facet_wrap(  ~ strain, nrow = 1) + 
  theme(      strip.background = element_blank() ) + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# outplot
Panel3 <- outplot
ggsave("FIGURE-Panel3.svg", outplot, width=15, height=5, units="cm")



###########################################################################
# Panel 3bmi Correlation between participant BMI and IC50 across variants #
###########################################################################

relevantData <- dt %>% 
  filter(COVID_vaccStatus %in% c(2)) %>%
  filter(sampleOrderInVaccStatus == 1) %>%
  pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50")
relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)



dataset <- c("D614G", "Kent", "SAfrica", "India2")
breakList <- c(5, 10, 64, 256, 1024, 5120)

relevantData$ic50_log2 <- log2(relevantData$ic50)
plotBreaks <- c(5, 10, 64, 256, 1024, 5120)
plotBreaks_log2 <- log2(plotBreaks)
plotBreaks_minor <- c(32, 128, 512, 2048)
plotBreaks_minor_log2 <- log2(plotBreaks_minor)


outplot<-ggscatter(relevantData, x="BMI", y="ic50_log2", color = "strain", alpha=0.5,
                   add = "reg.line", 
                   add.params = list(color = "black", fill = "grey"), # Customize reg. line
                   conf.int = TRUE, # Add confidence interval
) + 
  scale_y_continuous(
    breaks=plotBreaks_log2, 
    labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
    minor_breaks=plotBreaks_minor_log2
  )+
  scale_x_continuous(
    trans = "log2",
    breaks = c(16, 24, 32, 40, 48)
  ) +
  xlab("BMI")+
  ylab(bquote('Virus Neutralisation, '~IC[50]~ '')) +
  stat_cor(method = "spearman", label.sep="\n", label.y.npc = 0.2) + 
  theme_bw(base_family = "Helvetica Neue Thin") +
  theme(legend.position="none")+
  scale_colour_muted() +
  facet_wrap(  ~ strain, nrow = 1) + 
  theme(      strip.background = element_blank() ) + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

 # outplot
Panel3bmi <- outplot
ggsave("FIGURE-Panel3bmi.svg", outplot, width=15, height=5, units="cm")


###########################################################################
# Panel 3sex Correlation between participant BMI and IC50 across variants #
###########################################################################

relevantData <- dt %>% 
  filter(COVID_vaccStatus %in% c(2)) %>%
  filter(sampleOrderInVaccStatus == 1) %>%
  pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50")
relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)



dataset <- c("D614G", "Kent", "SAfrica", "India2")
breakList <- c(5, 10, 64, 256, 1024, 5120)

relevantData$ic50_log2 <- log2(relevantData$ic50)
plotBreaks <- c(5, 10, 64, 256, 1024, 5120)
plotBreaks_log2 <- log2(plotBreaks)
plotBreaks_minor <- c(32, 128, 512, 2048)
plotBreaks_minor_log2 <- log2(plotBreaks_minor)


outplot<-ggviolin(relevantData, x="sex", y="ic50_log2", color = "strain", trim=TRUE) +
  geom_jitter(shape=20, position=position_jitter(0.2), alpha=0.3) + 
  scale_y_continuous(
    breaks=plotBreaks_log2, 
    labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
    minor_breaks=plotBreaks_minor_log2
  )+
  xlab("Sex")+
  ylab(bquote('Virus Neutralisation, '~IC[50]~ '')) +
  stat_compare_means( method = "wilcox.test" )+
  theme_bw(base_family = "Helvetica Neue Thin") +
  theme(legend.position="none")+
  scale_colour_muted() +
  facet_wrap(  ~ strain, nrow = 1) + 
  theme(      strip.background = element_blank() ) + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

# outplot
Panel3sex <- outplot
ggsave("FIGURE-Panel3sex.svg", outplot, width=15, height=5, units="cm")


########################################################################
#   Panel 3sexBMI COMBINE variables with no correlation between       #
########################################################################

outplot <- plot_grid(Panel3sex, Panel3bmi, align = "h", nrow=2, labels = "AUTO")
Panel3sexBMI <- outplot
ggsave("FIGURE-Panel3sexBMI.svg", outplot, width=18, height=18, units="cm")


########################################################################
#  Panel 4  Corr. between days since dose for repeat visits           #
########################################################################


relevantData <- dt %>% 
  filter(COVID_vaccStatus %in% c(2)) %>%
  drop_na(Wuhan1_ic50) %>%
  group_by(participantAndVaccStatus) %>%
  mutate(maxVisitsInStatus = n()) %>%
  filter(maxVisitsInStatus >= 2) %>%
  mutate(latestVisit = max(collection_date)) %>% 
  mutate(earliestVisit = min(collection_date)) %>% 
  mutate(longestInterval = difftime(latestVisit, earliestVisit)) %>%
  slice(c(1,n())) %>%
  filter(longestInterval >= 51.65)  # 95%ile of 2nd jab sera collected is 51.65, so take > that

relevantData <- relevantData %>% pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50")
relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)

dataset <- c("D614G", "Kent", "SAfrica", "India2")
breakList <- c(5, 10, 64, 256, 1024, 5120)

relevantData$ic50_log2 <- log2(relevantData$ic50)
plotBreaks <- c(5, 10, 64, 256, 1024, 5120)
plotBreaks_log2 <- log2(plotBreaks)
plotBreaks_minor <- c(32, 128, 512, 2048)
plotBreaks_minor_log2 <- log2(plotBreaks_minor)


outplot<-ggscatter(relevantData, x="COVID_daysSinceJab2", y="ic50_log2", color = "strain", alpha=0.5) + 
  geom_line(aes(group = bc_participant_id))+
  scale_y_continuous(
    breaks=plotBreaks_log2, 
    labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
    minor_breaks=plotBreaks_minor_log2
  )+
  scale_x_continuous(
    limits = c(7, 125),
    breaks = c(25, 100)
  ) +
  xlab("Days since Dose 2")+
  ylab(bquote('Virus Neutralisation, '~IC[50]~ '')) +
  theme_bw(base_family = "Helvetica Neue Thin") +
  theme(legend.position="none")+
  scale_colour_muted() +
  facet_wrap(  ~ strain, nrow = 1) + 
  theme(      strip.background = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()
  )
# outplot
Panel4 <- outplot
ggsave("FIGURE-Panel4.svg", outplot, width=15, height=5, units="cm")

### Stats for Panel 4 ###
nrow(relevantData)/5

rt <- pTtest12(relevantData, "Wuhan1")
rt <- pTtest12(relevantData, "D614G", rt)
rt <- pTtest12(relevantData, "Kent", rt)
rt <- pTtest12(relevantData, "SAfrica", rt)
rt <- pTtest12(relevantData, "India2", rt)

rt <- pWtest12(relevantData, "Wuhan1", rt)
rt <- pWtest12(relevantData, "D614G", rt)
rt <- pWtest12(relevantData, "Kent", rt)
rt <- pWtest12(relevantData, "SAfrica", rt)
rt <- pWtest12(relevantData, "India2", rt)
rt




########################################################################
#   Panel 3p4 COMBINE Correlation of IC50 between virus variants      #
########################################################################

outplot <- plot_grid(Panel3, Panel4, align = "h", nrow=2)
Panel3p4 <- outplot
ggsave("FIGURE-Panel3p4.svg", outplot, width=10, height=10, units="cm")




########################################################################
#   Panel 5. Vaccine responses per strain vs. 1 dose Violinplot        #
########################################################################


relevantData <- dt %>% 
  filter(COVID_vaccStatus %in% c(1,2)) %>%
  filter(sampleOrderInVaccStatus == 1) %>%
  pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50")
relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)

relevantData$COVID_vaccStatus_pretty <- c("1", "2")[relevantData$COVID_vaccStatus]

relevantData %>% group_by(strain) %>% summarise(counts = n(), qtile = c(0.25, 0.5, 0.75),  value = quantile(ic50, c(0.25, 0.5, 0.75), na.rm=TRUE))




outplot <- ggplot(relevantData, aes(x=COVID_vaccStatus_pretty, y=ic50, color = strain, label = sample_barcode)) + 
  scale_colour_muted() +
  # scale_color_brewer(palette="Set1") + 
  geom_hline(yintercept = referenceIC50,linetype = 3 ) +
  geom_violin(trim=TRUE) + 
  # stat_summary(fun.data = give.n, geom = "text", fun.y = median,
  #              position = position_dodge(width = 0.75)) +
  # scale_y_log10(limits = c(0.1, 2^12))
  scale_y_continuous(trans='log2', 
                     breaks=c(5, 10, 64, 256, 1024, 5120), 
                     labels=c("[0]", "[<40]", "64", "256", "1024", "[>2560]"),
                     minor_breaks=c(32, 128, 512, 2048)) +
  ylab(bquote('Virus Neutralisation, '~IC[50]~ '')) +
  xlab("Vaccine Doses")+
  geom_jitter(shape=20, position=position_jitter(0.2), alpha=0.3) + 
  facet_wrap(~ strain, nrow = 1) +
  stat_summary(fun=median, geom = "point", color="black",  shape=5, size=1, stroke=1) + 
  theme_bw(base_family = "Helvetica Neue Thin") +
  theme(legend.position="none")+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_blank()
  )

# outplot

Panel5 <- outplot

ggsave("FIGURE-Panel5.svg", outplot, width=50, height=15, units="cm")

# ggplotly(outplot)

### Stats for Panel 5 ###
# grouped <- relevantData %>% group_by(strain)
# grouped %>% summarise(n())
# testTable <- relevantData %>% ungroup %>% select(sample_barcode, strain, ic50) 
# pairwise.wilcox.test(testTable$ic50, testTable$strain, paired = TRUE)
# 
# 
# if (runBootstrapStats == TRUE){
#   CIbootstrap(testTable, c("Wuhan1", "D614G"))
#   CIbootstrap(testTable, c("Wuhan1", "Kent"))
#   CIbootstrap(testTable, c("Wuhan1", "SAfrica"))
#   CIbootstrap(testTable, c("Wuhan1", "India2"))
#   
#   CIbootstrap(testTable, c("D614G", "Kent"))
#   CIbootstrap(testTable, c("SAfrica", "India2"))
# }


########################################################################
#   Panel 7. Vaccine responses  1 dose vs age Violinplot        #
########################################################################

dtToPlot <- dt 
relevantData <- dtToPlot %>% 
  filter(COVID_vaccStatus %in% c(1,2)) %>%
  filter(sampleOrderInVaccStatus == 1) %>%
  pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50") %>% drop_na(ic50) %>%
  mutate(ic50Range = cut(ic50, breaks=c(0,40,9999), right=FALSE, label=c("<40", ">40")))
relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)

relevantData$COVID_vaccStatus_pretty <- c("1 Dose", "2 Doses")[relevantData$COVID_vaccStatus]

relevantData %>% group_by(strain) %>% summarise(counts = n(), qtile = c(0.25, 0.5, 0.75),  value = quantile(ic50, c(0.25, 0.5, 0.75), na.rm=TRUE))

outplot <- ggplot(relevantData, aes(x=ic50Range, y=age, color = strain, label = sample_barcode)) + 
  scale_colour_muted() +
  # scale_color_brewer(palette="Set1") + 
  geom_violin(trim=TRUE) + 
  # stat_summary(fun.data = give.n, geom = "text", fun.y = median,
  #              position = position_dodge(width = 0.75)) +
  # scale_y_log10(limits = c(0.1, 2^12))
  scale_y_continuous(#trans='log2',
                     breaks=c(35, 50, 65) ,
                     labels=c("35", "", "65"),
                     limits = c(15, 80)) +
  xlab(bquote('Virus Neutralisation, '~IC[50]~ '')) +
  ylab("Participant Age")+
  geom_jitter(shape=20, position=position_jitter(0.2), alpha=0.05) + 
  facet_grid( COVID_vaccStatus_pretty ~ strain, switch = "y") +
  stat_summary(fun=median, geom = "point", color="black",  shape=5, size=1, stroke=1) + 
  theme_bw(base_family = "Helvetica Neue Thin") +
  theme(legend.position="none")+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside"
    
  ) + rotate()

# outplot
ggplotly(outplot)
Panel7 <- outplot

ggsave("FIGURE-Panel7.jpg", Panel7, width=18, height=6, units="cm")
ggsave("FIGURE-Panel7.svg", Panel7, width=50, height=15, units="cm")

# ggplotly(outplot)

### Stats for Panel 7 ###
grouped <- relevantData %>% group_by(strain, ic50Range) %>% filter(COVID_vaccStatus == 1)
grouped %>% summarise(n())
testTable <- relevantData %>% ungroup() %>% select(strain, ic50Range, age) %>% group_by(strain)

strainToTest <- "Kent"
low <- relevantData %>% filter(COVID_vaccStatus == 1, strain == strainToTest, ic50Range == "<40") %>% pull("age")
high <- relevantData %>% filter(COVID_vaccStatus == 1, strain == strainToTest, ic50Range == ">40") %>% pull("age")
t.test(low, high, paired = FALSE)


if (runBootstrapStats == TRUE){
  CIbootstrap(testTable, c("Wuhan1", "D614G"))
  CIbootstrap(testTable, c("Wuhan1", "Kent"))
  CIbootstrap(testTable, c("Wuhan1", "SAfrica"))
  CIbootstrap(testTable, c("Wuhan1", "India2"))

  CIbootstrap(testTable, c("D614G", "Kent"))
  CIbootstrap(testTable, c("SAfrica", "India2"))
}

########################################################################
#   Panel 6 - Neutralisation barplot across first and second doses     #
########################################################################


dtToPlot <- dt 
relevantData <- dtToPlot %>% 
  filter(COVID_vaccStatus %in% c(1,2)) %>%
  filter(sampleOrderInVaccStatus == 1) %>%
  pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50") %>% drop_na(ic50) %>%
  mutate(ic50Range = cut(ic50, breaks=c(0,40,256,9999), right=FALSE, label=c("<40", "40-256", ">256")))
relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)

relevantData$COVID_vaccStatus_pretty <- c("1 Dose", "2 Doses")[relevantData$COVID_vaccStatus]

outplot <- ggplot(relevantData, aes(x=ic50Range))  + geom_bar(aes(fill=strain, alpha=0.75)) +scale_fill_muted() + coord_flip() + 
  facet_grid( COVID_vaccStatus_pretty ~ strain, switch = "y") +
  theme_bw(base_family = "Helvetica Neue Thin") +
  ylab("Participants") +
  xlab( bquote('Virus Neutralisation, '~IC[50]~ '') ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position="none"   ,
    strip.background = element_blank(),
    strip.placement = "outside"
    )
# outplot
Panel6 <- outplot
ggsave("FIGURE-Panel6.svg", outplot, width=50, height=15, units="cm")

# ggplotly(outplot)

########################################################################
#   Panel 6stats - Ordered Logistical regression  for dose 1           #
########################################################################

strainOrder <- c("Wuhan1", "D614G", "Kent", "SAfrica", "India2")

dtToPlot <- dt
relevantData <- dose1cohort %>%
  pivot_longer(cols = ends_with("ic50"), names_to = "strain", values_to = "ic50") %>% drop_na(ic50) %>%
  mutate(ic50Range = cut(ic50, breaks=c(0,40,256,9999), right=FALSE, label=c("<40", "40-256", ">256")))
relevantData$strain <- str_replace_all(relevantData$strain, pattern = "_ic50", replacement = "")
relevantData$strain <- fct_relevel(relevantData$strain, strainOrder)

olrTable <- relevantData %>% ungroup() %>% select(bc_participant_id , age, strain, ic50Range)
nrow(olrTable)
save(olrTable, file ="olrTable-IDshuffled.rda")

# Ordered logistical regression model required
fit_first <- rms::lrm( ic50Range ~ strain , data=olrTable , x=TRUE, y=TRUE)
fit_first
# exact p-values rather than pretty table
1-pchisq(coef(fit_first)^2/diag(vcov(fit_first)),1)
my.valid <- validate(fit_first, method="boot", B=1000)
my.calib <- calibrate(fit_first, method="boot", B=1000)
par(bg="white", las=1)
plot(my.calib, las=1)


# Ordered logistical regression model WITH AGE  --- if I look at 'C' this
fit_first_age <- rms::lrm( ic50Range ~ strain*age , data=olrTable , x=TRUE, y=TRUE)
fit_first_age
anova(fit_first_age)
# exact p-values rather than pretty table
1-pchisq(coef(fit_first_age)^2/diag(vcov(fit_first_age)),1)
# boostrap to check...
my.valid <- validate(fit_first_age, method="boot", B=1000)
my.calib <- calibrate(fit_first_age, method="boot", B=1000)
par(bg="white", las=1)
plot(my.calib, las=1)


#### Now just consider age alone
fit_age <- rms::lrm( ic50Range ~ age , data=olrTable , x=TRUE, y=TRUE)
fit_age
1-pchisq(coef(fit_age)^2/diag(vcov(fit_age)),1)
my.valid <- validate(fit_age, method="boot", B=1000)
my.calib <- calibrate(fit_age, method="boot", B=1000)
par(bg="white", las=1)
plot(my.calib, las=1)



########################################################################
#                  OUTPUT MAIN COMBINED FIGURE                         #
########################################################################
# 
# col1 <- plot_grid(NULL, Panel1, Panel5, Panel6, ncol=1, rel_heights = c(0.065, 0.8, 0.6, 0.4) )
# col2 <- plot_grid(Panel2ab, Panel3, Panel4, ncol=1, rel_heights = c(1,1,1) )
# block12 <- plot_grid(col1, col2, rel_widths = c(1, 1.3), align = "h", axis='tb')
# 
# ggsave("FIGURE-Block12_1p25x.svg", block12, width=25, height=22, units="cm")
# ggsave("FIGURE-Block12_1p25x.jpg", block12, width=25, height=22, units="cm")
# 


########################################################################
#                  OUTPUT MAIN COMBINED FIGURE                         #
########################################################################

col1 <- plot_grid(NULL, Panel1, Panel5, Panel6, Panel7, ncol=1, rel_heights = c(0.065, 0.8, 0.6, 0.5, 0.5) )
col2 <- plot_grid(Panel2ab, Panel3, Panel4, ncol=1, rel_heights = c(1,1,1) )
block12 <- plot_grid(col1, col2, rel_widths = c(1, 1.3), align = "h", axis='tb')



ggsave("FIGURE-Block12_1p25x.svg", block12, width=25, height=25, units="cm")
ggsave("FIGURE-Block12_1p25x.jpg", block12, width=25, height=25, units="cm", dpi = 600)

