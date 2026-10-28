# DNAmTL analysis
#
# Purpose:
#   Calculate DNA methylation-based telomere length (DNAmTL), compare DNAmTL
#   with TelSeq-estimated leukocyte telomere length (TelSeq_LTL), and examine associations
#   with age, sex, and genetically inferred ancestry.
#
# Analysis overview:
#   1. Load packages and retrieve DNAmTL clock coefficients.
#   2. Estimate DNAmTL from methylation data and merge analysis datasets.
#   3. Create analysis subsets and assess DNAmTL/TelSeq agreement.
#   4. Generate plots and correlations overall, by sex, ancestry, and age group.
#   5. Fit regression models. 
#

# 1. Load required packages and initialize the analysis environment ----
library(methylclock)
library(methylclockData)
library(ExperimentHub)
library(tidyverse)
library(ggpubr)
library(dplyr)


# 1.1 Set up ExperimentHub cache ----
# Only needed when running locally in an environment that requires a configured ExperimentHub cache.
setExperimentHubOption("CACHE", "/path/to/ExperimentHub/")
ExperimentHub(localHub = T)

# 1.2 Retrieve CpG probes used by the DNAmTL clock ----
DNAmTL <- get_coefTL()
head(DNAmTL)
DNAmTL_probe_ids <- DNAmTL$CpGmarker[-1]


# 2.1 Load methylation data and calculate DNAmTL ----
# Load MVP methylation data
# Dataset with a subset of DNAmTL probes only
mvp_meth <- read_csv("/path/to/DNAm_data") 
head(mvp_meth)
dim(mvp_meth)

# Check overlap between available and required CpG probes
length(which(mvp_meth$CpG %in% DNAmTL_probe_ids)) #113/140 ~ 80% in MVP

# Estimate DNAmTL using the telomere-length clock and save estimates
DNAmTL <- DNAmAge(mvp_meth, clocks = "TL")
head(DNAmTL)
write.table(DNAmTL, "DNAmTL.csv", sep = ",", row.names = F, col.names = T)

# 2.2 Load DNAmTL, TelSeq, Age, Sex, and Genetic Ancestry data ----
DNAmTL_dat <- read_csv("/path/to/DNAmTL")
telseq <- read_csv("/path/to/telseq_LTL")
age_sex <- read_csv("/path/to/age_sex") 
gia_dat <- read_csv("/path/to/genetically_inferred_ancestry")

DNAmTL_telseq <- merge(DNAmTL_dat, telseq, by="participant_id")
DNAmTL_telseq_pheno <- merge(DNAmTL_telseq, age_sex, by = "participant_id")
DNAmTL_telseq_pheno <- merge(DNAmTL_telseq_pheno, gia_dat, by = "participant_id")

#DNAmTL_telseq_pheno is a dataset with the following columns: participant_id, DNAmTL, telseqLTL, age, SEX, gia.ancestry

# 3. Define analysis populations ----
# Full cohort analyses
# Includes all participants with overlapping methylation and WGS data across ancestry groups and ages. 
# Exclude participants missing sex or age
DNAmTL_telseq_pheno_all <- subset(DNAmTL_telseq_pheno, !(is.na(DNAmTL_telseq_pheno$SEX)) & 
                              !(is.na(DNAmTL_telseq_pheno$age)))

# Primary analyses
# Restricted to participants <80 years with EUR, AFR, or AMR ancestry.
DNAmTL_telseq_pheno_primary <- subset(DNAmTL_telseq_pheno_all, age < 80 & 
                                    (gia.ancestry == "EUR" |  gia.ancestry == "AFR" |
                                      gia.ancestry == "AMR"))

# 3.1 Correlation plot of DNAmTL with TelSeq in the primary analysis population (Figure 4) ----
# Primary analyses: participants <80 years with EUR, AFR, or AMR ancestry.

ggplot(DNAmTL_telseq_pheno_primary, aes(x=DNAmTL, y=telseq_LTL)) + geom_point(color="grey", alpha=0.5, size=2.5) +
  geom_smooth(method="lm", color="red", se=T)+
  xlab("DNAm-based LTL (kbp)")+ ylab("TelSeq estimated LTL (kbp)") + ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=12),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))

cor.test(DNAmTL_telseq_pheno_primary$DNAmTL, DNAmTL_telseq_pheno_primary$telseq_LTL)

# 3.2 Correlation between DNAmTL and age in the primary analysis population ----
cor.test(DNAmTL_telseq_pheno_primary$DNAmTL, DNAmTL_telseq_pheno_primary$age)

# 3.3 Assess agreement between DNAmTL and TelSeq using a Bland–Altman plot (Figure 4) ----

DNAmTL_telseq_pheno_primary <- DNAmTL_telseq_pheno_primary %>%
  mutate(
    mean_dnamtl_telseq = (DNAmTL + telseq_LTL)/2,
    diff_dnamtl_telseq = DNAmTL - telseq_LTL
  )

ggplot(DNAmTL_telseq_pheno_primary, aes(x=mean_dnamtl_telseq, y=diff_dnamtl_telseq))+
 geom_point(alpha = 0.5, color="grey", shape=19, size=2) + 
  geom_hline(yintercept = mean(DNAmTL_telseq_pheno_primary$diff_dnamtl_telseq, na.rm =T), color="blue", linetype = "dashed", size=0.6) +
  geom_hline(yintercept = mean(DNAmTL_telseq_pheno_primary$diff_dnamtl_telseq, na.rm =T) + 1.96* sd(DNAmTL_telseq_pheno_primary$diff_dnamtl_telseq, na.rm =T),
             color="red", linetype="dashed", size=0.6) +
  geom_hline(yintercept = mean(DNAmTL_telseq_pheno_primary$diff_dnamtl_telseq, na.rm =T) - 1.96* sd(DNAmTL_telseq_pheno_primary$diff_dnamtl_telseq, na.rm =T),
             color="red", linetype="dashed", size=0.6) + 
  geom_hline(yintercept = 0,color="black",size=0.5)+
  labs(x = "Mean of DNAmTL and TelSeq LTL", y= "Difference between DNAmTL and TelSeq LTL") +
  theme_bw() + xlim(3.17, 6.5) + ylim(-0.5, 6.5) +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=12),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color="black", linewidth = 0.2))

mean(DNAmTL_telseq_pheno_primary$diff_dnamtl_telseq, na.rm =T)
mean(DNAmTL_telseq_pheno_primary$diff_dnamtl_telseq, na.rm =T) + 1.96 * sd(DNAmTL_telseq_pheno_primary$diff_dnamtl_telseq, na.rm =T)
mean(DNAmTL_telseq_pheno_primary$diff_dnamtl_telseq, na.rm =T) - 1.96 * sd(DNAmTL_telseq_pheno_primary$diff_dnamtl_telseq, na.rm =T)

# 4.1 DNAmTL-age correlations by sex (Supplementary Figure 6) ----
# Use the primary analysis subset (n = 22,156)
DNAmTL_telseq_pheno_primary <- DNAmTL_telseq_pheno_primary %>% 
  mutate(point_color = case_when((SEX == "Female") ~ "darkgoldenrod",
                                (SEX == "Male") ~ "steelblue"),
         line_color = case_when((SEX == "Female") ~ "peru",
                                (SEX == "Male") ~ "steelblue4"))

DNAmTL_telseq_pheno_primary$point_color <- factor(DNAmTL_telseq_pheno_primary$point_color, levels=c("darkgoldenrod", "steelblue"))
DNAmTL_telseq_pheno_primary$line_color <- factor(DNAmTL_telseq_pheno_primary$line_color, levels=c("peru", "steelblue4"))
levels(as.factor(DNAmTL_telseq_pheno_primary$point_color))
levels(as.factor(DNAmTL_telseq_pheno_primary$line_color))

ggplot(DNAmTL_telseq_pheno_primary, aes(x=age, y=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color),method="lm", se=T)+
  ylab("DNAm-based LTL (kbp)")+ xlab("Age (y)") +  ylim(5.17,8.5) +
  scale_x_continuous(breaks=seq(0,100, by = 20), limits = c(5,105)) +
  theme_bw() +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=12),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

# Calculate the DNAmTL–age correlation among female participants
DNAmTL_telseq_pheno_primary.f<- DNAmTL_telseq_pheno_primary[DNAmTL_telseq_pheno_primary$SEX == "Female",]
nrow(DNAmTL_telseq_pheno_primary.f)
cor.test(DNAmTL_telseq_pheno_primary.f$DNAmTL, as.numeric(DNAmTL_telseq_pheno_primary.f$age))

# Calculate the DNAmTL–age correlation among male participants
DNAmTL_telseq_pheno_primary.m<- DNAmTL_telseq_pheno_primary[DNAmTL_telseq_pheno_primary$SEX == "Male",]
nrow(DNAmTL_telseq_pheno_primary.m)
cor.test(DNAmTL_telseq_pheno_primary.m$DNAmTL, as.numeric(DNAmTL_telseq_pheno_primary.m$age))

# Full cohort analyses
# Includes all eligible participants across ancestry groups and ages.
DNAmTL_telseq_pheno_all <- DNAmTL_telseq_pheno_all %>% 
  mutate(point_color = case_when((SEX == "Female") ~ "darkgoldenrod",
                                 (SEX == "Male") ~ "steelblue"),
         line_color = case_when((SEX == "Female") ~ "peru",
                                (SEX == "Male") ~ "steelblue4"))

DNAmTL_telseq_pheno_all$point_color <- factor(DNAmTL_telseq_pheno_all$point_color, levels=c("darkgoldenrod", "steelblue"))
DNAmTL_telseq_pheno_all$line_color <- factor(DNAmTL_telseq_pheno_all$line_color, levels=c("peru", "steelblue4"))
levels(as.factor(DNAmTL_telseq_pheno_all$point_color))
levels(as.factor(DNAmTL_telseq_pheno_all$line_color))

ggplot(DNAmTL_telseq_pheno_all, aes(x=age, y=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color),method="lm", se=T)+
  ylab("DNAm-based LTL (kbp)")+ xlab("Age (y)") +  ylim(5.17,8.5) +
  scale_x_continuous(breaks=seq(0,100, by = 20), limits = c(5,105)) +
  theme_bw() +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=12),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

# Calculate the DNAmTL–age correlation among female participants
DNAmTL_telseq_pheno_all.f<- DNAmTL_telseq_pheno_all[DNAmTL_telseq_pheno_all$SEX == "Female",]
nrow(DNAmTL_telseq_pheno_all.f)
cor.test(DNAmTL_telseq_pheno_all.f$DNAmTL, as.numeric(DNAmTL_telseq_pheno_all.f$age))

# Calculate the DNAmTL–age correlation among male participants
DNAmTL_telseq_pheno_all.m<- DNAmTL_telseq_pheno_all[DNAmTL_telseq_pheno_all$SEX == "Male",]
nrow(DNAmTL_telseq_pheno_all.m)
cor.test(DNAmTL_telseq_pheno_all.m$DNAmTL, as.numeric(DNAmTL_telseq_pheno_all.m$age))

# 4.2 DNAmTL-TelSeq_LTL correlations by sex (Supplementary Figure 6) ----

# Use the primary analysis subset (n = 22,156)
DNAmTL_telseq_pheno_primary <- DNAmTL_telseq_pheno_primary %>% 
  mutate(point_color = case_when((SEX == "Female") ~ "darkgoldenrod",
                                 (SEX == "Male") ~ "steelblue"),
         line_color = case_when((SEX == "Female") ~ "peru",
                                (SEX == "Male") ~ "steelblue4"))

DNAmTL_telseq_pheno_primary$point_color <- factor(DNAmTL_telseq_pheno_primary$point_color, levels=c("darkgoldenrod", "steelblue"))
DNAmTL_telseq_pheno_primary$line_color <- factor(DNAmTL_telseq_pheno_primary$line_color, levels=c("peru", "steelblue4"))
levels(as.factor(DNAmTL_telseq_pheno_primary$point_color))
levels(as.factor(DNAmTL_telseq_pheno_primary$line_color))

ggplot(DNAmTL_telseq_pheno_primary, aes(x=DNAmTL, y=telseq_LTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color),method="lm", se=T)+
  xlab("DNAm-based LTL (kbp)")+ ylab("TelSeq estimated LTL (kbp)") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=12),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

# Calculate the DNAmTL–TelSeq correlation among female participants
DNAmTL_telseq_pheno_primary.f<- DNAmTL_telseq_pheno_primary[DNAmTL_telseq_pheno_primary$SEX == "Female",]
nrow(DNAmTL_telseq_pheno_primary.f)
cor.test(DNAmTL_telseq_pheno_primary.f$DNAmTL, as.numeric(DNAmTL_telseq_pheno_primary.f$telseq_LTL))

# Calculate the DNAmTL–TelSeq correlation among male participants
DNAmTL_telseq_pheno_primary.m<- DNAmTL_telseq_pheno_primary[DNAmTL_telseq_pheno_primary$SEX == "Male",]
nrow(DNAmTL_telseq_pheno_primary.m)
cor.test(DNAmTL_telseq_pheno_primary.m$DNAmTL, as.numeric(DNAmTL_telseq_pheno_primary.m$telseq_LTL))

# Full cohort analyses
# Includes all eligible participants across ancestry groups.
DNAmTL_telseq_pheno_all <- DNAmTL_telseq_pheno_all %>% 
  mutate(point_color = case_when((SEX == "Female") ~ "darkgoldenrod",
                                 (SEX == "Male") ~ "steelblue"),
         line_color = case_when((SEX == "Female") ~ "peru",
                                (SEX == "Male") ~ "steelblue4"))

DNAmTL_telseq_pheno_all$point_color <- factor(DNAmTL_telseq_pheno_all$point_color, levels=c("darkgoldenrod", "steelblue"))
DNAmTL_telseq_pheno_all$line_color <- factor(DNAmTL_telseq_pheno_all$line_color, levels=c("peru", "steelblue4"))
levels(as.factor(DNAmTL_telseq_pheno_all$point_color))
levels(as.factor(DNAmTL_telseq_pheno_all$line_color))

ggplot(DNAmTL_telseq_pheno_all, aes(x=DNAmTL, y=telseq_LTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color),method="lm", se=T)+
  xlab("DNAm-based LTL (kbp)")+ ylab("TelSeq estimated LTL (kbp)") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=10), axis.title = element_text(size=12),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

# Calculate the DNAmTL–TelSeq correlation among female participants
DNAmTL_telseq_pheno_all.f<- DNAmTL_telseq_pheno_all[DNAmTL_telseq_pheno_all$SEX == "Female",]
nrow(DNAmTL_telseq_pheno_all.f)
cor.test(DNAmTL_telseq_pheno_all.f$DNAmTL, as.numeric(DNAmTL_telseq_pheno_all.f$telseq_LTL))

# Calculate the DNAmTL–TelSeq correlation among male participants
DNAmTL_telseq_pheno_all.m<- DNAmTL_telseq_pheno_all[DNAmTL_telseq_pheno_all$SEX == "Male",]
nrow(DNAmTL_telseq_pheno_all.m)
cor.test(DNAmTL_telseq_pheno_all.m$DNAmTL, as.numeric(DNAmTL_telseq_pheno_all.m$telseq_LTL))

# 4.3 DNAmTL-age correlations across age groups and genetic ancestry (Supplementary Figure 7) ----

# Full cohort analyses 
# Categorize participants by age 
DNAmTL_telseq_pheno_all <- DNAmTL_telseq_pheno_all %>% 
  mutate(age_categories = case_when((age <= 40) ~ 1, #20-40 y
                                    (age > 40 & age <= 60)~2,  #40-60 y
                                    (age > 60 & age <= 80)~3, #60-80 y
                                    (age > 80)~4), # >80 y
         age_categories = as.factor(age_categories))

# Create ancestry-specific analysis datasets
DNAmTL_telseq_pheno.eur <- subset(DNAmTL_telseq_pheno_all, gia.ancestry == "EUR")

DNAmTL_telseq_pheno.eur <- DNAmTL_telseq_pheno.eur %>% 
  mutate(point_color = case_when((age_categories == 1) ~ "lightsteelblue", 
                                 (age_categories == 2) ~ "steelblue",  
                                 (age_categories == 3) ~ "royalblue",
                                 (age_categories == 4) ~ "midnightblue"),
         line_color = case_when((age_categories == 1) ~ "red", 
                                (age_categories == 2) ~ "red1", 
                                (age_categories == 3) ~ "red2",
                                (age_categories == 4) ~ "red3"))

DNAmTL_telseq_pheno.afr <- subset(DNAmTL_telseq_pheno_all, gia.ancestry == "AFR")

DNAmTL_telseq_pheno.afr <- DNAmTL_telseq_pheno.afr %>% 
  mutate(point_color = case_when((age_categories == 1) ~ "thistle", 
                                 (age_categories == 2) ~ "plum",  
                                 (age_categories == 3) ~ "darkorchid",
                                 (age_categories == 4) ~ "#4B0082"),
         line_color = case_when((age_categories == 1) ~ "red", 
                                (age_categories == 2) ~ "red1", 
                                (age_categories == 3) ~ "red2",
                                (age_categories == 4) ~ "red3"))

DNAmTL_telseq_pheno.amr <- subset(DNAmTL_telseq_pheno_all, gia.ancestry == "AMR")

DNAmTL_telseq_pheno.amr <- DNAmTL_telseq_pheno.amr %>% 
  mutate(point_color = case_when((age_categories == 1) ~ "peachpuff", 
                                 (age_categories == 2) ~ "orange",  
                                 (age_categories == 3) ~ "peru",
                                 (age_categories == 4) ~ "sienna"),
         line_color = case_when((age_categories == 1) ~ "red", 
                                (age_categories == 2) ~ "red1", 
                                (age_categories == 3) ~ "red2",
                                (age_categories == 4) ~ "red3"))

DNAmTL_telseq_pheno.eas <- subset(DNAmTL_telseq_pheno_all, gia.ancestry == "EAS")

DNAmTL_telseq_pheno.eas <- DNAmTL_telseq_pheno.eas %>% 
  mutate(point_color = case_when((age_categories == 1) ~ "lightsalmon", 
                                 (age_categories == 2) ~ "tomato",  
                                 (age_categories == 3) ~ "firebrick",
                                 (age_categories == 4) ~ "maroon"),
         line_color = case_when((age_categories == 1) ~ "red", 
                                (age_categories == 2) ~ "red1", 
                                (age_categories == 3) ~ "red2",
                                (age_categories == 4) ~ "red3"))

DNAmTL_telseq_pheno.other <- subset(DNAmTL_telseq_pheno_all, gia.ancestry == "Other")

DNAmTL_telseq_pheno.other <- DNAmTL_telseq_pheno.other %>% 
  mutate(point_color = case_when((age_categories == 1) ~ "#C0C0C0", 
                                 (age_categories == 2) ~ "lightgray",  
                                 (age_categories == 3) ~ "dimgray",
                                 (age_categories == 4) ~ "black"),
         line_color = case_when((age_categories == 1) ~ "red", 
                                (age_categories == 2) ~ "red1", 
                                (age_categories == 3) ~ "red2",
                                (age_categories == 4) ~ "red3"))

DNAmTL_telseq_pheno.sas <- subset(DNAmTL_telseq_pheno_all, gia.ancestry == "SAS")

DNAmTL_telseq_pheno.sas <- DNAmTL_telseq_pheno.sas %>% 
  mutate(point_color = case_when((age_categories == 1) ~ "lightgreen", 
                                 (age_categories == 2) ~ "limegreen",  
                                 (age_categories == 3) ~ "forestgreen",
                                 (age_categories == 4) ~ "darkolivegreen"),
         line_color = case_when((age_categories == 1) ~ "red", 
                                (age_categories == 2) ~ "red1", 
                                (age_categories == 3) ~ "red2",
                                (age_categories == 4) ~ "red3"))

# Plot DNAmTL by age within each ancestry group
dnamtl_eur <- ggplot(DNAmTL_telseq_pheno.eur, aes(x=age, y=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + 
  ylab("")+ xlab("") + ylim(5.17,8.5) +
  scale_x_continuous(breaks=seq(0,100, by = 20), limits = c(5,105)) +
  theme_bw()+
  theme(axis.text = element_text(size=5), axis.title = element_text(size=12), plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

dnamtl_afr <- ggplot(DNAmTL_telseq_pheno.afr, aes(x=age, y=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("")  + ylim(5.17,8.5) +
  scale_x_continuous(breaks=seq(0,100, by = 20), limits = c(5,105)) +
  theme_bw()+
  theme(axis.text = element_text(size=5), axis.title = element_text(size=12),plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

dnamtl_amr <- ggplot(DNAmTL_telseq_pheno.amr, aes(x=age, y=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("")  + ylim(5.17,8.5) +
  scale_x_continuous(breaks=seq(0,100, by = 20), limits = c(5,105)) +
  theme_bw()+
  theme(axis.text = element_text(size=5), axis.title = element_text(size=12),plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

dnamtl_eas <- ggplot(DNAmTL_telseq_pheno.eas, aes(x=age, y=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("")  + ylim(5.17,8.5) +
  scale_x_continuous(breaks=seq(0,100, by = 20), limits = c(5,105)) +
  theme_bw()+
  theme(axis.text = element_text(size=5), axis.title = element_text(size=12),plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

dnamtl_sas <- ggplot(DNAmTL_telseq_pheno.sas, aes(x=age, y=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") + ylim(5.17,8.5) +
  scale_x_continuous(breaks=seq(0,100, by = 20), limits = c(5,105)) +
  theme_bw()+
  theme(axis.text = element_text(size=5), axis.title = element_text(size=12),plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

dnamtl_other <- ggplot(DNAmTL_telseq_pheno.other, aes(x=age, y=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") + ylim(5.17,8.5) +
  scale_x_continuous(breaks=seq(0,100, by = 20), limits = c(5,105)) +
  theme_bw()+
  theme(axis.text = element_text(size=5), axis.title = element_text(size=12),plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

ggarrange(dnamtl_eur, dnamtl_afr, dnamtl_amr, dnamtl_eas, dnamtl_sas, dnamtl_other,
          ncol=1, nrow=6)

# Calculate DNAmTL–age correlations within age groups
# DNAmTL–age association among participants with EUR ancestry
summary(as.factor(DNAmTL_telseq_pheno.eur$age_categories))
correlation_df_eur <- DNAmTL_telseq_pheno.eur %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, age, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)
print(correlation_df_eur)

# DNAmTL–age association among participants with AFR ancestry
summary(DNAmTL_telseq_pheno.afr$age_categories)
correlation_df_afr <- DNAmTL_telseq_pheno.afr %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, age, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_afr)

# DNAmTL–age association among participants with AMR ancestry
summary(DNAmTL_telseq_pheno.amr$age_categories)
correlation_df_amr <- DNAmTL_telseq_pheno.amr %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, age, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_amr)

# DNAmTL–age association among participants with EAS ancestry
summary(DNAmTL_telseq_pheno.eas$age_categories)
correlation_df_eas <- DNAmTL_telseq_pheno.eas %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, age, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_eas)

# DNAmTL–age association among participants with SAS ancestry
summary(DNAmTL_telseq_pheno.sas$age_categories)
correlation_df_sas <- DNAmTL_telseq_pheno.sas %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, age, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_sas)

# DNAmTL–age association among participants with other ancestry
summary(DNAmTL_telseq_pheno.other$age_categories)
correlation_df_other <- DNAmTL_telseq_pheno.other %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, age, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_other)

# Primary analysis subset
# Categorize participants by age
DNAmTL_telseq_pheno_primary <- DNAmTL_telseq_pheno_primary %>% 
  mutate(age_categories = case_when((age <= 40) ~ 1, #20-40 y
                                    (age > 40 & age <= 60)~2,  #40-60 y
                                    (age > 60 & age <= 80)~3, #60-80 y
                                    (age > 80)~4), # >80 y
         age_categories = as.factor(age_categories))


# Create ancestry-specific analysis datasets
DNAmTL_telseq_pheno_primary.eur <- subset(DNAmTL_telseq_pheno_primary, gia.ancestry == "EUR")

DNAmTL_telseq_pheno_primary.eur <- DNAmTL_telseq_pheno_primary.eur %>% 
  mutate(point_color = case_when((age_categories == 1) ~ "lightsteelblue", 
                                 (age_categories == 2) ~ "steelblue",  
                                 (age_categories == 3) ~ "royalblue",
                                 (age_categories == 4) ~ "midnightblue"),
         line_color = case_when((age_categories == 1) ~ "red", 
                                (age_categories == 2) ~ "red1", 
                                (age_categories == 3) ~ "red2",
                                (age_categories == 4) ~ "red3"))

DNAmTL_telseq_pheno_primary.afr <- subset(DNAmTL_telseq_pheno_primary, gia.ancestry == "AFR")

DNAmTL_telseq_pheno_primary.afr <- DNAmTL_telseq_pheno_primary.afr %>% 
  mutate(point_color = case_when((age_categories == 1) ~ "thistle", 
                                 (age_categories == 2) ~ "plum",  
                                 (age_categories == 3) ~ "darkorchid",
                                 (age_categories == 4) ~ "#4B0082"),
         line_color = case_when((age_categories == 1) ~ "red", 
                                (age_categories == 2) ~ "red1", 
                                (age_categories == 3) ~ "red2",
                                (age_categories == 4) ~ "red3"))

DNAmTL_telseq_pheno_primary.amr <- subset(DNAmTL_telseq_pheno_primary, gia.ancestry == "AMR")

DNAmTL_telseq_pheno_primary.amr <- DNAmTL_telseq_pheno_primary.amr %>% 
  mutate(point_color = case_when((age_categories == 1) ~ "peachpuff", 
                                 (age_categories == 2) ~ "orange",  
                                 (age_categories == 3) ~ "peru",
                                 (age_categories == 4) ~ "sienna"),
         line_color = case_when((age_categories == 1) ~ "red", 
                                (age_categories == 2) ~ "red1", 
                                (age_categories == 3) ~ "red2",
                                (age_categories == 4) ~ "red3"))

# Plot DNAmTL by age within each ancestry group
dnamtl_eur_subset <- ggplot(DNAmTL_telseq_pheno_primary.eur, aes(x=age, y=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + 
  ylab("")+ xlab("") + ylim(5.17,8.5) +
  scale_x_continuous(breaks=seq(0,100, by = 20), limits = c(5,105)) +
  theme_bw()+
  theme(axis.text = element_text(size=5), axis.title = element_text(size=12), plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

dnamtl_afr_subset <- ggplot(DNAmTL_telseq_pheno_primary.afr, aes(x=age, y=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("")  + ylim(5.17,8.5) +
  scale_x_continuous(breaks=seq(0,100, by = 20), limits = c(5,105)) +
  theme_bw()+
  theme(axis.text = element_text(size=5), axis.title = element_text(size=12),plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

dnamtl_amr_subset <- ggplot(DNAmTL_telseq_pheno_primary.amr, aes(x=age, y=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("")  + ylim(5.17,8.5) +
  scale_x_continuous(breaks=seq(0,100, by = 20), limits = c(5,105)) +
  theme_bw()+
  theme(axis.text = element_text(size=5), axis.title = element_text(size=12),plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

ggarrange(dnamtl_eur_subset, dnamtl_afr_subset, dnamtl_amr_subset, 
                              ncol=1, nrow=3)

# Calculate DNAmTL–age correlations within age groups
# DNAmTL–age association among participants with EUR ancestry
summary(as.factor(DNAmTL_telseq_pheno_primary.eur$age_categories))
correlation_df_eur_subset<- DNAmTL_telseq_pheno_primary.eur %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, age, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)
print(correlation_df_eur_subset)

# DNAmTL–age association among participants with AFR ancestry
summary(DNAmTL_telseq_pheno_primary.afr$age_categories)
correlation_df_afr <- DNAmTL_telseq_pheno_primary.afr %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, age, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_afr)

# DNAmTL–age association among participants with AMR ancestry
summary(DNAmTL_telseq_pheno_primary.amr$age_categories)
correlation_df_amr <- DNAmTL_telseq_pheno_primary.amr %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, age, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_amr)

# 4.4 DNAmTL-TelSeq correlations across age groups and genetic ancestry (Supplementary Figure 8) ----

#Full cohort analyses
# DNAmTL–TelSeq association among participants with EUR ancestry
#<40y 
p_age1.eur <- ggplot(subset(DNAmTL_telseq_pheno.eur,age_categories==1), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
 ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

#40-60y 
p_age2.eur <-ggplot(subset(DNAmTL_telseq_pheno.eur,age_categories==2), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=6, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

#60-80y 
p_age3.eur <-ggplot(subset(DNAmTL_telseq_pheno.eur,age_categories==3), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
 ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

#>80y
p_age4.eur <-ggplot(subset(DNAmTL_telseq_pheno.eur,age_categories==4), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

# Calculate DNAmTL–TelSeq correlations within age groups
summary(DNAmTL_telseq_pheno.eur$age_categories)
correlation_df_eur <- DNAmTL_telseq_pheno.eur %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, telseq_LTL, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_eur)

# DNAmTL–TelSeq association among participants with AFR ancestry
p_age1.afr <- ggplot(subset(DNAmTL_telseq_pheno.afr,age_categories==1), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
 ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age2.afr <-ggplot(subset(DNAmTL_telseq_pheno.afr,age_categories==2), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age3.afr <-ggplot(subset(DNAmTL_telseq_pheno.afr,age_categories==3), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age4.afr <-ggplot(subset(DNAmTL_telseq_pheno.afr,age_categories==4), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

# Calculate DNAmTL–TelSeq correlations within age groups
summary(DNAmTL_telseq_pheno.afr$age_categories)
correlation_df_afr <- DNAmTL_telseq_pheno.afr %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, telseq_LTL, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_afr)

# DNAmTL–TelSeq association among participants with AMR ancestry
p_age1.amr <- ggplot(subset(DNAmTL_telseq_pheno.amr,age_categories==1), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age2.amr <-ggplot(subset(DNAmTL_telseq_pheno.amr,age_categories==2), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age3.amr <-ggplot(subset(DNAmTL_telseq_pheno.amr,age_categories==3), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age4.amr <-ggplot(subset(DNAmTL_telseq_pheno.amr,age_categories==4), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

# Calculate DNAmTL–TelSeq correlations within age groups
summary(DNAmTL_telseq_pheno.amr$age_categories)
correlation_df_amr <- DNAmTL_telseq_pheno.amr %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, telseq_LTL, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_amr)

# DNAmTL–TelSeq association among participants with EAS ancestry
p_age1.eas <- ggplot(subset(DNAmTL_telseq_pheno.eas,age_categories==1), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age2.eas <-ggplot(subset(DNAmTL_telseq_pheno.eas,age_categories==2), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age3.eas <-ggplot(subset(DNAmTL_telseq_pheno.eas,age_categories==3), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age4.eas <-ggplot(subset(DNAmTL_telseq_pheno.eas,age_categories==4), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

# Calculate DNAmTL–TelSeq correlations within age groups
summary(DNAmTL_telseq_pheno.eas$age_categories)
correlation_df_eas <- DNAmTL_telseq_pheno.eas %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, telseq_LTL, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_eas)

# DNAmTL–TelSeq association among participants with SAS ancestry
p_age1.sas <- ggplot(subset(DNAmTL_telseq_pheno.sas,age_categories==1), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age2.sas <-ggplot(subset(DNAmTL_telseq_pheno.sas,age_categories==2), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age3.sas <-ggplot(subset(DNAmTL_telseq_pheno.sas,age_categories==3), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age4.sas <-ggplot(subset(DNAmTL_telseq_pheno.sas,age_categories==4), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

# Calculate DNAmTL–TelSeq correlations within age groups
summary(DNAmTL_telseq_pheno.sas$age_categories)
correlation_df_sas <- DNAmTL_telseq_pheno.sas %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, telseq_LTL, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_sas)

# DNAmTL–TelSeq association among participants with other ancestry
p_age1.other <- ggplot(subset(DNAmTL_telseq_pheno.other,age_categories==1), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age2.other <-ggplot(subset(DNAmTL_telseq_pheno.other,age_categories==2), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age3.other <-ggplot(subset(DNAmTL_telseq_pheno.other,age_categories==3), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

p_age4.other <-ggplot(subset(DNAmTL_telseq_pheno.other,age_categories==4), aes(y=telseq_LTL, x=DNAmTL)) + 
  geom_point(aes(color=point_color), alpha=0.5, size=1.5) +
  geom_smooth(aes(color=line_color), method="lm", formula = y ~ x, se=T) + theme_minimal()+
  ylab("")+ xlab("") +  ylim(-1.1,5)+ xlim(5.17,8.5) +
  theme_bw() +
  theme(axis.text = element_text(size=6), axis.title = element_text(size=12),
        plot.title=element_text(size=12, hjust=0.5),
        plot.margin = margin(0,0,0,0, "cm"),
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "black", linewidth=0.2))+
  scale_color_identity()

# Calculate DNAmTL–TelSeq correlations within age groups
summary(DNAmTL_telseq_pheno.other$age_categories)
correlation_df_other <- DNAmTL_telseq_pheno.other %>%
  group_by(age_categories) %>%
  summarise(cor_test = list(cor.test(DNAmTL, telseq_LTL, use = "complete.obs")), 
            .groups = "drop") %>%
  mutate(corr= sapply(cor_test, function(x) x$estimate), 
         p_val = sapply(cor_test, function(x) x$p.value)) %>%
  select(age_categories, corr, p_val)

print(correlation_df_other)

ggarrange(p_age1.eur, p_age2.eur, p_age3.eur, p_age4.eur,
                      p_age1.afr, p_age2.afr, p_age3.afr, p_age4.afr,
                      p_age1.amr, p_age2.amr, p_age3.amr, p_age4.amr,
                      p_age1.eas, p_age2.eas, p_age3.eas, p_age4.eas ,
                      p_age1.sas, p_age2.sas, p_age3.sas, p_age4.sas, 
                      p_age1.other, p_age2.other, p_age3.other, p_age4.other,
                      ncol=4, nrow=6, widths = c(0.1,0.1))

# 5.1 Model DNAmTL as a function of age, sex, and genetic ancestry (Table 2, Model 2) ----
DNAmTL_telseq_pheno_primary$gia.ancestry <- as.factor(DNAmTL_telseq_pheno_primary$gia.ancestry)
DNAmTL_telseq_pheno_primary$gia.ancestry <- relevel(DNAmTL_telseq_pheno_primary$gia.ancestry, ref = "EUR") # Set EUR as reference
levels(DNAmTL_telseq_pheno_primary$gia.ancestry)

summary(lm(DNAmTL ~ age + SEX + gia.ancestry, data=DNAmTL_telseq_pheno_primary))

# 5.2 Model age as a function of DNAmTL, TelSeq, sex, and genetic ancestry (Table 2, Model 3) ----
summary(lm(age ~ DNAmTL +telseq_LTL+ SEX + gia.ancestry, data=DNAmTL_telseq_pheno_primary))
