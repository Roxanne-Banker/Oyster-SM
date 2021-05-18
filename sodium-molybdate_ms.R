#rm(list = ls())

setwd("/Users/roxannebanker/bioinformatics/projects/oyster-oa/microbiome/")

# prepping for GIT

# February 2020
# Document to analyze oyster larvae/juvenile data from 
# Oyster multi-omics project run in Aug-Oct 2019
# Following this tutorial: https://benjjneb.github.io/dada2/tutorial.html
# and the Harbor_R_Notebook from Cassie Ettinger

# loading in required packages

library(phyloseq)
#packageVersion("phyloseq")
library(ggplot2)
#packageVersion("ggplot2")
library(rmarkdown)
library(dada2)
#packageVersion("dada2")
library(decontam)
#packageVersion("decontam")
library(phangorn)
library(DECIPHER)
library(vegan)
#packageVersion("vegan")
library(RColorBrewer)
library(reshape)
library(coin)
library(FSA)
#packageVersion("FSA")
library(Biostrings)
library(DESeq2)
library(decontam)
library(gridExtra)
library(seqinr)
library(viridis)
library(gtable)
library(grid)
library(ALDEx2)
library(dplyr)

set.seed(5311)


# Summary Function defined for later use
#################################################################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  
  return(datac)
}
#################################################################



# Alpha diversity plots
#################################################################
# reading in the phyloseq object that includes the whole dataset as a phyloseq object
data <- readRDS("data.rds")
sample.names <- readRDS("sample_names.rds")
# reading in mapping file data
# reading in mapping file
mapping <- read.csv("phyloseq-metadata.csv")
row.names(mapping) <- mapping$sample_names

# creating a list to subset phyloseq object to remove negative and positive controls
names <- sample_names(data)
sample_names(data)
# list of the samples we want to keep
names <- names[-c(1,2,12:26,29,35,36)]

# prune samples based on that list
ps.alpha <- prune_samples(names,data)

# check to make sure these are the ones we want
sample_names(ps.alpha)

# apply estimate richness function from phyloseq as observed asvs and using the shannon index
alpha.measure <- estimate_richness(ps.alpha, measures=c("Observed", "Shannon"))
# add new key column with sample names
alpha.measure$sample_name <- names 

# transform sample_names column from mapping to character class
mapping$sample_names <- as.character(mapping$sample_names)
# join these two dataframes so now we have observed and shannon index along with appropriate key columns for easy graphing
alpha.measure.final <- full_join(alpha.measure, mapping, by = c("sample_name" = "sample_names"))
# removes samples from the final dataframe that we are not dealing with here
alpha.measure.final <- na.omit(alpha.measure.final)


# treatment level summary
observed.treat <- summarySE(alpha.measure.final, measurevar="Observed", groupvars=c("treatment"))
shannon.treat <- summarySE(alpha.measure.final, measurevar="Shannon", groupvars=c("treatment"))


# don't think I will need this, can delete later
# factor ordering so the x axis is in the order i want
#observed.treat$treatment <- factor(observed.treat$treatment,levels = c("larvae", "control", "low-ph", "molybdate"))

alpha.cols <- c("#721F81FF", "#B63679FF", "#F1605DFF", "#FEAF77FF", "#FCFDBFFF")

#"#721F81FF", "#B63679FF", "#F1605DFF", "#FEAF77FF", "#FCFDBFFF"

observed.treat.plot <-
  ggplot(observed.treat, aes(x=treatment, y=Observed, fill=treatment))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se),
                width=.1, size=0.4,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "none") +
  scale_fill_manual(values=alpha.cols[3:4]) +
  scale_x_discrete(breaks=c("control","molybdate"),
                   labels=c("Control","Sodium \nMolybdate")) +
  xlab("") + ylab("Observed ASVs")

shannon.treat.plot <-
  ggplot(shannon.treat, aes(x=treatment, y=Shannon, fill=treatment))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=Shannon-se, ymax=Shannon+se),
                width=.1, size=0.4,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), legend.position = "none") +
  scale_fill_manual(values=alpha.cols[3:4]) +
  scale_x_discrete(breaks=c("control","molybdate"),
                   labels=c("Control","Sodium \nMolybdate")) +
  xlab("") + ylab("Shannon Index")



# bucket level summary
observed.bucket <- summarySE(alpha.measure.final, measurevar="Observed", groupvars=c("bucket.id","treatment"))
shannon.bucket <- summarySE(alpha.measure.final, measurevar="Shannon", groupvars=c("bucket.id","treatment"))


alpha.breaks <- as.character(observed.bucket$bucket.id)
alpha.ids <- c("C38-1",  "C38-2",  "C51", 
               "SM38-1", "SM38-2", "SM51")

# factor ordering so the x axis is in the order i want
observed.bucket$bucket.id <- factor(observed.bucket$bucket.id,levels = c("C1", "C2", "C3", "SM1", "SM2", "SM3"))


# bucket level graphs
observed.bucket.plot <- 
  ggplot(observed.bucket, aes(x=bucket.id, y=Observed, fill=treatment))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se),
                width=.1, size=0.4,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values=alpha.cols[3:4]) +
  scale_x_discrete(breaks=c("C1", "C2", "C3", "SM1", "SM2", "SM3"),
                   labels=c("C-1", "C-2", "C-3", "SM-1", "SM-2", "SM-3")) +
  xlab("") + ylab("Observed ASVs")
  
  

# factor ordering so the x axis is in the order i want
shannon.bucket$bucket.id <- factor(shannon.bucket$bucket.id,levels = c("C1", "C2", "C3","SM1", "SM2", "SM3"))

shannon.bucket.plot <- 
  ggplot(shannon.bucket, aes(x=bucket.id, y=Shannon, fill=treatment))+
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin=Shannon-se, ymax=Shannon+se),
                width=.1, size=0.4,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values=alpha.cols[3:4], labels=c("Control"," Sodium \nMolybdate")) +
  scale_x_discrete(breaks=c("C1", "C2", "C3", "SM1", "SM2", "SM3"),
                   labels=c("C-1", "C-2", "C-3", "SM-1", "SM-2", "SM-3")) +
  xlab("") + ylab("Shannon Index")


# get onto one image, both plots same size
g1a <- ggplotGrob(observed.treat.plot)
g1b <- ggplotGrob(shannon.treat.plot)
g1 <- rbind(g1a, g1b, size = "first")

g2a <- ggplotGrob(observed.bucket.plot) # if this function does not work at first, reload the ggplot2 package then run.
g2b <- ggplotGrob(shannon.bucket.plot)
g2 <- rbind(g2a, g2b, size = "first")

g3 <- cbind(g1, g2)

#g3$widths <- unit.pmax(g1$widths, g2$widths)
grid.newpage()
grid.draw(g3)




# Alpha Diversity Statistics

# bucket ID
kruskal.test(Shannon ~ bucket.id, data=alpha.measure.final)
#kruskal_test(Shannon ~ bucket.id, distribution = approximate(nresample = 9999), data=alpha.measure.final) # this call isnot working, back to the original that actually works
# Kruskal-Wallis chi-squared = 10.978, df = 5, p-value = 0.05182

# shannon.bucket.dt <- dunnTest(Shannon ~ bucket.id, data = alpha.measure.full, method = "bonferroni")
# shannon.bucket.dt

kruskal.test(Observed ~ bucket.id, data=alpha.measure.final)
# Kruskal-Wallis chi-squared = 8.3897, df = 5, p-value = 0.136

# observed.bucket.dt <- dunnTest(Observed ~ bucket.id, data = alpha.measure.full, method = "bonferroni")
# observed.bucket.dt

# nothing came back significant here, so no need for post hocs


# treatment
kruskal.test(Shannon ~ treatment, data=alpha.measure.final)
# Kruskal-Wallis chi-squared = 0.33894, df = 1, p-value = 0.5604

kruskal.test(Observed ~ treatment, data=alpha.measure.final)
# Kruskal-Wallis chi-squared = 0.22689, df = 1, p-value = 0.6338




#################################################################



# Normalizing for sample depth via variance stabalizing transformation
#################################################################
# Will be using Deseq2 to variance stabilizing transformation, instead of the classic rarefaction or
# turning counts into proportions. This is recommended by the McMurdie and Holmes 2014 Plos Computational Biology paper, Waste not Want not.
# using this tutorial: https://astrobiomike.github.io/amplicon/dada2_workflow_ex#analysis-in-r

# first we need to make a DESeq2 object
# the design I set, bucket.id, will create groups based on both treatment and time (point sampled)
data_deseq <- phyloseq_to_deseq2(data, ~ bucket.id)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(data_deseq), 1, gm_mean)
data_deseq = estimateSizeFactors(data_deseq, geoMeans = geoMeans, type="poscount")
#data_deseq = DESeq(data_deseq, test="Wald", fitType="parametric")
# great this worked!

data_deseq_vst <- varianceStabilizingTransformation(data_deseq, blind = FALSE, fitType = "parametric")

# and here is pulling out our transformed table
vst_transform <- assay(data_deseq_vst)

# okay. the VST log-like transformation produces negative values
# for counts that are <1
# therefore I am going to set these to zero so that distance metrics can be applied
# https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#negative-numbers-in-my-transformed-data-table

vst_transform[vst_transform < 0.0] <- 0.0

# and calculating our Euclidean distance matrix
euc_dist <- assay(data_deseq_vst)
euc_dist <- dist(t(vst_transform))

#################################################################


# New Phyloseq object with adjusted ASV counts via DESEQ2
#################################################################
# reading in the vst transformed sequence counts so I do not need to do that step again (takes time)
vst_transform <- readRDS("vst_transform_table.rds")

# making our phyloseq object with transformed table
vst_transform_phy <- otu_table(vst_transform, taxa_are_rows=T)

# creating phyloseq object
ps.vst <- phyloseq(vst_transform_phy, mapping, taxa)

# Calculate RA per sample
ps.vst.RA <- transform_sample_counts(ps.vst, function(x) 100 * x/sum(x))

#################################################################


# Ordinatation
#################################################################

# read in the data here so I don't have to go through the long vst steps
ps.all <- readRDS("ps_all.rds")

ps.control <- subset_samples(ps.all, sample_data(ps.all)$treatment == "control")
ps.molybdate <- subset_samples(ps.all, sample_data(ps.all)$treatment == "molybdate")

ps.final <- merge_phyloseq(ps.control,ps.molybdate)

# creating a factor so I can use a discrete scale later for figure
sample_data(ps.final)$age <- as.factor(sample_data(ps.final)$age)

# Ordination between control and molybdate samples
ord.nmds.bray <- ordinate(ps.final, method="NMDS", distance="bray")
ord.nmds.bray


bucket.labs <- c("C-1", "C-2", "C-3", "SM-1", "SM-2", "SM-3")

bucket.shapes <-c(19, 19, 17,
                  19, 19, 17)

plot_ordination(ps.final, ord.nmds.bray, shape="bucket.id", color="bucket.id") +
  geom_point(size = 4) +
  theme(legend.title=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_colour_manual(name = "Bucket ID",
                      labels = bucket.labs,
                      values = magma(9)[3:8]) +   
  scale_shape_manual(name = "Bucket ID",
                     labels = bucket.labs,
                     values = bucket.shapes) +
  annotate(geom="text", x=.33, y=-0.25, label="Stress = 0.089", color="black")
  


# Ordination Statistics
# pulling out sample data for later
sample.data <- data.frame(sample_data(ps.final))

# Control - Molybate - Treatment 

# distance call
Dist.cntrl.sm = phyloseq::distance(ps.final, method = "bray", type="samples")

# betadisper: dispersion test
cntrl.sm.beta <- betadisper(Dist.cntrl.sm, sample.data$treatment)
permutest(cntrl.sm.beta)
#           Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.000007 0.0000068 0.002    999  0.974
# Residuals 14 0.046494 0.0033210    

# Adonis test
adonis(Dist.cntrl.sm ~ treatment, data = sample.data)
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
# treatment  1   0.22794 0.227944  3.1197 0.18223  0.008 **
# Residuals 14   1.02291 0.073065         0.81777          
# Total     15   1.25086                  1.00000 



# Control - Molybdate - Treatment/Time

# betadisper: dispersion test
cntrl.sm.beta.group <- betadisper(Dist.cntrl.sm, sample.data$group)
permutest(cntrl.sm.beta.group)
#            Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     3 0.0061976 0.0020659 1.2804    999   0.37
# Residuals 12 0.0193617 0.0016135    

# Adonis test
adonis(Dist.cntrl.sm ~ group, data = sample.data)
#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# group      3   0.78909 0.263029  6.8354 0.63084  0.001 ***
# Residuals 12   0.46177 0.038481         0.36916           
# Total     15   1.25086                  1.00000 


# Pairwise tests between treatment-time points (groups column in sample.data)

ps.c38 <- subset_samples(ps.all, group=="C38")
ps.c51 <- subset_samples(ps.all, group=="C51")
ps.sm38 <- subset_samples(ps.all, group=="SM38")
ps.sm51 <- subset_samples(ps.all, group=="SM51")


# C38:C51
ps.c38.c51 <- merge_phyloseq(ps.c38,ps.c51)
# make a data frame from the sample_data
data.c38.c51 <- data.frame(sample_data(ps.c38.c51))
Dist.c38.c51 = phyloseq::distance(ps.c38.c51, method = "bray", type="samples")

# betadisper: dispersion test
c38.c51.beta <- betadisper(Dist.c38.c51, data.c38.c51$group)
permutest(c38.c51.beta)
#            Df    Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.0012086 0.0012086 0.492    999  0.497
# Residuals  7 0.0171944 0.0024563             

# Adonis test
adonis(Dist.c38.c51 ~ group, data = data.c38.c51)
#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# group      1   0.27857 0.278570  6.4619 0.48002  0.018 *
# Residuals  7   0.30177 0.043109         0.51998         
# Total      8   0.58034                  1.00000  


# C38:SM38

ps.c38.sm38 <- merge_phyloseq(ps.c38,ps.sm38)
# make a data frame from the sample_data
data.c38.sm38 <- data.frame(sample_data(ps.c38.sm38))
Dist.c38.sm38 = phyloseq::distance(ps.c38.sm38, method = "bray", type="samples")

# betadisper: dispersion test
c38.sm38.beta <- betadisper(Dist.c38.sm38, data.c38.sm38$group)
permutest(c38.sm38.beta)
#            Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.0015648 0.00156480 4.0024    999  0.063 .
# Residuals  8 0.0031277 0.00039097

# Adonis test
adonis(Dist.c38.sm38 ~ group, data = data.c38.sm38)
#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
# group      1   0.25207 0.252071  6.2971 0.44045  0.008 **
# Residuals  8   0.32024 0.040029         0.55955          
# Total      9   0.57231                  1.00000   


# C51:SM51

ps.c51.sm51 <- merge_phyloseq(ps.c51,ps.sm51)
# make a data frame from the sample_data
data.c51.sm51 <- data.frame(sample_data(ps.c51.sm51))
Dist.c51.sm51 = phyloseq::distance(ps.c51.sm51, method = "bray", type="samples")

# betadisper: dispersion test
c51.sm51.beta <- betadisper(Dist.c51.sm51, data.c51.sm51$group)
permutest(c51.sm51.beta)
#           Df    Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.0013718 0.0013718 0.338    719 0.7014
# Residuals  4 0.0162339 0.0040585

# Adonis test
adonis(Dist.c51.sm51 ~ group, data = data.c51.sm51)
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# group      1   0.12766 0.127661   3.608 0.47423    0.1
# Residuals  4   0.14153 0.035383         0.52577       
# Total      5   0.26919                  1.00000 


# SM38:SM51

ps.sm38.sm51 <- merge_phyloseq(ps.sm38,ps.sm51)
# make a data frame from the sample_data
data.sm38.sm51 <- data.frame(sample_data(ps.sm38.sm51))
Dist.sm38.sm51 = phyloseq::distance(ps.sm38.sm51, method = "bray", type="samples")

# betadisper: dispersion test
sm38.sm51.beta <- betadisper(Dist.sm38.sm51, data.sm38.sm51$group)
permutest(sm38.sm51.beta)
#            Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.0014706 0.00147057 3.3928    999  0.138
# Residuals  5 0.0021672 0.00043344     

# Adonis test
adonis(Dist.sm38.sm51 ~ group, data = data.sm38.sm51)
#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# group      1   0.28257 0.282573  8.8303 0.63847  0.028 *
# Residuals  5   0.16000 0.032001         0.36153         
# Total      6   0.44258                  1.00000    

#################################################################



# Taxonomy tests

# All Treatments
# Comparison on means, KW, Dunn 
#################################################################

#define standard error function
se <- function(x) sqrt(var(x)/length(x))

#collapse ASVs at family level
AvgRA_o = tax_glom(ps.final, taxrank="Order", NArm = FALSE)

# comparison between treatments
#filter out taxa that don't vary > 0.2 %
AvgRA99_O = filter_taxa(AvgRA_o, function(x) var(x) > .002, TRUE)
df_o <- psmelt(AvgRA99_O)
#group and calculate mean, sd and se for different taxonomic levels
grouped_o <- group_by(df_o, treatment, Phylum, Class, Order)
avgs_o <- dplyr::summarise(grouped_o, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))

avgs_filt <- filter(avgs_o, mean > 0)
# I want to read out this table for sup data
write.csv(avgs_filt, "SM_treatment_fam_avg_sd_se.csv")


avgs_filt$ST <- as.factor(avgs_filt$treatment)
avgs_filt_st <- melt(data.frame(avgs_filt), id.vars=c("Order", "treatment"),
                     measure.vars=c("mean"))

# What is this stuff? Junk code?
# kruskal_test(mean ~ ST, distribution=approximate(nresample=9999), data=avgs_filt)
# 
# for (i in 1:length(levels(df_o$Family))) {
#   new_df <- subset(df_o, df_o$Family == levels(df_o$Family)[i])
#   new_df$ST <- as.factor(new_df$treatment)
#   print(new_df$Genus[1])
#   print(kruskal_test(Abundance ~ ST, distribution=approximate(nresample=9999), data=new_df))
#   print( dunnTest(Abundance ~ ST, data=new_df, method ="bonferroni"))}

#start
chisq = NULL
pvals = NULL
listofcats_sig = NULL
DataSet = df_o
taxa = as.factor(unique(df_o$Order))


for (cat in taxa) {
  new_df <- subset(DataSet, DataSet$Order == cat)
  new_df$ST <- as.factor(new_df$treatment)
  kw = kruskal.test(Abundance ~ ST, data=new_df)
  chisq = c(kw$statistic, chisq)
  pvals = c(kw$p.value, pvals)
  if (kw$p.value <= 0.05) {
    listofcats_sig = c(cat, listofcats_sig)
  }
}

pvals.bonf = p.adjust(pvals, method="bonferroni")
df_taxa = data.frame(rev(taxa), chisq, pvals, pvals.bonf)
write.table(df_taxa, 'SM_Treatment_VST_Mean_FAM_KW.txt', sep="\t")



# Got the Dunn to work, but don't need because I am not doing multiple comparisons between treatment groups
# just Control vs Sodium Molybdate

# pvals.dunn = NULL
# pvals.dunn.bonf = NULL
# Zsc.dunn = NULL
# comparison.dunn = NULL
# cats.dunn = NULL
# 
# for (cat in listofcats_sig){
#   new_df <- subset(DataSet, DataSet$Family == cat)
#   new_df$ST <- as.factor(new_df$treatment)
#   dT <- dunnTest(Abundance ~ ST, data = new_df, method = "bonferroni")
#   for (i in 1:length(dT$res$Comparison)) {
#     print(cat)
#     pvals.dunn <- c(dT$res$P.unadj[i], pvals.dunn)
#     pvals.dunn.bonf <- c(dT$res$P.adj[i], pvals.dunn.bonf)
#     Zsc.dunn <- c(dT$res$Z[i], Zsc.dunn)
#     comparison.dunn <- c(dT$res$Comparison, comparison.dunn)
#     cats.dunn <- c(cat, cats.dunn)
#   }
# }
# 
# df_taxa = data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bonf)
# tax_oi <- filter(df_taxa, pvals.dunn.bonf < 0.05)
# 
# write.csv(df_taxa, 'SM_Treatment_VST_Mean_FAM_KW_Dunn_OI.csv')


#################################################################


# Treatments/Age
# Comparison on means, KW, Dunn 
#################################################################


#collapse ASVs at family level
AvgRA_o = tax_glom(ps.final, taxrank="Order", NArm = FALSE)

#filter out taxa that don't vary > 0.2 %
AvgRA99_O = filter_taxa(AvgRA_o, function(x) var(x) > .002, TRUE)
df_o <- psmelt(AvgRA99_O)
#group and calculate mean, sd and se for different taxonomic levels
grouped_o <- group_by(df_o, group, Phylum, Class, Order)
avgs_o <- dplyr::summarise(grouped_o, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))

avgs_filt <- filter(avgs_o, mean > 0)
# I want to read out this table for sup data
write.csv(avgs_filt, "SM_treat-time_fam_avg_sd_se .csv")


avgs_filt$ST <- as.factor(avgs_filt$group)
avgs_filt_st <- melt(data.frame(avgs_filt), id.vars=c("Order", "group"),
                     measure.vars=c("mean"))

#start
chisq = NULL
pvals = NULL
listofcats_sig = NULL
DataSet = df_o
taxa = as.factor(unique(df_o$Order))


for (cat in taxa) {
  new_df <- subset(DataSet, DataSet$Order == cat)
  new_df$ST <- as.factor(new_df$group)
  kw = kruskal.test(Abundance ~ ST, data=new_df)
  chisq = c(kw$statistic, chisq)
  pvals = c(kw$p.value, pvals)
  if (kw$p.value <= 0.05) {
    listofcats_sig = c(cat, listofcats_sig)
  }
}

pvals.bonf = p.adjust(pvals, method="bonferroni")
df_taxa = data.frame(rev(taxa), chisq, pvals, pvals.bonf)
write.csv(df_taxa, 'SM_treat-time_VST_Mean_FAM_KW.csv')

pvals.dunn = NULL
pvals.dunn.bonf = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$Order == cat)
  new_df$ST <- as.factor(new_df$group)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bonferroni")
  for (i in 1:length(dT$res$Comparison)) {
    print(cat)
    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bonf = c(dT$res$P.adj[i], pvals.dunn.bonf)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(dT$res$Comparison[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)
  }
}

df_taxa = data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bonf)

# df_taxa_red <- filter(df_taxa, comparison.dunn != c( "C38 - SM51","C51 - SM38"))
# 
# df_taxa_red$pvals.dunn.bonf <-  p.adjust(df_taxa_red$pvals.dunn, method="bonferroni")

tax_oi <- filter(df_taxa, pvals.dunn.bonf < 0.05)

write.csv(df_taxa, 'SM_treat-time_VST_Mean_FAM_KW_Dunn_OI.csv')



#################################################################





# Taxonomic bar chart figure
#################################################################

#group and calculate mean, sd and se for different taxonomic levels
groups <- group_by(df_o, group, Phylum, Class, Order)
groups_avgs <- dplyr::summarise(groups, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))

avgs_0p <- filter(groups_avgs, mean > 0)


# new df, change families
avgs_f <- avgs_0p
avgs_f$Order <- factor(avgs_f$Order)
avgs_f <- na.omit(avgs_f)
# remove families with %<5
avgs_5 <- filter(avgs_f, mean > 5)
avgs_3 <- filter(avgs_f, mean > 3)
avgs_2 <- filter(avgs_f, mean > 2)
avgs_1.5 <- filter(avgs_f, mean > 1.5)
avgs_1 <- filter(avgs_f, mean > 1)



# New facet label names for treatment variable
ad_labeller <- as_labeller( c("control" = "Control",
                              "larvae" = "Larva",
                              "low-ph" = "Low-pH",
                              "molybdate" = "Sodium Molybdate") )

# 12 colors
col_12 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')


# 16 colors
col_16 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f',
            '#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928','#999999','#4d4d4d',
            '#c994c7','#df65b0')


alpha.breaks <- unique(avgs_f$group)

alpha.ids <- c("[C-1+C-2]",  "[C-3]",
               "[SM-1+SM-2]", "[SM-3]")

# All treatments

ggplot(avgs_1.5, aes(x=group, y=mean, fill=Order))+
  geom_bar(stat="identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)),
                width=.4, position=position_dodge(.9)) +
  theme(legend.title=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"),
        legend.position="bottom") +
  xlab("Treatment-Age") + ylab("Mean Percent Abundance") +
  scale_fill_manual(values = col_16) +
  coord_cartesian(ylim = c(1.5,32)) +
  scale_x_discrete(breaks=alpha.breaks,
                   labels=alpha.ids) 


#################################################################




# Sulfate-reducing bacteria 
#################################################################

ps.sulf <- subset_taxa(ps.final, Class==c("Deltaproteobacteria"))

AvgRA_sulf_f = tax_glom(ps.sulf, taxrank="Family", NArm = FALSE)

df_sulf <- psmelt(AvgRA_sulf_f)

#group and calculate mean, sd and se for different taxonomic levels
grouped_sulf <- group_by(df_sulf, group, Phylum, Class, Order, Family)
avgs_sulf <- dplyr::summarise(grouped_sulf, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))

avgs_filt_sulf <- filter(avgs_sulf, mean > 0)


ggplot(avgs_filt_sulf, aes(x=group, y=mean, fill=Family))+
  geom_bar(stat="identity", position = "dodge", color="black") +
  geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)),
                width=.4, position=position_dodge(.9)) +
  theme(legend.title=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="bold"),
        legend.position="bottom") +
  scale_fill_manual(values = col_16) +
  xlab("Treatment-Age") + ylab("Mean Percent Abundance") +
  coord_cartesian(ylim = c(0.07,1.6)) +
  scale_x_discrete(breaks=alpha.breaks,
                   labels=alpha.ids)


# some exploratory things

ex1 <- subset_taxa(ps.final, Order==c("Desulfarculales"))
ex2 <- subset_taxa(ps.final, Order==c("Desulfovibrionales"))
ex3 <- subset_taxa(ps.final, Order==c("Desulfobacterales"))
ex4 <- subset_taxa(ps.final, Order==c("Desulfobulbaceae"))
ex5 <- subset_taxa(ps.final, Order==c("Campylobacterales"))

ex6 <- subset_taxa(ps.final, Class==c("Deltaproteobacteria"))
plot_bar(ps.sulf, fill="Family")


ex4 <- merge_phyloseq(ex1,ex2,ex3,ex5)

quartz()
plot_bar(ex4, fill="Family")



#################################################################




# Shell analysis 
#################################################################

setwd("C:/Users/Xanne/Desktop/Oyster microbes/shell-data/")

shell.dat <-read.csv("oyster_seed_shell_data.csv")
names(shell.dat)
unique(shell.dat$treatment)

shell.dat.red <- filter(shell.dat, shell.dat$treatment != "low-ph")



# reshaping data and initial statistical tests
# using filter to take out all rows that have "whole" as the entry in the stage column
whole <- filter(shell.dat.red, stage=="whole")
larval <- filter(shell.dat.red, stage=="larval")

# renaming larval column names
names(larval) <- c("image.id", "treatment", "bucket.no", "age", "stage", "larval.area.um2", "larval.perimeter.um")

# renaming whole column names
names(whole) <- c("image.id", "treatment", "bucket.no", "age", "stage", "whole.area.um2", "whole.perimeter.um")

# taking columns 1-4, 6, and 7 from larval and 6 and 7 from whole
final <- as.data.frame(c(larval[,c(1:4,6:7)], whole[,6:7]))

names(final)
# note that we no longer need the "stage" key column, instead larval and whole measurements have their own columns

# creating a column for just the juvenile part of the shell by subtracting the larval shell meaurements from the whole
final$juvenile.area.um2 <- final$whole.area.um2 - final$larval.area.um2
# new column that are treatment-age groups
final$group <- paste(final$treatment,final$age)


# going to use analysis of variance (ANOVA) to evaluate whether or not shell areas are statistically different between treatment groups

# calls aov() from base R stats comparing larval.area.um2 between treatments
larv.group.aov <- aov(larval.area.um2 ~ group, data = final)
summary(larv.group.aov)
#            Df    Sum Sq   Mean Sq F value  Pr(>F)   
# group        3 2.141e+09 713787732   5.586 0.00201 **
# Residuals   56 7.155e+09 127772402  

# Now do to the pairwise comparisons, also known as a posthoc test
larv.thsd <- TukeyHSD(larv.group.aov)
larv.thsd
#                             diff         lwr       upr     p adj
# control 52-control 38     11102.6716    920.3346 21285.009 0.0274539
# molybdate 38-control 38    9056.5176    742.6743 17370.361 0.0276535
# molybdate 52-control 38   10314.7421    132.4051 20497.079 0.0459900
# molybdate 38-control 52   -2046.1540 -12228.4910  8136.183 0.9508939
# molybdate 52-control 52    -787.9295 -12545.4795 10969.621 0.9979915
# molybdate 52-molybdate 38  1258.2245  -8924.1125 11440.561 0.9877760



# Now we repeat for the juvenile shell measurements
juvenile.group.aov <- aov(juvenile.area.um2 ~ group, data = final)
summary(juvenile.group.aov)
#               Df    Sum Sq   Mean Sq F value   Pr(>F)    
# group        3 3.359e+12 1.120e+12   43.78 1.06e-14 ***
# Residuals   56 1.432e+12 2.558e+10    

juvenile.thsd <- TukeyHSD(juvenile.group.aov)
juvenile.thsd
#                            diff        lwr        upr     p adj
# control 52-control 38      461362.1  297345.04 625379.244 0.0000000
# molybdate 38-control 38    289637.7  155718.27 423557.080 0.0000025
# molybdate 52-control 38    662279.4  498262.25 826296.463 0.0000000
# molybdate 38-control 52   -171724.5 -335741.57  -7707.361 0.0368055
# molybdate 52-control 52    200917.2   11526.58 390307.858 0.0335486
# molybdate 52-molybdate 38  372641.7  208624.58 536658.789 0.0000008



# # Okay, now we know that shell area is significantly different between the contorl and low-ph treatments, 
# # and between the molybdate and control groups
# whole.group.aov <- aov(whole.area.um2 ~ group, data = final)
# summary(whole.group.aov)
# #           Df    Sum Sq   Mean Sq F value   Pr(>F)    
# # group        3 3.480e+12 1.160e+12   44.51 7.64e-15 ***
# # Residuals   56 1.459e+12 2.606e+10         
# 
# whole.thsd <- TukeyHSD(whole.group.aov)
# whole.thsd
# #                            diff         lwr        upr     p adj
# # control 52-control 38      472464.8  306906.374 638023.249 0.0000000
# # molybdate 38-control 38    298694.2  163516.294 433872.091 0.0000016
# # molybdate 52-control 38    672594.1  507035.664 838152.539 0.0000000
# # molybdate 38-control 52   -173770.6 -339329.057  -8212.182 0.0361717
# # molybdate 52-control 52    200129.3    8958.873 391299.707 0.0368368
# # molybdate 52-molybdate 38  373899.9  208341.471 539458.347 0.0000010




# Shell Figure

shell.dat.fig <-read.csv("shell_dat_fig.csv")
shell.dat.fig <- filter(shell.dat.fig, stage != "whole")

shell.dat.fig$group <- paste(shell.dat.fig$treatment,shell.dat.fig$age)
shell.dat.fig$area.mm2 <- shell.dat.fig$area.um2/1000

# summary function below called here
area.mm.summary <- summarySE(shell.dat.fig, measurevar="area.mm2", groupvars=c("treatment","stage","age"))
area.mm.summary


# # New facet label names for treatment variable
# ad_labeller_shell <- as_labeller( c("larval" = "Larval",
#                               "juv" = "Juvenile",
#                               "whole" = "Whole",
#                               "38" = "38",
#                               "52" = "52") )

ggplot(shell.dat.fig, aes(x=treatment, y=area.mm2))+
  geom_boxplot() + geom_jitter() +
  theme(legend.title=element_blank()) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  #coord_cartesian(ylim = c(30,650)) +
  xlab("Treatment") + ylab(bquote('Area (mm' ^2*')')) +
  facet_grid(stage~age, scales="free_y") +
  stat_summary(fun=mean, geom="point", shape=20, size=7, color="red", fill="red")+
  scale_x_discrete(breaks=c("control","molybdate"),
                   labels=c("Control","Molybdate"))


# getting weird results from this so going to try something else


shell.dat.larv <- filter(shell.dat.fig, stage == "larval")
shell.dat.juv <- filter(shell.dat.fig, stage == "juvenile")

ggplot(shell.dat.larv, aes(x=treatment, y=area.mm2)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter()





  
################################################################












