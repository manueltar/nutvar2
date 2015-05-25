
# Set the working directory

setwd("~/Dropbox/Proyecto_NutVar2/nutvar2/data/training")

#install.packages("ggplot2")
library(ggplot2)

# Open the initial table of data

Matrix_of_data<-read.table(file="Matrix_snpeff_CCDS_plus_tags_added_gene_based_scores_NEW_RESULTS_DECEMBER_2014.txt" ,header=TRUE,sep ="\t",na.strings="NA")

str(Matrix_of_data)

# Select the initial subsets of common (NonPathogenic) and Pathogenic truncations. The Matrix_of_data$Location_Tag field indicates
# the origin of the truncation; when it is 1 it comes from the pool of variants from Exome projects, when it is 2 it comes
# from ClinVar and when it is 3 it is present in both groups. 

# We accept as Non Pathogenic variants the ones with a Matrix_of_data$Location_Tag equal to 1 
# and a Matrix_of_data$Credible_Tag..1.rare..0.not_credible.1.common. equal to 1

# We accept as Pathogenic variants the ones with a Matrix_of_data$Pathogenic_Tag >=1 
# and a Matrix_of_data$Credible_Tag..1.rare..0.not_credible.1.common. different from 1 (NaN, -1 or 0)

NonPathogenic<-subset(Matrix_of_data,Matrix_of_data$Location_Tag==1 & Matrix_of_data$Credible_Tag..1.rare..0.not_credible.1.common.==1)

Pathogenic1<-subset(Matrix_of_data,Matrix_of_data$Pathogenicity_Tag>=1 & is.na(Matrix_of_data$Credible_Tag..1.rare..0.not_credible.1.common.))
Pathogenic2<-subset(Matrix_of_data,Matrix_of_data$Pathogenicity_Tag>=1 & Matrix_of_data$Credible_Tag..1.rare..0.not_credible.1.common. !=1)

Pathogenic<-rbind(Pathogenic1,Pathogenic2)

check<-Pathogenic$Credible_Tag..1.rare..0.not_credible.1.common.
check2<-NonPathogenic$Credible_Tag..1.rare..0.not_credible.1.common.
check3<-NonPathogenic$Pathogenicity_Tag

# Now we subdivide the Pathogenic and NonPathogenic subsets in stop_gains (ratio_stop>0), frameshifts (ratio_frameshift>0) and  splice (ratio_splice>0) subsets

NonPathogenic.stop<-subset(NonPathogenic, NonPathogenic$ratioAffectedIsoforms_stop.gained >0)
NonPathogenic.frameshift<-subset(NonPathogenic, NonPathogenic$ratioAffectedIsoforms_frameshift >0)
NonPathogenic.splice<-subset(NonPathogenic, NonPathogenic$ratioAffectedIsoforms_splice >0)

Pathogenic.stop<-subset(Pathogenic, Pathogenic$ratioAffectedIsoforms_stop.gained >0)
Pathogenic.frameshift<-subset(Pathogenic, Pathogenic$ratioAffectedIsoforms_frameshift >0)
Pathogenic.splice<-subset(Pathogenic, Pathogenic$ratioAffectedIsoforms_splice >0)

# Stop, frameshift and splice subsets 

check.stop<-rbind(NonPathogenic.stop,Pathogenic.stop)
check.frameshift<-rbind(NonPathogenic.frameshift,Pathogenic.frameshift)
check.splice<-rbind(NonPathogenic.splice,Pathogenic.splice)

# Create the label factor with NonPathogenic/Pathogenic on a per type of truncation basis

  # First we create a factor grp with the labels Non-Pathogenic and Pathogenic 
  # encompasing all types of truncations

  # STOPs
  
  group_NP_stop<-NULL
  group_P_stop<-NULL
  
  for(i in 1:length(NonPathogenic.stop$DomainINFOAvailable))
  {
    group_NP_stop[i]<-"Non_Pathogenic_stop" 
  }
  for(i in 1:length(Pathogenic.stop$DomainINFOAvailable))
  {
    group_P_stop[i]<-"Pathogenic_stop" 
  }
  
  # frameshifts

  group_NP_frameshift<-NULL
  group_P_frameshift<-NULL
  
  for(i in 1:length(NonPathogenic.frameshift$DomainINFOAvailable))
  {
    group_NP_frameshift[i]<-"Non_Pathogenic_frameshift" 
  }
  for(i in 1:length(Pathogenic.frameshift$DomainINFOAvailable))
  {
    group_P_frameshift[i]<-"Pathogenic_frameshift" 
  }
  
  # splice

  group_NP_splice<-NULL
  group_P_splice<-NULL
  
  for(i in 1:length(NonPathogenic.splice$DomainINFOAvailable))
  {
    group_NP_splice[i]<-"Non_Pathogenic_splice" 
  }
  for(i in 1:length(Pathogenic.splice$DomainINFOAvailable))
  {
    group_P_splice[i]<-"Pathogenic_splice" 
  }
  
  # Now we combine all the labels

  group<-c(group_NP_stop,group_P_stop,group_NP_frameshift,group_P_frameshift,group_NP_splice,group_P_splice)
  grp<-factor(group,levels=c("Non_Pathogenic_stop", "Pathogenic_stop", 
                              "Non_Pathogenic_frameshift", "Pathogenic_frameshift", 
                              "Non_Pathogenic_splice", "Pathogenic_splice"))
summary(grp)  
levels(grp)
                 
# BOXPLOT


# NMD

check_NMD_stop<-check.stop$ratioAffectedIsoformsTargetedbyNMD
check_NMD_frameshift<-check.frameshift$ratioAffectedIsoformsTargetedbyNMD
check_NMD_splice<-check.splice$ratioAffectedIsoformsTargetedbyNMD

# For frameshifts NMD is imputed to 0

check_NMD_frameshift_imputed<-0*check_NMD_frameshift


NMD<-c(check_NMD_stop,check_NMD_frameshift_imputed,check_NMD_splice)

  # Here we calculate the mean for each of the groups

  means.NMD <- aggregate(NMD~grp, data=NULL, mean)
  means.NMD.rounded <-round(means.NMD$NMD,2)
  
  colour.factor<-factor(c("NP","P","NP"))
  
# Here we boxplot the data
  
  pdf(file="boxplot_NMD.pdf")
  ggplot(data=NULL, aes(x=grp, y=NMD, fontsize=6)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 1.5, label = means.NMD.rounded, size=5, colour="darkred")
  dev.off()


# NMD_derived

check_NMD_derived_stop<-check.stop$ratioAffectedIsoformsTargetedby_derived_NMD
check_NMD_derived_frameshift<-check.frameshift$ratioAffectedIsoformsTargetedby_derived_NMD
check_NMD_derived_splice<-check.splice$ratioAffectedIsoformsTargetedby_derived_NMD

# For stop and splice NMD_derived is imputed to 0

check_NMD_derived_stop_imputed<-0*check_NMD_stop
check_NMD_derived_splice_imputed<-0*check_NMD_splice


NMD_derived<-c(check_NMD_derived_stop_imputed,check_NMD_derived_frameshift,check_NMD_derived_splice_imputed)

  # Here we calculate the mean for each of the groups
  
  means.NMD_derived <- aggregate(NMD_derived~grp, data=NULL, mean)

  means.NMD_derived.rounded <-round(means.NMD_derived$NMD_derived,2)
  
  # Here we boxplot the data
  
  pdf(file="boxplot_NMD_derived.pdf")
  ggplot(data=NULL, aes(x=grp, y=NMD_derived)) + geom_boxplot(fill="steelblue") +
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=4,show_guide = FALSE) + 
    annotate("text", x = 1:6, y = 2.5, label = means.NMD_derived.rounded, size=5, colour="darkred")
  dev.off()
  

# Longest

check_Longest_stop<-check.stop$IsWithinLongestCCDS
check_Longest_frameshift<-check.frameshift$IsWithinLongestCCDS
check_Longest_splice<-check.splice$IsWithinLongestCCDS


Longest<-c(check_Longest_stop,check_Longest_frameshift,check_Longest_splice)

  # Here we calculate the mean for each of the groups
  
  means.Longest <- aggregate(Longest~grp, data=NULL, mean)
  means.Longest.rounded <-round(means.Longest$Longest,2)
  
  # Here we boxplot the data
  
  pdf(file="boxplot_Longest.pdf")
  ggplot(data=NULL, aes(x=grp, y=Longest)) + geom_boxplot( fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 2.5, label = means.Longest.rounded, size=5, colour="darkred")
  dev.off()


# APPRIS

check_APPRIS_stop<-check.stop$IsPrincipalIsoformAffected
check_APPRIS_frameshift<-check.frameshift$IsPrincipalIsoformAffected
check_APPRIS_splice<-check.splice$IsPrincipalIsoformAffected

APPRIS<-c(check_APPRIS_stop,check_APPRIS_frameshift,check_APPRIS_splice)

# Impute APPRIS to Longest isoform when NaN

APPRIS_imputed<-APPRIS

for(i in 1:length(APPRIS_imputed))
{
  if(is.na(APPRIS_imputed[i]))
  {
    APPRIS_imputed[i]<-Longest[i]
  }  
}

  # Here we calculate the mean for each of the groups
  
  means.APPRIS <- aggregate(APPRIS~grp, data=NULL, mean)
  means.APPRIS.rounded <-round(means.APPRIS$APPRIS,2)
  
  # Here we boxplot the data
  
  pdf(file="boxplot_APPRIS.pdf")
  ggplot(data=NULL, aes(x=grp, y=APPRIS)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 2.5, label = means.APPRIS.rounded, size=5, colour="darkred")
  dev.off()
  
  # Here we calculate the mean for each of the groups
  
  means.APPRIS_imputed <- aggregate(APPRIS_imputed~grp, data=NULL, mean)
  means.APPRIS_imputed.rounded <-round(means.APPRIS_imputed$APPRIS_imputed,2)
  
  # Here we boxplot the data
  
  pdf(file="boxplot_APPRIS_imputed.pdf")
  ggplot(data=NULL, aes(x=grp, y=APPRIS_imputed)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 2.5, label = means.APPRIS_imputed.rounded, size=5, colour="darkred")
  dev.off()


# Pervasive

check_Pervasive_stop<-check.stop$IsWithinPervasiveIsoform
check_Pervasive_frameshift<-check.frameshift$IsWithinPervasiveIsoform
check_Pervasive_splice<-check.splice$IsWithinPervasiveIsoform

Pervasive<-c(check_Pervasive_stop,check_Pervasive_frameshift,check_Pervasive_splice)

# Impute Pervasive to Longest isoform when NaN

Pervasive_imputed<-Pervasive

for(i in 1:length(Pervasive_imputed))
{
  if(is.na(Pervasive_imputed[i]))
  {
    Pervasive_imputed[i]<-Longest[i]
  }  
}

  # Here we calculate the mean for each of the groups
  
  means.Pervasive <- aggregate(Pervasive~grp, data=NULL, mean)
  means.Pervasive.rounded <-round(means.Pervasive$Pervasive,2)
  
  # Here we boxplot the data
  
  pdf(file="boxplot_Pervasive.pdf")
  ggplot(data=NULL, aes(x=grp, y=Pervasive)) + geom_boxplot( fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 2.5, label = means.Pervasive.rounded, size=5, colour="darkred")
  dev.off()
  
  # Here we calculate the mean for each of the groups
  
  means.Pervasive_imputed <- aggregate(Pervasive_imputed~grp, data=NULL, mean)
  means.Pervasive_imputed.rounded <-round(means.Pervasive_imputed$Pervasive_imputed,2)
  
  # Here we boxplot the data
  
  pdf(file="boxplot_Pervasive_imputed.pdf")
  ggplot(data=NULL, aes(x=grp, y=Pervasive_imputed)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 2.5, label = means.Pervasive_imputed.rounded, size=5, colour="darkred")
  dev.off()

# Longest_length

check_Longest_length_stop<-check.stop$LongestCCDSLength
check_Longest_length_frameshift<-check.frameshift$LongestCCDSLength
check_Longest_length_splice<-check.splice$LongestCCDSLength


Longest_length<-c(check_Longest_length_stop,check_Longest_length_frameshift,check_Longest_length_splice)

  # Here we calculate the mean for each of the groups
  
  means.Longest_length <- aggregate(Longest_length~grp, data=NULL, mean)
  means.Longest_length.rounded <-round(means.Longest_length$Longest_length,2)
  
  # Here we boxplot the data
  
    pdf(file="boxplot_Longest_length.pdf")
    ggplot(data=NULL, aes(x=grp, y=Longest_length)) + geom_boxplot(fill="steelblue") +
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=4,show_guide = FALSE) + 
    annotate("text", x = 1:6, y = 35000, label = means.Longest_length.rounded, size=5, colour="darkred")
    dev.off()

# Percentage_seq

check_Percentage_seq_stop<-check.stop$PercentagePrincipalOrLongestCCDSAffected
check_Percentage_seq_frameshift<-check.frameshift$PercentagePrincipalOrLongestCCDSAffected
check_Percentage_seq_splice<-check.splice$PercentagePrincipalOrLongestCCDSAffected

# Issue Splice variants for which percentage is NaN investigate why this is so, switching off this imputation

Percentage_seq<-c(check_Percentage_seq_stop,check_Percentage_seq_frameshift,check_Percentage_seq_splice)
Percentage_seq_imputed<-Percentage_seq

Lebron<-is.na(Percentage_seq_imputed)
summary(Lebron)
for(i in 1:length(Percentage_seq_imputed))
{
  if(is.na( Percentage_seq_imputed[i]))
  {
    Percentage_seq_imputed[i]<-0
  }  
}
Lebron<-is.na(Percentage_seq_imputed)
summary(Lebron)

  # Here we calculate the mean for each of the groups
  
  means.Percentage_seq <- aggregate(Percentage_seq~grp, data=NULL, mean)
  means.Percentage_seq.rounded <-round(means.Percentage_seq$Percentage_seq,2)
  
  # Here we boxplot the data
  
  pdf(file="boxplot_Percentage_seq.pdf")
  ggplot(data=NULL, aes(x=grp, y=Percentage_seq)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 110, label = means.Percentage_seq.rounded, size=5, colour="darkred")
  dev.off()

  # Here we calculate the mean for each of the groups
  

  # Here we calculate the mean for each of the groups
  
  means.Percentage_seq_imputed <- aggregate(Percentage_seq_imputed~grp, data=NULL, mean)
  means.Percentage_seq_imputed.rounded <-round(means.Percentage_seq_imputed$Percentage_seq_imputed,2)
  
  # Here we boxplot the data
  
  pdf(file="boxplot_Percentage_seq_imputed.pdf")
  ggplot(data=NULL, aes(x=grp, y=Percentage_seq_imputed)) + geom_boxplot(fill="steelblue") +
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=4,show_guide = FALSE) + 
    annotate("text", x = 1:6, y = 110, label = means.Percentage_seq_imputed.rounded, size=5, colour="darkred")  
  dev.off()


# Domain_INFO

check_Domain_INFO_stop<-check.stop$DomainINFOAvailable
check_Domain_INFO_frameshift<-check.frameshift$DomainINFOAvailable
check_Domain_INFO_splice<-check.splice$DomainINFOAvailable


Domain_INFO<-c(check_Domain_INFO_stop,check_Domain_INFO_frameshift,check_Domain_INFO_splice)

  # Here we calculate the mean for each of the groups
  
  means.Domain_INFO <- aggregate(Domain_INFO~grp, data=NULL, mean)
  means.Domain_INFO.rounded <-round(means.Domain_INFO$Domain_INFO,2)
  
  # Here we boxplot the data
  
  pdf(file="boxplot_Domain_INFO.pdf")
  ggplot(data=NULL, aes(x=grp, y=Domain_INFO)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 1.5, label = means.Domain_INFO.rounded, size=5, colour="darkred")
  dev.off()

# Domain_percentage

check_Domain_percentage_stop<-check.stop$PercentageOfDomainPositionsAffected
check_Domain_percentage_frameshift<-check.frameshift$PercentageOfDomainPositionsAffected
check_Domain_percentage_splice<-check.splice$PercentageOfDomainPositionsAffected

Domain_percentage<-c(check_Domain_percentage_stop,check_Domain_percentage_frameshift,check_Domain_percentage_splice)

# Impute Domain_percentage to seq Percetange when missing

Domain_percentage_imputed<-Domain_percentage
Harden<-is.na(Domain_percentage_imputed)
summary(Harden)
for(i in 1:length(Domain_percentage_imputed))
{
  if(is.na(Domain_percentage_imputed[i]))
  {
    Domain_percentage_imputed[i]<-Percentage_seq_imputed[i]
  }  
}
Harden<-is.na(Domain_percentage_imputed)
summary(Harden)

  # Here we calculate the mean for each of the groups
  
  means.Domain_percentage <- aggregate(Domain_percentage~grp, data=NULL, mean)
  means.Domain_percentage.rounded <-round(means.Domain_percentage$Domain_percentage,2)
  
  # Here we boxplot the data
  
  pdf(file="boxplot_Domain_percentage.pdf")
  ggplot(data=NULL, aes(x=grp, y=Domain_percentage)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 110, label = means.Domain_percentage.rounded, size=5, colour="darkred")
  dev.off()
  
  # Here we calculate the mean for each of the groups
  
  means.Domain_percentage_imputed <- aggregate(Domain_percentage_imputed~grp, data=NULL, mean)
  means.Domain_percentage_imputed.rounded <-round(means.Domain_percentage_imputed$Domain_percentage_imputed,2)

  # Here we boxplot the data
  
  pdf(file="boxplot_Domain_percentage_imputed.pdf")
  ggplot(data=NULL, aes(x=grp, y=Domain_percentage_imputed)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 110, label = means.Domain_percentage_imputed.rounded, size=5, colour="darkred")  
  dev.off()

# MaxPerDomain

check_MaxPerDomain_stop<-check.stop$maxPercDomainAffected
check_MaxPerDomain_frameshift<-check.frameshift$maxPercDomainAffected
check_MaxPerDomain_splice<-check.splice$maxPercDomainAffected

MaxPerDomain<-c(check_MaxPerDomain_stop,check_MaxPerDomain_frameshift,check_MaxPerDomain_splice)

# Impute MaxPerDomain to 100% when NaN

MaxPerDomain_imputed<-MaxPerDomain

for(i in 1:length(MaxPerDomain_imputed))
{
  if(is.na(MaxPerDomain_imputed[i]))
  {
    MaxPerDomain_imputed[i]<-100
  }  
}

  # Here we calculate the mean for each of the groups
  
  means.MaxPerDomain <- aggregate(MaxPerDomain~grp, data=NULL, mean)
  means.MaxPerDomain.rounded <-round(means.MaxPerDomain$MaxPerDomain,2)

  
  # Here we boxplot the data
  
  pdf(file="boxplot_MaxPerDomain.pdf")
  ggplot(data=NULL, aes(x=grp, y=MaxPerDomain)) + geom_boxplot(fill="steelblue") +
    stat_summary(fun.y=mean, colour="darkred", geom="point", 
                 shape=18, size=4,show_guide = FALSE) + 
    annotate("text", x = 1:6, y = 110, label = means.MaxPerDomain.rounded, size=5, colour="darkred")  
  dev.off()
  
  # Here we calculate the mean for each of the groups
  
  means.MaxPerDomain_imputed <- aggregate(MaxPerDomain_imputed~grp, data=NULL, mean)
  means.MaxPerDomain_imputed.rounded <-round(means.MaxPerDomain_imputed$MaxPerDomain_imputed,2)

  
  # Here we boxplot the data
  
  pdf(file="boxplot_MaxPerDomain_imputed.pdf")
  ggplot(data=NULL, aes(x=grp, y=MaxPerDomain_imputed)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 110, label = means.MaxPerDomain_imputed.rounded, size=5, colour="darkred") 
  dev.off()


# number_of_100_damaged

check_number_of_100_damaged_stop<-check.stop$NumberOfDomains100Damage
check_number_of_100_damaged_frameshift<-check.frameshift$NumberOfDomains100Damage
check_number_of_100_damaged_splice<-check.splice$NumberOfDomains100Damage

number_of_100_damaged<-c(check_number_of_100_damaged_stop,check_number_of_100_damaged_frameshift,check_number_of_100_damaged_splice)

# Impute number_of_100_damaged to 1 when NaN

number_of_100_damaged_imputed<-number_of_100_damaged

for(i in 1:length(number_of_100_damaged_imputed))
{
  if(is.na(number_of_100_damaged_imputed[i]))
  {
    number_of_100_damaged_imputed[i]<-1
  }  
}

  # Here we calculate the mean for each of the groups
  
  means.number_of_100_damaged <- aggregate(number_of_100_damaged~grp, data=NULL, mean)
  means.number_of_100_damaged.rounded <-round(means.number_of_100_damaged$number_of_100_damaged,2)

  # Here we boxplot the data
  
  pdf(file="boxplot_number_of_100_damaged.pdf")
  ggplot(data=NULL, aes(x=grp, y=number_of_100_damaged)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 400, label = means.number_of_100_damaged.rounded, size=5, colour="darkred") 
  dev.off()
  
  # Here we calculate the mean for each of the groups
  
  means.number_of_100_damaged_imputed <- aggregate(number_of_100_damaged_imputed~grp, data=NULL, mean)
  means.number_of_100_damaged_imputed.rounded <-round(means.number_of_100_damaged_imputed$number_of_100_damaged_imputed,2)

  # Here we boxplot the data
  
  pdf(file="boxplot_number_of_100_damaged_imputed.pdf")
  ggplot(data=NULL, aes(x=grp, y=number_of_100_damaged_imputed)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 400, label = means.number_of_100_damaged_imputed.rounded, size=5, colour="darkred") 
  dev.off()

# Domain_matched

check_Domain_matched_stop<-check.stop$DomainMatched
check_Domain_matched_frameshift<-check.frameshift$DomainMatched
check_Domain_matched_splice<-check.splice$DomainMatched

Domain_matched<-c(check_Domain_matched_stop,check_Domain_matched_frameshift,check_Domain_matched_splice)

# Impute Domain_matched to 1 when NaN

Domain_matched_imputed<-Domain_matched

for(i in 1:length(Domain_matched_imputed))
{
  if(is.na(Domain_matched_imputed[i]))
  {
    Domain_matched_imputed[i]<-1
  }  
}

  # Here we calculate the mean for each of the groups
  
  means.Domain_matched <- aggregate(Domain_matched~grp, data=NULL, mean)
  means.Domain_matched.rounded <-round(means.Domain_matched$Domain_matched,2)

  
  # Here we boxplot the data
  
  pdf(file="boxplot_Domain_matched.pdf")
  ggplot(data=NULL, aes(x=grp, y=Domain_matched)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 1.5, label = means.Domain_matched.rounded, size=5, colour="darkred") 

  dev.off()
  
  # Here we calculate the mean for each of the groups
  
  means.Domain_matched_imputed <- aggregate(Domain_matched_imputed~grp, data=NULL, mean)
  means.Domain_matched_imputed.rounded <-round(means.Domain_matched_imputed$Domain_matched_imputed,2)

  
  # Here we boxplot the data
  
  pdf(file="boxplot_Domain_matched_imputed.pdf")
  ggplot(data=NULL, aes(x=grp, y=Domain_matched_imputed)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 1.5, label = means.Domain_matched_imputed.rounded, size=5, colour="darkred") 
  dev.off()

########################################### SITE  ####################################################

# Site_INFO

check_Site_INFO_stop<-check.stop$SiteINFOAvailable
check_Site_INFO_frameshift<-check.frameshift$SiteINFOAvailable
check_Site_INFO_splice<-check.splice$SiteINFOAvailable


Site_INFO<-c(check_Site_INFO_stop,check_Site_INFO_frameshift,check_Site_INFO_splice)

# Here we calculate the mean for each of the groups

means.Site_INFO <- aggregate(Site_INFO~grp, data=NULL, mean)
means.Site_INFO.rounded <-round(means.Site_INFO$Site_INFO,2)

# Here we boxplot the data

pdf(file="boxplot_Site_INFO.pdf")
ggplot(data=NULL, aes(x=grp, y=Site_INFO)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 1.5, label = means.Site_INFO.rounded, size=5, colour="darkred")
dev.off()

# Site_percentage

check_Site_percentage_stop<-check.stop$PercentageOfSitePositionsAffected
check_Site_percentage_frameshift<-check.frameshift$PercentageOfSitePositionsAffected
check_Site_percentage_splice<-check.splice$PercentageOfSitePositionsAffected

Site_percentage<-c(check_Site_percentage_stop,check_Site_percentage_frameshift,check_Site_percentage_splice)

# Impute Site_percentage to seq Percetange when missing

Site_percentage_imputed<-Site_percentage

for(i in 1:length(Site_percentage_imputed))
{
  if(is.na(Site_percentage_imputed[i]))
  {
    Site_percentage_imputed[i]<-Percentage_seq_imputed[i]
  }  
}

# Here we calculate the mean for each of the groups

means.Site_percentage <- aggregate(Site_percentage~grp, data=NULL, mean)
means.Site_percentage.rounded <-round(means.Site_percentage$Site_percentage,2)

# Here we boxplot the data

pdf(file="boxplot_Site_percentage.pdf")
ggplot(data=NULL, aes(x=grp, y=Site_percentage)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 110, label = means.Site_percentage.rounded, size=5, colour="darkred")
dev.off()

# Here we calculate the mean for each of the groups

means.Site_percentage_imputed <- aggregate(Site_percentage_imputed~grp, data=NULL, mean)
means.Site_percentage_imputed.rounded <-round(means.Site_percentage_imputed$Site_percentage_imputed,2)

# Here we boxplot the data

pdf(file="boxplot_Site_percentage_imputed.pdf")
ggplot(data=NULL, aes(x=grp, y=Site_percentage_imputed)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 110, label = means.Site_percentage_imputed.rounded, size=5, colour="darkred")  
dev.off()

# MaxPerSite

check_MaxPerSite_stop<-check.stop$maxPercSiteAffected
check_MaxPerSite_frameshift<-check.frameshift$maxPercSiteAffected
check_MaxPerSite_splice<-check.splice$maxPercSiteAffected

MaxPerSite<-c(check_MaxPerSite_stop,check_MaxPerSite_frameshift,check_MaxPerSite_splice)

# Impute MaxPerSite to 100% when NaN

MaxPerSite_imputed<-MaxPerSite

for(i in 1:length(MaxPerSite_imputed))
{
  if(is.na(MaxPerSite_imputed[i]))
  {
    MaxPerSite_imputed[i]<-100
  }  
}

# Here we calculate the mean for each of the groups

means.MaxPerSite <- aggregate(MaxPerSite~grp, data=NULL, mean)
means.MaxPerSite.rounded <-round(means.MaxPerSite$MaxPerSite,2)


# Here we boxplot the data

pdf(file="boxplot_MaxPerSite.pdf")
ggplot(data=NULL, aes(x=grp, y=MaxPerSite)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 110, label = means.MaxPerSite.rounded, size=5, colour="darkred")  
dev.off()

# Here we calculate the mean for each of the groups

means.MaxPerSite_imputed <- aggregate(MaxPerSite_imputed~grp, data=NULL, mean)
means.MaxPerSite_imputed.rounded <-round(means.MaxPerSite_imputed$MaxPerSite_imputed,2)


# Here we boxplot the data

pdf(file="boxplot_MaxPerSite_imputed.pdf")
ggplot(data=NULL, aes(x=grp, y=MaxPerSite_imputed)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 110, label = means.MaxPerSite_imputed.rounded, size=5, colour="darkred") 
dev.off()


# number_of_100_damaged_Site

check_number_of_100_damaged_Site_stop<-check.stop$NumberOfSites100Damage
check_number_of_100_damaged_Site_frameshift<-check.frameshift$NumberOfSites100Damage
check_number_of_100_damaged_Site_splice<-check.splice$NumberOfSites100Damage

number_of_100_damaged_Site<-c(check_number_of_100_damaged_Site_stop,check_number_of_100_damaged_Site_frameshift,check_number_of_100_damaged_Site_splice)

# Impute number_of_100_damaged_Site to 1 when NaN

number_of_100_damaged_Site_imputed<-number_of_100_damaged_Site

for(i in 1:length(number_of_100_damaged_Site_imputed))
{
  if(is.na(number_of_100_damaged_Site_imputed[i]))
  {
    number_of_100_damaged_Site_imputed[i]<-1
  }  
}

# Here we calculate the mean for each of the groups

means.number_of_100_damaged_Site <- aggregate(number_of_100_damaged_Site~grp, data=NULL, mean)
means.number_of_100_damaged_Site.rounded <-round(means.number_of_100_damaged_Site$number_of_100_damaged_Site,2)

# Here we boxplot the data

pdf(file="boxplot_number_of_100_damaged_Site.pdf")
ggplot(data=NULL, aes(x=grp, y=number_of_100_damaged_Site)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 75, label = means.number_of_100_damaged_Site.rounded, size=5, colour="darkred") 
dev.off()

# Here we calculate the mean for each of the groups

means.number_of_100_damaged_Site_imputed <- aggregate(number_of_100_damaged_Site_imputed~grp, data=NULL, mean)
means.number_of_100_damaged_Site_imputed.rounded <-round(means.number_of_100_damaged_Site_imputed$number_of_100_damaged_Site_imputed,2)

# Here we boxplot the data

pdf(file="boxplot_number_of_100_damaged_Site_imputed.pdf")
ggplot(data=NULL, aes(x=grp, y=number_of_100_damaged_Site_imputed)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 75, label = means.number_of_100_damaged_Site_imputed.rounded, size=5, colour="darkred") 
dev.off()


# Site_matched

check_Site_matched_stop<-check.stop$SiteMatched
check_Site_matched_frameshift<-check.frameshift$SiteMatched
check_Site_matched_splice<-check.splice$SiteMatched

Site_matched<-c(check_Site_matched_stop,check_Site_matched_frameshift,check_Site_matched_splice)

# Impute Site_matched to 0 when NaN

Site_matched_imputed<-Site_matched

for(i in 1:length(Site_matched_imputed))
{
  if(is.na(Site_matched_imputed[i]))
  {
    Site_matched_imputed[i]<-0
  }  
}

# Here we calculate the mean for each of the groups

means.Site_matched <- aggregate(Site_matched~grp, data=NULL, mean)
means.Site_matched.rounded <-round(means.Site_matched$Site_matched,2)


# Here we boxplot the data

pdf(file="boxplot_Site_matched.pdf")
ggplot(data=NULL, aes(x=grp, y=Site_matched)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 1.5, label = means.Site_matched.rounded, size=5, colour="darkred") 

dev.off()

# Here we calculate the mean for each of the groups

means.Site_matched_imputed <- aggregate(Site_matched_imputed~grp, data=NULL, mean)
means.Site_matched_imputed.rounded <-round(means.Site_matched_imputed$Site_matched_imputed,2)


# Here we boxplot the data

pdf(file="boxplot_Site_matched_imputed.pdf")
ggplot(data=NULL, aes(x=grp, y=Site_matched_imputed)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 1.5, label = means.Site_matched_imputed.rounded, size=5, colour="darkred") 
dev.off()

# IImunity

check_IImunity_stop<-check.stop$IsInnateImmunity
check_IImunity_frameshift<-check.frameshift$IsInnateImmunity
check_IImunity_splice<-check.splice$IsInnateImmunity


IImunity<-c(check_IImunity_stop,check_IImunity_frameshift,check_IImunity_splice)

# Here we calculate the mean for each of the groups

means.IImunity <- aggregate(IImunity~grp, data=NULL, mean)
means.IImunity_matched.rounded <-round(means.IImunity$IImunity,2)

# Here we boxplot the data

pdf(file="boxplot_IImunity.pdf")
ggplot(data=NULL, aes(x=grp, y=IImunity)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 1.5, label =  means.IImunity_matched.rounded, size=5, colour="darkred") 

dev.off()

# Antiviral

check_Antiviral_stop<-check.stop$IsAntiviral
check_Antiviral_frameshift<-check.frameshift$IsAntiviral
check_Antiviral_splice<-check.splice$IsAntiviral


Antiviral<-c(check_Antiviral_stop,check_Antiviral_frameshift,check_Antiviral_splice)

# Here we calculate the mean for each of the groups

means.Antiviral <- aggregate(Antiviral~grp, data=NULL, mean)
means.Antiviral_matched.rounded <-round(means.Antiviral$Antiviral,2)

# Here we boxplot the data

pdf(file="boxplot_Antiviral.pdf")
ggplot(data=NULL, aes(x=grp, y=Antiviral)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 1.5, label =  means.Antiviral_matched.rounded, size=5, colour="darkred") 

dev.off()

# ISG

check_ISG_stop<-check.stop$IsISG
check_ISG_frameshift<-check.frameshift$IsISG
check_ISG_splice<-check.splice$IsISG


ISG<-c(check_ISG_stop,check_ISG_frameshift,check_ISG_splice)

# Here we calculate the mean for each of the groups

means.ISG <- aggregate(ISG~grp, data=NULL, mean)
means.ISG_matched.rounded <-round(means.ISG$ISG,2)

# Here we boxplot the data

pdf(file="boxplot_ISG.pdf")
ggplot(data=NULL, aes(x=grp, y=ISG)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 1.5, label =  means.ISG_matched.rounded, size=5, colour="darkred") 

dev.off()

# OMIMrecessive

check_OMIMrecessive_stop<-check.stop$IsOMIMrecessive
check_OMIMrecessive_frameshift<-check.frameshift$IsOMIMrecessive
check_OMIMrecessive_splice<-check.splice$IsOMIMrecessive


OMIMrecessive<-c(check_OMIMrecessive_stop,check_OMIMrecessive_frameshift,check_OMIMrecessive_splice)

# Here we calculate the mean for each of the groups

means.OMIMrecessive <- aggregate(OMIMrecessive~grp, data=NULL, mean)
means.OMIMrecessive_matched.rounded <-round(means.OMIMrecessive$OMIMrecessive,2)

# Here we boxplot the data

pdf(file="boxplot_OMIMrecessive.pdf")
ggplot(data=NULL, aes(x=grp, y=OMIMrecessive)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 1.5, label =  means.OMIMrecessive_matched.rounded, size=5, colour="darkred") 

dev.off()

# pRDG

check_pRDG_stop<-check.stop$pRDG_score
check_pRDG_frameshift<-check.frameshift$pRDG_score
check_pRDG_splice<-check.splice$pRDG_score


pRDG<-c(check_pRDG_stop,check_pRDG_frameshift,check_pRDG_splice)

# Here we calculate the mean for each of the groups

means.pRDG <- aggregate(pRDG~grp, data=NULL, mean)
means.pRDG_matched.rounded <-round(means.pRDG$pRDG,2)

# Here we boxplot the data

pdf(file="boxplot_pRDG.pdf")
ggplot(data=NULL, aes(x=grp, y=pRDG)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 1.5, label =  means.pRDG_matched.rounded, size=5, colour="darkred") 

dev.off()

# RVIS

check_RVIS_stop<-check.stop$RVIS_score
check_RVIS_frameshift<-check.frameshift$RVIS_score
check_RVIS_splice<-check.splice$RVIS_score


RVIS<-c(check_RVIS_stop,check_RVIS_frameshift,check_RVIS_splice)

# Here we calculate the mean for each of the groups

means.RVIS <- aggregate(RVIS~grp, data=NULL, mean)
means.RVIS_matched.rounded <-round(means.RVIS$RVIS,2)

# Here we boxplot the data

pdf(file="boxplot_RVIS.pdf")
ggplot(data=NULL, aes(x=grp, y=RVIS)) + geom_boxplot(fill="steelblue") +
  stat_summary(fun.y=mean, colour="darkred", geom="point", 
               shape=18, size=4,show_guide = FALSE) + 
  annotate("text", x = 1:6, y = 20, label =  means.RVIS_matched.rounded, size=5, colour="darkred") 

dev.off()

# Once imputed NaN, we export the table

result<-cbind(NMD,NMD_derived,Longest,APPRIS,Pervasive,Longest_length,Percentage_seq,Domain_INFO,Domain_percentage,MaxPerDomain,number_of_100_damaged,Domain_matched,Site_INFO,Site_percentage,MaxPerSite,number_of_100_damaged_Site,Site_matched,IImunity,Antiviral,ISG,OMIMrecessive,pRDG,RVIS,group)

result_imputed<-cbind(NMD,NMD_derived,Longest,APPRIS_imputed,Pervasive_imputed,Longest_length,Percentage_seq_imputed,Domain_INFO,Domain_percentage_imputed,MaxPerDomain_imputed,number_of_100_damaged_imputed,Domain_matched_imputed,Site_INFO,Site_percentage_imputed,MaxPerSite_imputed,number_of_100_damaged_Site_imputed,Site_matched_imputed,IImunity,Antiviral,ISG,OMIMrecessive,pRDG,RVIS,group)
Lebron<-is.na(Percentage_seq_imputed)
summary(Lebron)
write.table(result,"Table.txt")
write.table(result_imputed,"Table_imputed.txt")
