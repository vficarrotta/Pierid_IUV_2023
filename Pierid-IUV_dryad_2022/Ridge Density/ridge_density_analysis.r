# Ridge Density and associated analyses
# Random effects linear model
# anova
# plot of ridge Density

################################################
##### Direct measures of ridge Density dataset #
################################################

master <- read.csv('D:/Pierid_Phylogenetics/ridge_dens_master_v2.csv')

###################
#  Summary Stats  #
###################
# Each species cover and ground scales, separately
library('rapportools')

# wrangling
master$sp_sc <- paste0(master$species, '_', master$type)
redu <- data.frame(master$density, master$sp_sc)
datsum <- split(redu, f=redu$master.sp_sc)
sumsum <- lapply(datsum, FUN=summary)

############
# st dev   #
############
stdat <- unlist(lapply(datsum, function(i) sd(i[[1]])))
write.csv(stdat, 'D:/Pierid_Phylogenetics/ridge density stuff/stdevRD.csv')

############
#  Range   #
############
datrange <- unlist(lapply(datsum, function(i) range(i[[1]]))) # method of specifying column of dataframe when using lapply

datrange1 <- read.csv('D:/Pierid_Phylogenetics/ridge density stuff/rangesRD.csv')

mean(datrange1$Cover)/mean(datrange1$Ground)
  # 2.876077 times the range of cover scales density vs ground scales density

#############################
#  tree of SEM IUV species  #
#############################

library('phytools'); library('ggplot2'); library('lme4')

tree <- ladderize(read.tree('C:/Users/vfica/Documents/Pierid_Phylogenetics/Tree_5-2020/Pieridae_7114_final_names_trim_dates.tre'))
keepers <- unique(master$name)
keepers <- gsub('  ', ' ', keepers); keepers <- gsub(' ', '_', keepers)
tree$tip.label <- gsub('Eronia_leda', 'Afrodryas_leda', tree$tip.label) # swap Afrodryas for Eronia generic name
IUV.tree <- keep.tip(tree, keepers)

#########################
#       Ridge Plots     #
#########################
library('ggplot2'); library('dplyr')

# Plot: cover vs ground scale by Species
ggplot(master) +
  geom_violin(scale='width', trim=F, adjust=1.3, width=0.3, aes(x=as.factor(type), y=density)) +
  # geom_boxplot(width=0.1, aes(x=as.factor(name_st), y=density)) +
  geom_point(position=position_jitter(w=0.05), size=0.5, aes(x=as.factor(type), y=density, color=name)) +
  xlab("Scale Layer") +
  ylab('Ridge Density (ridges/micrometer)') +
  theme_classic()

# Plot: ratio of cover to ground scale Density
spec.order <- c('Eurema lisa', 'Dercas verhuelli', 'Gonepteryx mahaguru', 'Phoebis sennae', 'Zerene cesonia', 'Colias eurytheme', 'Afrodryas leda', 'Hebomoia glaucippe', 'Zegris  fausti', 'Ixias pyrene', 'Colotis danae')

ggplot(master, aes(x=name, y=dens_rat, color=dens_rat)) +
    geom_violin(scale='width', trim=F, adjust=1.3, width=0.4, fill='black') +
    geom_point(position=position_jitter(w=0.1), size=0.9, alpha=0.6) +
    scale_color_gradientn(colors=c('dark green', high='yellow'), limits=c(1, max(master$dens_rat))) +  # manually declare colors
    xlab('') +
    ylab('Density Ratio (Cover/Ground)') +
    theme_light() +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=-0.01)) + # modify x axis tick text
    scale_x_discrete(limits=spec.order)

########################################
# random effects model ANOVA analysis  #
########################################

library(lme4) # https://towardsdatascience.com/random-effects-in-linear-models-15845b2540ac

  # random effects can explain the amount of varaiance among species.  Phylogenetic structuring can influence
  # species level measurements and therefore attempts to account for this effect.

null.model <- lm(density ~ type, data=master) # null model

# model.1 <- lmer(density ~ type + (1|name), data=master) # random effects model

model.1b <- aov(density ~ type + name, data=master) # random effects model

anova(null.model)
  #Call:
    #aov(formula = null.model)
    #
    #Terms:
    #type Residuals
    #Sum of Squares  30.162786  6.035581
    #Deg. of Freedom         1       218
    #
    #Residual standard error: 0.1663916
    #Estimated effects may be unbalanced

anova(model.1b)
    #Analysis of Variance Table

    #Response: density
    #           Df  Sum Sq Mean Sq  F value    Pr(>F)
    #type        1 30.1628 30.1628 1651.885 < 2.2e-16 ***
    #name       10  2.2376  0.2238   12.254 < 2.2e-16 ***
    #Residuals 208  3.7980  0.0183
    #---
    #Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
