# IUV scoring QC
### 'C:\Program Files\R\R-4.1.0\bin'

###################################################################################
# Checking state distribution ratios between full data sets and working data sets #
###################################################################################

library('ggplot2')

### Data
uvg <- read.csv('D:/Pierid_Phylogenetics/uvg_7-7-2022.csv')

### Plot - point sizes are frequency of uv+ species of genus
ggplot(uvg, aes(x = Genus, y = factor(number_species), color = matrix))+
    ylab('Species Count')+
    geom_point(alpha=.5, aes(size=Ratio))+
    geom_point(shape=1, color='black', aes(size=Ratio))+
    theme_classic()+
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.1))+
    scale_color_manual('', values = c('cyan2', 'yellow2'))+
    labs(title = 'Frequency of UV+ Species for Genera', size = 'UV+ Species/Species Count')

################################################################
# Chi-squared test of state representation in reduced dataset  #
################################################################

  # full dataset
scales <- read.csv('D:/Pierid_phylogenetics/uv_mat_syn-edit-V2_2-22-22.csv')
  # genera with IUV from full dataset
genus_uv <- read.csv('D:/Pierid_Phylogenetics/genus_uv_7-7-2022.csv')
  # species overlap between tree and IUV dataset
overlap_data <- read.csv('D:/Pierid_Phylogenetics/data-tree_overlap_7-7-2022.csv')
  # IUV genera overlap between tree adn IUV dataset
genus_uv_overlap <- read.csv('D:/Pierid_Phylogenetics/genus_uv_overlap_7-7-2022.csv')

   # taxa in dataset with IUV
d_uvp <- scales[which(scales$male_ultrastructure_dorsal == 1),] # 0.1755102
   # taxa in dataset without IUV
d_uva <- scales[which(scales$male_ultrastructure_dorsal == 0),] # 0.8244898

gwuvp <- unique(d_uvp$Genus)

genus_uv <- scales[scales$Genus %in% gwuvp, ]


# Full dataset percentages of IUV and non-IUV
fuvpr <- length(d_uvp[,1]) / (length(d_uvp[,1]) + length(d_uva[,1])) # 0.1755102
fuvar <- 1 - (length(d_uvp[,1]) / (length(d_uvp[,1]) + length(d_uva[,1]))) # 0.8244898

# Overlap dataset
  # taxa in dataset with uv
o_uvp <- overlap_data[overlap_data$IUV == 1,]
  # taxa in dataset without uv
o_vap <- overlap_data[overlap_data$IUV == 0,]

ouvpr <- length(o_uvp[,1]) / (length(o_uvp[,1]) + length(o_vap[,1])) # 0.1908602
ouvar <- 1 - (length(o_uvp[,1]) / (length(o_uvp[,1]) + length(o_vap[,1]))) # 0.8091398

# Generate expected full IUV data
data_chi <- as.data.frame(rbind(c(
  # counts expected ovelap
        round(fuvpr * (length(o_uvp[,1]) + length(o_vap[,1]))), # 62
        round(fuvar * (length(o_uvp[,1]) + length(o_vap[,1]))) # 307
    ),
  # counts observed overlap
  c(
        length(o_uvp[,1]), # 69
        length(o_vap[,1]) # 300
    )))

colnames(data_chi) <- c('IUV', 'species_count')

chisq.test(data_chi)

      # Pearson's Chi-squared test with Yates' continuity correction

      # data:  data_chi
      # X-squared = 0.1428, df = 1, p-value = 0.7055

### Chi-squared on genus with UV
sfiuvg <- split(genus_uv, f=genus_uv$Genus)
soiuvg <- split(genus_uv_overlap, f=genus_uv_overlap$Genus)

counts_sfiuvg <- unlist(lapply(sfiuvg, function(x) length(x$IUV)))
counts_soiuvg <- unlist(lapply(soiuvg, function(x) length(x$IUV)))

sums_sfiuvg <- unlist(lapply(sfiuvg, function(x) sum(x$IUV)))
sums_soiuvg <- unlist(lapply(soiuvg, function(x) sum(x$IUV)))

exp_n_IUV <- round((sums_sfiuvg/counts_sfiuvg) * counts_soiuvg) # expected number of IUV in each UV+ genera

o_data_chi <- cbind(exp_n_IUV, counts_soiuvg)

chisq.test(o_data_chi)

  #       Pearson's Chi-squared test

  # data:  o_data_chi
  # X-squared = 2.9945, df = 16, p-value = 0.9998
