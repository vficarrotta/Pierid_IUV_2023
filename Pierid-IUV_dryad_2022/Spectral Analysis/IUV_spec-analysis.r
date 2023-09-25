### Spectral Analysis of Iridescent Ultraviolet Wing Patterns
# 'D:/Pierid_phylogenetics/Spec Data/IUV raw text files'

### Libraries
library('pavo'); library('ggplot2')

### Data
data <- getspec('D:/Pierid_phylogenetics/Spec Data/IUV raw text files', ext='txt')
    # test plot
    plot(data[,2]~data[,1])

wlx <- data$wl      # vector of x axis wavelengths
data1 <- data[-1]   # data frame of measure only

### splitting trials by species
zfau <- data1[1:9]
cdan <- data1[10:18]
gmah <- data1[19:27]
ipyr <- data1[28:36]
dver <- data1[37:45]
hgla <- data1[46:54]
psen <- data1[55:63]
ceur <- data1[64:72]
aled <- data1[73:81]
enis <- data1[82:90]
zces <- data1[90:99]

### Average all trials/replicates per species
zfau_a <- apply(zfau,1,FUN=mean)
cdan_a <- apply(cdan,1,FUN=mean)
gmah_a <- apply(gmah,1,FUN=mean)
ipyr_a <- apply(ipyr,1,FUN=mean)
dver_a <- apply(dver,1,FUN=mean)
hgla_a <- apply(hgla,1,FUN=mean)
psen_a <- apply(psen,1,FUN=mean)
ceur_a <- apply(ceur,1,FUN=mean)
aled_a <- apply(aled,1,FUN=mean)
enis_a <- apply(enis,1,FUN=mean)
zces_a <- apply(zces,1,FUN=mean)

### compiled averages
data2 <- data.frame(zfau_a, cdan_a, gmah_a, ipyr_a, dver_a, hgla_a, psen_a, ceur_a, aled_a, enis_a, zces_a)

### normalize readings
data3 <- procspec(data2, op='min'); data3[1] <- wlx # minimum normalization
  # test
  plot(data3[,3]~data3[,1], type='l')

### Plots of spectra per species
setwd('D:/Pierid_Phylogenetics/Spec Data/plots/3.5x3.25 plots')

ggplot() +
  geom_line(data3, mapping=aes(x=wl, y=zfau_a)) +
  theme_classic() +
  scale_y_continuous(limits=c(0,100)) +
  theme(legend.position='none') +
  ggtitle('Z. fausti') +
  xlab('Wavelength (nm)') +
  ylab('Percent Reflectance')
ggsave('Z.fausti_norm.pdf', device = 'pdf', dpi = 'retina', height=3.5, width=3.25, units='in')

ggplot() +
  geom_line(data3, mapping=aes(x=data3$wl, y=data3$cdan_a)) +
  theme_classic() +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,100)) +
  ggtitle('C. danae') +
  xlab('Wavelength (nm)') +
  ylab('Percent Reflectance')
ggsave('C.danae_norm.pdf', device = 'pdf', dpi = 'retina', height=3.5, width=3.25, units='in')

ggplot() +
  geom_line(data3, mapping=aes(x=data3$wl, y=data3$gmah_a)) +
  theme_classic() +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,100)) +
  ggtitle('G. mahaguru') +
  xlab('Wavelength (nm)') +
  ylab('Percent Reflectance')
ggsave('G.mahaguru_norm.pdf', device = 'pdf', dpi = 'retina', height=3.5, width=3.25, units='in')

ggplot() +
  geom_line(data3, mapping=aes(x=data3$wl, y=data3$ipyr_a)) +
  theme_classic() +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,100)) +
  ggtitle('I. pyrene') +
  xlab('Wavelength (nm)') +
  ylab('Percent Reflectance')
ggsave('I.pyrene_norm.pdf', device = 'pdf', dpi = 'retina', height=3.5, width=3.25, units='in')

ggplot() +
  geom_line(data3, mapping=aes(x=data3$wl, y=data3$dver_a)) +
  theme_classic() +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,100)) +
  ggtitle('D. verhuelli') +
  xlab('Wavelength (nm)') +
  ylab('Percent Reflectance')
ggsave('D.verhuelli_norm.pdf', device = 'pdf', dpi = 'retina', height=3.5, width=3.25, units='in')

ggplot() +
  geom_line(data3, mapping=aes(x=data3$wl, y=data3$hgla_a)) +
  theme_classic() +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,100)) +
  ggtitle('H. glaucippe') +
  xlab('Wavelength (nm)') +
  ylab('Percent Reflectance')
ggsave('H.glaucippe_norm.pdf', device = 'pdf', dpi = 'retina', height=3.5, width=3.25, units='in')

ggplot() +
  geom_line(data3, mapping=aes(x=data3$wl, y=data3$psen_a)) +
  theme_classic() +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,100)) +
  ggtitle('P. sennae') +
  xlab('Wavelength (nm)') +
  ylab('Percent Reflectance')
ggsave('P.sennae_norm.pdf', device = 'pdf', dpi = 'retina', height=3.5, width=3.25, units='in')

ggplot() +
  geom_line(data3, mapping=aes(x=data3$wl, y=data3$ceur_a)) +
  theme_classic() +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,100)) +
  ggtitle('C. eurytheme') +
  xlab('Wavelength (nm)') +
  ylab('Percent Reflectance')
ggsave('C.eurytheme_norm.pdf', device = 'pdf', dpi = 'retina', height=3.5, width=3.25, units='in')

ggplot() +
  geom_line(data3, mapping=aes(x=data3$wl, y=data3$aled_a)) +
  theme_classic() +
  theme(legend.position='none') +
  ggtitle('A. leda') +
  scale_y_continuous(limits=c(0,100)) +
  xlab('Wavelength (nm)') +
  ylab('Percent Reflectance')
ggsave('A.leda_norm.pdf', device = 'pdf', dpi = 'retina', height=3.5, width=3.25, units='in')

ggplot() +
  geom_line(data3, mapping=aes(x=data3$wl, y=data3$enis_a)) +
  theme_classic() +
  theme(legend.position='none') +
  ggtitle('E. nise') +
  scale_y_continuous(limits=c(0,100)) +
  xlab('Wavelength (nm)') +
  ylab('Percent Reflectance')
ggsave('E.nise_norm.pdf', device = 'pdf', dpi = 'retina', height=3.5, width=3.25, units='in')

ggplot() +
  geom_line(data3, mapping=aes(x=data3$wl, y=data3$zces_a)) +
  theme_classic() +
  theme(legend.position='none') +
  scale_y_continuous(limits=c(0,100)) +
  ggtitle('Z. cesonia') +
  xlab('Wavelength (nm)') +
  ylab('Percent Reflectance')
ggsave('Z.cesonia_norm.pdf', device = 'pdf', dpi = 'retina', height=3.5, width=3.25, units='in')
