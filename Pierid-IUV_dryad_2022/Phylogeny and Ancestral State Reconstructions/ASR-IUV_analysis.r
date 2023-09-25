### ASR of IUV in Pieridae
### 'C:\Program Files\R\R-4.1.0\bin'


### Libraries
library('ape');library('geiger');library('phytools');library('corHMM')

### Data Files
scales <- read.csv('D:/Pierid_phylogenetics/uv_mat_syn-edit-V2_2-22-22.csv')
tree <- ladderize(read.tree('D:/Pierid_Phylogenetics/Pieridae_7114_final_names_trim_dates_2-22-2022.tre'))

##########################
### Data Wrangle         #
##########################
# Data format
uv <- scales[,4:6];rownames(uv) <- scales[,3]           #male/female presence/absence
uv.data <- treedata(tree, uv, sort = T, warnings = T)   #match tree and data
phy <- uv.data$phy                                      #trimmed tree
data <- uv.data$data                                    #trimmed data

# Final data
a <- scales[match(phy$tip.label,scales[,3]),c(3,4)]     #male
b <- scales[match(phy$tip.label,scales[,3]),c(3,5)]     #female
c <- scales[match(phy$tip.label,scales[,3]),c(3,6)]     #multistate
d <- data.frame(phy$tip.label,male=a[,2],female=b[,2])  #binary

# verify position of 3rd node for fixing
testclade <- extract.clade(phy, length(phy$tip.label) + 3)
plot(testclade) # first split is between Coliadinae and Pierinae

### corHMM

# Note that the rate classes are ordered from slowest (R1) to fastest (Rn) with respect to state 0
# 1,R1 = slow absent
# 2,R1 = slow present
# 1,R2 = fast absent
# 2,R2 = fast present
# Transitions between absent and present in the slow category in either direction are 100.

############################
# Presence/Absence Dataset #
############################

male.ER <- corHMM(phy,a,rate.cat=1,model="ER")
male.ARD <- corHMM(phy,a,rate.cat=1)
male.ER.2cat <- corHMM(phy,a,rate.cat=2,model="ER",get.tip.states = T)
male.ARD.2cat <- corHMM(phy,a,rate.cat=2)

# Appending node labels for fixing 3rd node as present
phy$node.label <- rep(times=phy$Nnode, NA)
phy$node.label[[3]] <- 2 # 1 is absence, 2 is presence

male.ER.fix2 <- corHMM(phy,a,rate.cat=1,model="ER", fixed.nodes=T)
male.ARD.fix2 <- corHMM(phy,a,rate.cat=1,fixed.nodes=T)
male.ER.2cat.fix2 <- corHMM(phy,a,rate.cat=2,model='ER',get.tip.states=T,fixed.nodes=T)
male.ARD.2cat.fix2 <- corHMM(phy,a,rate.cat=2,get.tip.states=T,fixed.nodes=T)

# Appending node labels for fixing 3rd node as absent
phy$node.label[[3]] <- 1 # 1 is absence, 2 is presence

male.ER.fix1 <- corHMM(phy,a,rate.cat=1,model="ER", fixed.nodes=T)
male.ARD.fix1 <- corHMM(phy,a,rate.cat=1,fixed.nodes=T)
male.ER.2cat.fix1 <- corHMM(phy,a,rate.cat=2,model='ER',get.tip.states=T,fixed.nodes=T)
male.ARD.2cat.fix1 <- corHMM(phy,a,rate.cat=2,get.tip.states=T,fixed.nodes=T)

male.scores <- data.frame(rbind(male.ER[1:3],
                                male.ARD[1:3],
                                male.ER.2cat[1:3],
                                male.ARD.2cat[1:3],
                                male.ER.2cat.fix1[1:3],
                                male.ARD.2cat.fix1[1:3],
                                male.ER.fix1[1:3],
                                male.ARD.fix1[1:3],
                                male.ER.2cat.fix2[1:3],
                                male.ARD.2cat.fix2[1:3],
                                male.ER.fix2[1:3],
                                male.ARD.fix2[1:3]))
rownames(male.scores) <- c("ER","ARD","ER2","ARD2","ER2F1","ARD2F1","ERF1","ARDF1","ER2F2","ARD2F2","ERF2","ARDF2")
male.scores$AICw <- round(exp(0.5*(min(unlist(male.scores$AICc))-unlist(male.scores$AICc)))/
                            sum(exp(0.5*(min(unlist(male.scores$AICc))-unlist(male.scores$AICc)))),2)
male.scores

#################################
# Mutual Ornamentation Dataset  #
#################################

# Multistate
multi.ER <- corHMM(phy,c,rate.cat=1,model = "ER", nstarts = 10)
multi.SYM <- corHMM(phy,c,rate.cat=1,model = "SYM", nstarts = 10)
multi.ARD <- corHMM(phy,c,rate.cat=1, nstarts = 10)
binary.ER <- corHMM(phy,d,rate.cat=1,model = "ER", nstarts = 10)
binary.SYM <- corHMM(phy,d,rate.cat=1,model="SYM", nstarts = 10)
binary.ARD <- corHMM(phy,d,rate.cat=1, nstarts = 10)

multi.ER.2state <- corHMM(phy,c,rate.cat=2, model="ER", nstarts = 10)
multi.SYM.2state <- corHMM(phy,c,rate.cat=2, model="SYM", nstarts = 10)
multi.ARD.2state <- corHMM(phy,c,rate.cat=2, nstarts = 10)
binary.ER.2state <- corHMM(phy,d,rate.cat=2, model="ER", nstarts = 10)
binary.SYM.2state <- corHMM(phy,d,rate.cat=2, model="SYM", nstarts = 10)
binary.ARD.2state <- corHMM(phy,d,rate.cat=2, nstarts = 10)

# Appending node labels for fixing 3rd node as absent
phy$node.label[[3]] <- 1 # 1 is absence, 2 is males, 3 is both

multi.ER.fix1 <- corHMM(phy,c,rate.cat=1,model = "ER", nstarts = 10, fixed.nodes=T)
multi.SYM.fix1 <- corHMM(phy,c,rate.cat=1,model = "SYM", nstarts = 10, fixed.nodes=T)
multi.ARD.fix1 <- corHMM(phy,c,rate.cat=1,nstarts = 10,fixed.nodes=T)
binary.ER.fix1 <- corHMM(phy,d,rate.cat=1,model = "ER", nstarts = 10,fixed.nodes=T)
binary.SYM.fix1 <- corHMM(phy,d,rate.cat=1,model="SYM", nstarts = 10,fixed.node=T)
binary.ARD.fix1 <- corHMM(phy,d,rate.cat=1, nstarts = 10,fixed.node=T)

multi.ER.2state.fix1 <- corHMM(phy,c,rate.cat=2, model="ER", nstarts = 10,fixed.nodes=T)
multi.SYM.2state.fix1 <- corHMM(phy,c,rate.cat=2, model="SYM", nstarts = 10,fixed.nodes=T)
multi.ARD.2state.fix1 <- corHMM(phy,c,rate.cat=2, nstarts = 10,fixed.nodes=T)
binary.ER.2state.fix1 <- corHMM(phy,d,rate.cat=2, model="ER", nstarts = 10,fixed.nodes=T)
binary.SYM.2state.fix1 <- corHMM(phy,d,rate.cat=2, model="SYM", nstarts = 10,fixed.nodes=T)
binary.ARD.2state.fix1 <- corHMM(phy,d,rate.cat=2, nstarts = 10,fixed.nodes=T)

# Appending node labels for fixing 3rd node as absent
phy$node.label[[3]] <- 2 # 1 is absence, 2 is males, 3 is both

multi.ER.fix2 <- corHMM(phy,c,rate.cat=1,model = "ER", nstarts = 10, fixed.nodes=T)
multi.SYM.fix2 <- corHMM(phy,c,rate.cat=1,model = "SYM", nstarts = 10, fixed.nodes=T)
multi.ARD.fix2 <- corHMM(phy,c,rate.cat=1,nstarts = 10,fixed.nodes=T)
binary.ER.fix2 <- corHMM(phy,d,rate.cat=1,model = "ER", nstarts = 10,fixed.nodes=T)
binary.SYM.fix2 <- corHMM(phy,d,rate.cat=1,model="SYM", nstarts = 10,fixed.node=T)
binary.ARD.fix2 <- corHMM(phy,d,rate.cat=1, nstarts = 10,fixed.node=T)

multi.ER.2state.fix2 <- corHMM(phy,c,rate.cat=1,model="ER", nstarts = 10,fixed.node=T)
multi.SYM.2state.fix2 <- corHMM(phy,c,rate.cat=1,model="SYM", nstarts = 10,fixed.node=T)
multi.ARD.2state.fix2 <- corHMM(phy,c,rate.cat=1, nstarts = 10,fixed.node=T)
binary.ER.2state.fix2 <- corHMM(phy,d,rate.cat=2, model="ER", nstarts = 10, fixed.node=T)
binary.SYM.2state.fix2 <- corHMM(phy,d,rate.cat=2, model="SYM", nstarts = 10, fixed.node=T)
binary.ARD.2state.fix2 <- corHMM(phy,d,rate.cat=2, nstarts = 10, fixed.node=T)

# Appending node labels for fixing 3rd node as absent
phy$node.label[[3]] <- 3 # 1 is absence, 2 is males, 3 is both

multi.ER.fix3 <- corHMM(phy,c,rate.cat=1,model = "ER", nstarts = 10, fixed.nodes=T)
multi.SYM.fix3 <- corHMM(phy,c,rate.cat=1,model = "SYM", nstarts = 10, fixed.nodes=T)
multi.ARD.fix3 <- corHMM(phy,c,rate.cat=1,nstarts = 10,fixed.nodes=T)
binary.ER.fix3 <- corHMM(phy,d,rate.cat=1,model = "ER", nstarts = 10,fixed.nodes=T)
binary.SYM.fix3 <- corHMM(phy,d,rate.cat=1,model="SYM", nstarts = 10,fixed.node=T)
binary.ARD.fix3 <- corHMM(phy,d,rate.cat=1, nstarts = 10,fixed.node=T)

multi.ER.2state.fix3 <- corHMM(phy,c,rate.cat=1,model="ER", nstarts = 10,fixed.node=T)
multi.SYM.2state.fix3 <- corHMM(phy,c,rate.cat=1,model="SYM", nstarts = 10,fixed.node=T)
multi.ARD.2state.fix3 <- corHMM(phy,c,rate.cat=1, nstarts = 10,fixed.node=T)
binary.ER.2state.fix3 <- corHMM(phy,d,rate.cat=2, model="ER", nstarts = 10, fixed.node=T)
binary.SYM.2state.fix3 <- corHMM(phy,d,rate.cat=2, model="SYM", nstarts = 10, fixed.node=T)
binary.ARD.2state.fix3 <- corHMM(phy,d,rate.cat=2, nstarts = 10, fixed.node=T)

multi.scores <- data.frame(rbind(multi.ER[1:3],
                                 multi.SYM[1:3],
                                 multi.ARD[1:3],
                                 binary.ER[1:3],
                                 binary.SYM[1:3],
                                 binary.ARD[1:3],
                                 multi.ER.2state[1:3],
                                 multi.SYM.2state[1:3],
                                 multi.ARD.2state[1:3],
                                 binary.ER.2state[1:3],
                                 binary.SYM.2state[1:3],
                                 binary.ARD.2state[1:3],

                                 multi.ER.fix1[1:3],
                                 multi.SYM.fix1[1:3],
                                 multi.ARD.fix1[1:3],
                                 binary.ER.fix1[1:2],
                                 binary.SYM.fix1[1:3],
                                 binary.ARD.fix1[1:3],
                                 multi.ER.2state.fix1[1:3],
                                 multi.SYM.2state.fix1[1:3],
                                 multi.ARD.2state.fix1[1:3],
                                 binary.ER.2state.fix1[1:3],
                                 binary.SYM.2state.fix1[1:3],
                                 binary.ARD.2state.fix1[1:3],

                                 multi.ER.fix2[1:3],
                                 multi.SYM.fix2[1:3],
                                 multi.ARD.fix2[1:3],
                                 binary.ER.fix2[1:2],
                                 binary.SYM.fix2[1:3],
                                 binary.ARD.fix2[1:3],
                                 multi.ER.2state.fix2[1:3],
                                 multi.SYM.2state.fix2[1:3],
                                 multi.ARD.2state.fix2[1:3],
                                 binary.ER.2state.fix2[1:3],
                                 binary.SYM.2state.fix2[1:3],
                                 binary.ARD.2state.fix2[1:3],

                                 multi.ER.fix3[1:3],
                                 multi.SYM.fix3[1:3],
                                 multi.ARD.fix3[1:3],
                                 binary.ER.fix3[1:2],
                                 binary.SYM.fix3[1:3],
                                 binary.ARD.fix3[1:3],
                                 multi.ER.2state.fix3[1:3],
                                 multi.SYM.2state.fix3[1:3],
                                 multi.ARD.2state.fix3[1:3],
                                 binary.ER.2state.fix3[1:3],
                                 binary.SYM.2state.fix3[1:3],
                                 binary.ARD.2state.fix3[1:3]))

rownames(multi.scores) <- c("multi.ER",
                            "multi.SYM",
                            "multi.ARD",
                            "binary.ER",
                            "binary.SYM",
                            "binary.ARD",
                            "multi.ER2",
                            "multi.SMY2",
                            "multi.ARD2",
                            "binary.ER2",
                            "binary.SYM2",
                            "binary.ARD2",

                            "multi.ER.F1",
                            "multi.SYM.F1",
                            "multi.ARD.F1",
                            "binary.ER.F1",
                            "binary.SYM.F1",
                            "binary.ARD.F1",
                            "multi.ER2.F1",
                            "multi.SYM2.F1",
                            "multi.ARD2.F1",
                            "binary.ER2.F1",
                            "binary.SYM2.F1",
                            "binary.ARD2.F1",

                            "multi.ER.F2",
                            "multi.SYM.F2",
                            "multi.ARD.F2",
                            "binary.ER.F2",
                            "binary.SYM.F2",
                            "binary.ARD.F2",
                            "multi.ER2.F2",
                            "multi.SYM2.F2",
                            "multi.ARD2.F2",
                            "binary.ER2.F2",
                            "binary.SYM2.F2",
                            "binary.ARD2.F2",

                            "multi.ER.F3",
                            "multi.SYM.F3",
                            "multi.ARD.F3",
                            "binary.ER.F3",
                            "binary.SYM.F3",
                            "binary.ARD.F3",
                            "multi.ER2.F3",
                            "multi.SYM2.F3",
                            "multi.ARD2.F3",
                            "binary.ER2.F3",
                            "binary.SYM2.F3",
                            "binary.ARD2.F3")

multi.scores$AICw <- round(exp(0.5*(min(unlist(multi.scores$AICc))-unlist(multi.scores$AICc)))/
                             sum(exp(0.5*(min(unlist(multi.scores$AICc))-unlist(multi.scores$AICc)))),2)
multi.scores
