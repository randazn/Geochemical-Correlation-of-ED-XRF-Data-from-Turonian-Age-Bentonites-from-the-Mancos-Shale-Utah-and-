
#Call data from Excel file
library(readxl)
Bentonite <- read_excel("C:/Mancos Bentonite Data.xlsx")

#Limit the number of columns and rows to those containing the desried information
BentonitePCA <- Bentonite[1:243,1:58]

#Calculate the mean of each of the columns containing selected elemental data and round to 2 digits
round_mean_Fe <- round(mean(BentonitePCA$Fe), digits =2)
round_mean_Th <- round(mean(BentonitePCA$Th), digits =2)
round_mean_Ni <- round(mean(BentonitePCA$Ni), digits =2)
round_mean_Mn <- round(mean(BentonitePCA$Mn), digits =2)
round_mean_Nd <- round(mean(BentonitePCA$Nd), digits =2)
round_mean_Zr <- round(mean(BentonitePCA$Zr), digits =2)
round_mean_Y <- round(mean(BentonitePCA$Y), digits =2)
round_mean_S <- round(mean(BentonitePCA$S), digits =2)
round_mean_Ca <- round(mean(BentonitePCA$Ca), digits =2)
round_mean_Ce <- round(mean(BentonitePCA$Ce), digits =2)
round_mean_La <- round(mean(BentonitePCA$La), digits =2)

#Ground data based upon where the sampels were collected
SST <- BentonitePCA[1:23,1:58]

Shale <- BentonitePCA[28:52,1:58]

B6 <- BentonitePCA[53:77,1:58]

UpperGallup <- BentonitePCA[78:108,1:58]

FBase1 <- BentonitePCA[109:134,1:58]

FBase2 <- BentonitePCA[135:163,1:58]

FBase3 <- BentonitePCA[164:191,1:58]

EBase <- BentonitePCA[192:217,1:58]

Torrivio <- BentonitePCA[218:243,1:58]

#Prepare the data for principal component analysis
#Sandstone data
SST1<- SST[,apply(SST, 2, var, na.rm=TRUE) != 0]
SST.pr<- prcomp(SST1[c(1:20)], center = TRUE, scale = TRUE)
s <- summary(SST.pr)
biplot(SST.pr)
# Create groups
pch.group <- c(rep(21, times=8), rep(22, times=10), rep(24, times=30))
col.group <- c(rep("skyblue2", times=8), rep("gold", times=10), rep("green2", times=30))

#Plot individuals
plot(SST.pr$x[,1], SST.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)

#Shale data
Shale1<- Shale[,apply(Shale, 2, var, na.rm=TRUE) != 0]
Shale.pr<- prcomp(Shale1[c(1:54)], center = TRUE, scale = TRUE)
summary(Shale.pr)

#Juana Lopez data
B61<- B6[,apply(B6, 2, var, na.rm=TRUE) != 0]
B6.pr<- prcomp(B61[c(1:57)], center = TRUE, scale = TRUE)
summary(B6.pr)

#Upper Gallup Data
UpperGallup1<- UpperGallup[,apply(UpperGallup, 2, var, na.rm=TRUE) != 0]
UpperGallup.pr<- prcomp(UpperGallup1[c(1:58)], center = TRUE, scale = TRUE)
summary(UpperGallup.pr)

#LS4 1 Data
FBase1_1<- FBase1[,apply(FBase1, 2, var, na.rm=TRUE) != 0]
FBase1.pr<- prcomp(FBase1_1[c(1:57)], center = TRUE, scale = TRUE)
summary(FBase1.pr)

#LS4 2 Data
FBase2_1<- FBase2[,apply(FBase2, 2, var, na.rm=TRUE) != 0]
FBase2.pr<- prcomp(FBase2_1[c(1:57)], center = TRUE, scale = TRUE)
summary(FBase2.pr)

#LS4 3 Data
FBase3_1<- FBase3[,apply(FBase3, 2, var, na.rm=TRUE) != 0]
FBase3.pr<- prcomp(FBase3_1[c(1:56)], center = TRUE, scale = TRUE)
summary(FBase3.pr)

#LS5 1 Data
EBase1<- EBase[,apply(EBase, 2, var, na.rm=TRUE) != 0]
EBase.pr<- prcomp(EBase1[c(1:57)], center = TRUE, scale = TRUE)
summary(EBase.pr)

#Torrivio Data
Torrivio1<- Torrivio[,apply(Torrivio, 2, var, na.rm=TRUE) != 0]
Torrivio.pr<- prcomp(Torrivio1[c(1:55)], center = TRUE, scale = TRUE)
summary(Torrivio.pr)

wdbc.pr <- prcomp(BentonitePCA[c(2:58)], center = TRUE, scale = TRUE)
summary(wdbc.pr)


pairwise.comparison(wdbc.pr,group,members=NULL,spots=NULL,a.order=NULL,b.order=NULL,method="unlogged",logged=TRUE)


install.packages("psych")
library(psych)
library(readxl)

#Calculate the correlation matrices for first 31 elements

pairs.panels(SST[2:10])
pairs.panels(SST[11:21])
pairs.panels(SST[22:32])

#Subset the data based on redox and terrigenous elements and then calculate the correlation matrices for redox and terrigenous
SSTSubset <- SST[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Ce","La")]
pairs.panels(SSTSubset[1:11])
SSTTerra <- SST[,c("Al","K","Si","Ca","Ti","Cl","S","Zr","Fe")]
pairs.panels(SSTTerra[1:9])

#Calculate the rounded mean for Fe, Th, Ni, Mn, Nd, Zr, Y, S, Ca, Ce, and La
round_mean_SSTFe <- round(mean(SSTSubset$Fe), digits =2)
round_mean_SSTTh <- round(mean(SSTSubset$Th), digits =2)
round_mean_SSTNi <- round(mean(SSTSubset$Ni), digits =2)
round_mean_SSTMn <- round(mean(SSTSubset$Mn), digits =2)
round_mean_SSTNd <- round(mean(SSTSubset$Nd), digits =2)
round_mean_SSTZr <- round(mean(SSTSubset$Zr), digits =2)
round_mean_SSTY <- round(mean(SSTSubset$Y), digits =2)
round_mean_SSTS <- round(mean(SSTSubset$S), digits =2)
round_mean_SSTCa <- round(mean(SSTSubset$Ca), digits =2)
round_mean_SSTCe <- round(mean(SSTSubset$Ce), digits =2)
round_mean_SSTLa <- round(mean(SSTSubset$La), digits =2)

cbindFe <- cbind(round_mean_EBaseCa,round_mean_FBase1Ca)

#Calculate the correlation matrices for the first 31 elements
pairs.panels(Shale[2:10])
pairs.panels(Shale[11:21])
pairs.panels(Shale[22:32]) 

#Subset the data based on redox and terrigenous elements and then calculate the correlation matrices for redox and terrigenous
ShaleSubset <- Shale[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Ce", "La")]
pairs.panels(ShaleSubset[1:11])
ShaleTerra <- Shale[,c("Al","K","Si","Ca","Ti","Cl","S","Zr","Fe")]
pairs.panels(ShaleTerra[1:9])

#Calculate the rounded mean for Fe, Th, Ni, Mn, Nd, Zr, Y, S, Ca, Ce, and La
round_mean_ShaleFe <- round(mean(ShaleSubset$Fe), digits =2)
round_mean_ShaleTh <- round(mean(ShaleSubset$Th), digits =2)
round_mean_ShaleNi <- round(mean(ShaleSubset$Ni), digits =2)
round_mean_ShaleMn <- round(mean(ShaleSubset$Mn), digits =2)
round_mean_ShaleNd <- round(mean(ShaleSubset$Nd), digits =2)
round_mean_ShaleZr <- round(mean(ShaleSubset$Zr), digits =2)
round_mean_ShaleY <- round(mean(ShaleSubset$Y), digits =2)
round_mean_ShaleS <- round(mean(ShaleSubset$S), digits =2)
round_mean_ShaleCa <- round(mean(ShaleSubset$Ca), digits =2)
round_mean_ShaleCe <- round(mean(ShaleSubset$Ce), digits =2)
round_mean_ShaleLa <- round(mean(ShaleSubset$La), digits =2)


#Calculate the correlation matrices for all of the elemental data
pairs.panels(B6[2:10])
pairs.panels(B6[11:21])
pairs.panels(B6[22:32])
pairs.panels(B6[33:43])
pairs.panels(B6[44:54])
pairs.panels(B6[55:58])

#Subset the data based on redox and terrigenous elements and then calculate the correlation matrices for redox and terrigenous
B6Subset <- B6[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Ce","La")]
pairs.panels(B6Subset[1:11])
B6Terra <- B6[,c("Al","K","Si","Ca","Ti","Cl","S","Zr","Fe")]
pairs.panels(B6Terra[1:9])

#Calculate the rounded mean for Fe, Th, Ni, Mn, Nd, Zr, Y, S, Ca, Ce, and La
round_mean_B6Fe <- round(mean(B6Subset$Fe), digits =2)
round_mean_B6Th <- round(mean(B6Subset$Th), digits =2)
round_mean_B6Ni <- round(mean(B6Subset$Ni), digits =2)
round_mean_B6Mn <- round(mean(B6Subset$Mn), digits =2)
round_mean_B6Nd <- round(mean(B6Subset$Nd), digits =2)
round_mean_B6Zr <- round(mean(B6Subset$Zr), digits =2)
round_mean_B6Y <- round(mean(B6Subset$Y), digits =2)
round_mean_B6S <- round(mean(B6Subset$S), digits =2)
round_mean_B6Ca <- round(mean(B6Subset$Ca), digits =2)
round_mean_B6Ce <- round(mean(B6Subset$Ce), digits =2)
round_mean_B6La <- round(mean(B6Subset$La), digits =2)

#Calculate correlation matrices for all of the 57 analyzed elements
pairs.panels(UpperGallup[2:10])
pairs.panels(UpperGallup[11:21])
pairs.panels(UpperGallup[22:32])
pairs.panels(UpperGallup[33:43])
pairs.panels(UpperGallup[44:54])
pairs.panels(UpperGallup[55:58])

#Subset the data based on redox and terrigenous elements and then calculate the correlation matrices for redox and terrigenous
UpperGallupSubset <- UpperGallup[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Ce","La")]
pairs.panels(UpperGallupSubset[1:11])
UpperGallupTerra <- UpperGallup[,c("Al","K","Si","Ca","Ti","Cl","S","Zr","Fe")]
pairs.panels(UpperGallupTerra[1:9])

#Calculate the rounded mean for Fe, Th, Ni, Mn, Nd, Zr, Y, S, Ca, Ce, and La
round_mean_UpperGallupFe <- round(mean(UpperGallupSubset$Fe), digits =2)
round_mean_UpperGallupTh <- round(mean(UpperGallupSubset$Th), digits =2)
round_mean_UpperGallupNi <- round(mean(UpperGallupSubset$Ni), digits =2)
round_mean_UpperGallupMn <- round(mean(UpperGallupSubset$Mn), digits =2)
round_mean_UpperGallupNd <- round(mean(UpperGallupSubset$Nd), digits =2)
round_mean_UpperGallupZr <- round(mean(UpperGallupSubset$Zr), digits =2)
round_mean_UpperGallupY <- round(mean(UpperGallupSubset$Y), digits =2)
round_mean_UpperGallupS <- round(mean(UpperGallupSubset$S), digits =2)
round_mean_UpperGallupCa <- round(mean(UpperGallupSubset$Ca), digits =2)
round_mean_UpperGallupCe <- round(mean(UpperGallupSubset$Ce), digits =2)
round_mean_UpperGallupLa <- round(mean(UpperGallupSubset$La), digits =2)

Redox_Elements <- cbind(UpperGallup$Fe,B6$Fe,Shale$Fe,SST$Fe,FBase1$Fe,FBase2$Fe,FBase3$Fe,EBase$Fe,Torrivio$Fe)

Fe <- cbind(SST$Fe,FBase1$Fe,EBase$Fe)

pairs.panel(Fe[1-3])

sapply(UpperGallup[c("Fe","Th","Ni", "Nd", "Zr", "Y", "S", "Ca")], cor, UpperGallup["Cr"])

sapply(FBase1[c("Fe","Th","Ni", "Nd", "Zr", "Y", "S", "Ca")], cor, FBase1["Cr"])

pairs.panels(UpperGallup["Fe","Th","Ni"])


#Calculate the correlation matrices for all of the elemental data
pairs.panels(FBase1[2:10])
pairs.panels(FBase1[11:21])
pairs.panels(FBase1[22:32])
pairs.panels(FBase1[33:43])
pairs.panels(FBase1[44:54])
pairs.panels(FBase1[55:58])

#Subset the data based on redox and terrigenous elements and then calculate the correlation matrices for redox and terrigenous
FBase1Subset <- FBase1[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Ce","La")]
pairs.panels(FBase1Subset[1:11])
FBase1Terra <- FBase1[,c("Al","K","Si","Ca","Ti","Cl","S","Zr","Fe")]
pairs.panels(FBase1Terra[1:9])

#Calculate the rounded mean for Fe, Th, Ni, Mn, Nd, Zr, Y, S, Ca, Ce, and La
round_mean_FBase1Fe <- round(mean(FBase1Subset$Fe), digits =2)
round_mean_FBase1Th <- round(mean(FBase1Subset$Th), digits =2)
round_mean_FBase1Ni <- round(mean(FBase1Subset$Ni), digits =2)
round_mean_FBase1Mn <- round(mean(FBase1Subset$Mn), digits =2)
round_mean_FBase1Nd <- round(mean(FBase1Subset$Nd), digits =2)
round_mean_FBase1Zr <- round(mean(FBase1Subset$Zr), digits =2)
round_mean_FBase1Y <- round(mean(FBase1Subset$Y), digits =2)
round_mean_FBase1S <- round(mean(FBase1Subset$S), digits =2)
round_mean_FBase1Ca <- round(mean(FBase1Subset$Ca), digits =2)
round_mean_FBase1Ce <- round(mean(FBase1Subset$Ce), digits =2)
round_mean_FBase1La <- round(mean(FBase1Subset$La), digits =2)

#Calculate the correlation matrices for all of the elemental data
pairs.panels(FBase2[2:10])
pairs.panels(FBase2[11:21])
pairs.panels(FBase2[22:32])
pairs.panels(FBase2[33:43])
pairs.panels(FBase2[44:54])
pairs.panels(FBase2[55:58])

#Subset the data based on redox and terrigenous elements and then calculate the correlation matrices for redox and terrigenous
FBase2Subset <- FBase2[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Ce","La")]
pairs.panels(FBase2Subset[1:11])
FBase2Terra <- FBase2[,c("Al","K","Si","Ca","Ti","Cl","S","Zr","Fe")]
pairs.panels(FBase2Terra[1:9])

#Calculate the rounded mean for Fe, Th, Ni, Mn, Nd, Zr, Y, S, Ca, Ce, and La
round_mean_FBase2Fe <- round(mean(FBase2Subset$Fe), digits =2)
round_mean_FBase2Th <- round(mean(FBase2Subset$Th), digits =2)
round_mean_FBase2Ni <- round(mean(FBase2Subset$Ni), digits =2)
round_mean_FBase2Mn <- round(mean(FBase2Subset$Mn), digits =2)
round_mean_FBase2Nd <- round(mean(FBase2Subset$Nd), digits =2)
round_mean_FBase2Zr <- round(mean(FBase2Subset$Zr), digits =2)
round_mean_FBase2Y <- round(mean(FBase2Subset$Y), digits =2)
round_mean_FBase2S <- round(mean(FBase2Subset$S), digits =2)
round_mean_FBase2Ca <- round(mean(FBase2Subset$Ca), digits =2)
round_mean_FBase2Ce <- round(mean(FBase2Subset$Ce), digits =2)
round_mean_FBase2La <- round(mean(FBase2Subset$La), digits =2)

pairs.panels(FBase3[2:10])
pairs.panels(FBase3[11:21])
pairs.panels(FBase3[22:32])
pairs.panels(FBase3[33:43])
pairs.panels(FBase3[43:53])
pairs.panels(FBase3[54:58])

#Subset the data based on redox and terrigenous elements and then calculate the correlation matrices for redox and terrigenous
FBase3Subset <- FBase3[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Ce","La")]
pairs.panels(FBase3Subset[1:11])
FBase3Terra <- FBase3[,c("Al","K","Si","Ca","Ti","Cl","S","Zr","Fe")]

pairs.panels(FBase3Terra[1:9])

#Calculate the rounded mean for Fe, Th, Ni, Mn, Nd, Zr, Y, S, Ca, Ce, and La
round_mean_FBase3Fe <- round(mean(FBase3Subset$Fe), digits =2)
round_mean_FBase3Th <- round(mean(FBase3Subset$Th), digits =2)
round_mean_FBase3Ni <- round(mean(FBase3Subset$Ni), digits =2)
round_mean_FBase3Mn <- round(mean(FBase3Subset$Mn), digits =2)
round_mean_FBase3Nd <- round(mean(FBase3Subset$Nd), digits =2)
round_mean_FBase3Zr <- round(mean(FBase3Subset$Zr), digits =2)
round_mean_FBase3Y <- round(mean(FBase3Subset$Y), digits =2)
round_mean_FBase3S <- round(mean(FBase3Subset$S), digits =2)
round_mean_FBase3Ca <- round(mean(FBase3Subset$Ca), digits =2)
round_mean_FBase3Ce <- round(mean(FBase3Subset$Ce), digits =2)
round_mean_FBase3La <- round(mean(FBase3Subset$La), digits =2)

pairs.panels(EBase[2:10])
pairs.panels(EBase[11:21])
pairs.panels(EBase[22:32])
pairs.panels(EBase[33:43])
pairs.panels(EBase[44:54])
pairs.panels(EBase[55:58])

#Subset the data based on redox and terrigenous elements and then calculate the correlation matrices for redox and terrigenous
EBaseSubset <- EBase[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Ce","La")]
EBaseTerra <- EBase[,c("Al","K","Si","Ca","Ti","Cl","S","Zr", "Fe")]
pairs.panels(EBaseSubset[1:11])
pairs.panels(EBaseTerra[1:9])

#Calculate the rounded mean for Fe, Th, Ni, Mn, Nd, Zr, Y, S, Ca, Ce, and La
round_mean_EBaseFe <- round(mean(EBaseSubset$Fe), digits =2)
round_mean_EBaseTh <- round(mean(EBaseSubset$Th), digits =2)
round_mean_EBaseNi <- round(mean(EBaseSubset$Ni), digits =2)
round_mean_EBaseMn <- round(mean(EBaseSubset$Mn), digits =2)
round_mean_EBaseNd <- round(mean(EBaseSubset$Nd), digits =2)
round_mean_EBaseZr <- round(mean(EBaseSubset$Zr), digits =2)
round_mean_EBaseY <- round(mean(EBaseSubset$Y), digits =2)
round_mean_EBaseS <- round(mean(EBaseSubset$S), digits =2)
round_mean_EBaseCa <- round(mean(EBaseSubset$Ca), digits =2)
round_mean_EBaseCe <- round(mean(EBaseSubset$Ce), digits =2)
round_mean_EBaseLa <- round(mean(EBaseSubset$La), digits =2)

pairs.panels(Torrivio[2:10])
pairs.panels(Torrivio[11:21])
pairs.panels(Torrivio[22:32])
pairs.panels(Torrivio[33:43])
pairs.panels(Torrivio[44:54])
pairs.panels(Torrivio[55:58])

#Remove row with NA
Torrivio <- Torrivio[-c(30),]
#Subset the data based on redox and terrigenous elements and then calculate the correlation matrices for redox and terrigenous
TorrivioSubset <- Torrivio[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Ce","La")]
TorrivioTerra <- Torrivio[,c("Al","K","Si","Ca","Ti","Cl","S","Zr","Fe")]
pairs.panels(TorrivioSubset[1:11])
pairs.panels(TorrivioTerra[1:9])

#Calculate the rounded mean for Fe, Th, Ni, Mn, Nd, Zr, Y, S, Ca, Ce, and La
round_mean_TorrivioFe <- round(mean(TorrivioSubset$Fe), digits =2)
round_mean_TorrivioTh <- round(mean(TorrivioSubset$Th), digits =2)
round_mean_TorrivioNi <- round(mean(TorrivioSubset$Ni), digits =2)
round_mean_TorrivioMn <- round(mean(TorrivioSubset$Mn), digits =2)
round_mean_TorrivioNd <- round(mean(TorrivioSubset$Nd), digits =2)
round_mean_TorrivioZr <- round(mean(TorrivioSubset$Zr), digits =2)
round_mean_TorrivioY <- round(mean(TorrivioSubset$Y), digits =2)
round_mean_TorrivioS <- round(mean(TorrivioSubset$S), digits =2)
round_mean_TorrivioCa <- round(mean(TorrivioSubset$Ca), digits =2)
round_mean_TorrivioCe <- round(mean(TorrivioSubset$Ce), digits =2)
round_mean_TorrivioLa <- round(mean(TorrivioSubset$La), digits =2)

#Calculate the variance in each of the elements for all of th bentonite sample groupings (Fe, Th, Ni, Mn, Nd, Zr, Y, Ce, S, Ca, La)
FeVar <- (27 * (round_mean_SSTFe-round_mean_Fe)^2 + 25 * (round_mean_ShaleFe-round_mean_Fe)^2 + 28* (round_mean_FBase3Fe-round_mean_Fe)^2
+ 29 * (round_mean_FBase2Fe-round_mean_Fe)^2 + 26 * (round_mean_FBase1Fe-round_mean_Fe)^2 + 31 * (round_mean_UpperGallupFe-round_mean_Fe)^2 + 26 * (round_mean_EBaseFe-round_mean_Fe)^2
+ 25 * (round_mean_B6Fe-round_mean_Fe)^2 + 30 * (round_mean_TorrivioFe-round_mean_Fe)^2)/11

ThVar <- (27 * (round_mean_SSTTh-round_mean_Th)^2 + 25 * (round_mean_ShaleTh-round_mean_Th)^2 + 31 * (round_mean_UpperGallupTh-round_mean_Th)^2 + 28* (round_mean_FBase3Th-round_mean_Th)^2
+ 29 * (round_mean_FBase2Th-round_mean_Th)^2 + 26 * (round_mean_FBase1Th-round_mean_Th)^2 + 26 * (round_mean_EBaseTh-round_mean_Th)^2
+ 25 * (round_mean_B6Th-round_mean_Th)^2 + 30 * (round_mean_TorrivioTh-round_mean_Th)^2)/11

NiVar <- (27 * (round_mean_SSTNi-round_mean_Ni)^2 + 31 * (round_mean_UpperGallupNi-round_mean_Ni)^2 + 25 * (round_mean_ShaleNi-round_mean_Ni)^2 + 28* (round_mean_FBase3Ni-round_mean_Ni)^2
+ 29 * (round_mean_FBase2Ni-round_mean_Ni)^2 + 26 * (round_mean_FBase1Ni-round_mean_Ni)^2 + 26 * (round_mean_EBaseNi-round_mean_Ni)^2
+ 25 * (round_mean_B6Ni-round_mean_Ni)^2 + 30 * (round_mean_TorrivioNi-round_mean_Ni)^2)/11

MnVar <- (27 * (round_mean_SSTMn-round_mean_Mn)^2 + 31 * (round_mean_UpperGallupMn-round_mean_Mn)^2 + 25 * (round_mean_ShaleMn-round_mean_Mn)^2 + 28* (round_mean_FBase3Mn-round_mean_Mn)^2
+ 29 * (round_mean_FBase2Mn-round_mean_Mn)^2 + 26 * (round_mean_FBase1Mn-round_mean_Mn)^2 + 26 * (round_mean_EBaseMn-round_mean_Mn)^2
+ 25 * (round_mean_B6Mn-round_mean_Mn)^2 + 30 * (round_mean_TorrivioMn-round_mean_Mn)^2)/11

NdVar <- (27 * (round_mean_SSTNd-round_mean_Nd)^2 + 31 * (round_mean_UpperGallupNd-round_mean_Nd)^2 + 25 * (round_mean_ShaleNd-round_mean_Nd)^2 + 28* (round_mean_FBase3Nd-round_mean_Nd)^2
+ 29 * (round_mean_FBase2Nd-round_mean_Nd)^2 + 26 * (round_mean_FBase1Nd-round_mean_Nd)^2 + 26 * (round_mean_EBaseNd-round_mean_Nd)^2
+ 25 * (round_mean_B6Nd-round_mean_Nd)^2 + 30 * (round_mean_TorrivioNd-round_mean_Nd)^2)/11

ZrVar <- (27 * (round_mean_SSTZr-round_mean_Zr)^2 + 31 * (round_mean_UpperGallupZr-round_mean_Zr)^2 + 25 * (round_mean_ShaleZr-round_mean_Zr)^2 + 28* (round_mean_FBase3Zr-round_mean_Zr)^2
+ 29 * (round_mean_FBase2Zr-round_mean_Zr)^2 + 26 * (round_mean_FBase1Zr-round_mean_Zr)^2 + 26 * (round_mean_EBaseZr-round_mean_Zr)^2
+ 25 * (round_mean_B6Zr-round_mean_Zr)^2 + 30 * (round_mean_TorrivioZr-round_mean_Zr)^2)/11

YVar <- (27 * (round_mean_SSTY-round_mean_Y)^2 + 25 * (round_mean_ShaleY-round_mean_Y)^2 + 31 * (round_mean_UpperGallupY-round_mean_Y)^2 + 28* (round_mean_FBase3Y-round_mean_Y)^2
+ 29 * (round_mean_FBase2Y-round_mean_Y)^2 + 26 * (round_mean_FBase1Y-round_mean_Y)^2 + 26 * (round_mean_EBaseY-round_mean_Y)^2
+ 25 * (round_mean_B6Y-round_mean_Y)^2 + 30 * (round_mean_TorrivioY-round_mean_Y)^2)/11

CeVar <- (27 * (round_mean_SSTCe-round_mean_Ce)^2 + 25 * (round_mean_ShaleCe-round_mean_Ce)^2 + 31 * (round_mean_UpperGallupCe-round_mean_Ce)^2 + 28* (round_mean_FBase3Ce-round_mean_Ce)^2
+ 29 * (round_mean_FBase2Ce-round_mean_Ce)^2 + 26 * (round_mean_FBase1Ce-round_mean_Ce)^2 + 26 * (round_mean_EBaseCe-round_mean_Ce)^2
+ 25 * (round_mean_B6Ce-round_mean_Ce)^2 + 30 * (round_mean_TorrivioCe-round_mean_Ce)^2)/11

SVar <- (27 * (round_mean_SSTS-round_mean_S)^2 + 31 * (round_mean_UpperGallupS-round_mean_S)^2 + 25 * (round_mean_ShaleS-round_mean_S)^2 + 28* (round_mean_FBase3S-round_mean_S)^2
+ 29 * (round_mean_FBase2S-round_mean_S)^2 + 26 * (round_mean_FBase1S-round_mean_S)^2 + 26 * (round_mean_EBaseS-round_mean_S)^2
+ 25 * (round_mean_B6S-round_mean_S)^2 + 30 * (round_mean_TorrivioS-round_mean_S)^2)/11

CaVar <- (27 * (round_mean_SSTCa-round_mean_Ca)^2 + 31 * (round_mean_UpperGallupCa-round_mean_Ca)^2 + 25 * (round_mean_ShaleCa-round_mean_Ca)^2 + 28* (round_mean_FBase3Ca-round_mean_Ca)^2
+ 29 * (round_mean_FBase2Ca-round_mean_Ca)^2 + 26 * (round_mean_FBase1Ca-round_mean_Ca)^2 + 26 * (round_mean_EBaseCa-round_mean_Ca)^2
+ 25 * (round_mean_B6Ca-round_mean_Ca)^2 + 30 * (round_mean_TorrivioCa-round_mean_Ca)^2)/11

LaVar <- (27 * (round_mean_SSTLa-round_mean_La)^2 + 31 * (round_mean_UpperGallupLa-round_mean_La)^2 + 25 * (round_mean_ShaleLa-round_mean_La)^2 + 28* (round_mean_FBase3La-round_mean_La)^2
          + 29 * (round_mean_FBase2La-round_mean_La)^2 + 26 * (round_mean_FBase1La-round_mean_La)^2 + 26 * (round_mean_EBaseLa-round_mean_La)^2
          + 25 * (round_mean_B6La-round_mean_La)^2 + 30 * (round_mean_TorrivioLa-round_mean_La)^2)/11

#Calculate the variance for each element between the each of the sample groupings
FewithinVarSST <- sum((SSTSubset$Fe-round_mean_SSTFe)^2) 
FewithinVarShale <- sum((ShaleSubset$Fe-round_mean_ShaleFe)^2)
FewithinVarEBase <- sum((EBaseSubset$Fe-round_mean_EBaseFe)^2) 
FewithinVarFBase1 <- sum((FBase1Subset$Fe-round_mean_FBase1Fe)^2)
FewithinVarFBase2 <- sum((FBase2Subset$Fe-round_mean_FBase2Fe)^2)
FewithinVarFBase3 <- sum((FBase3Subset$Fe-round_mean_FBase3Fe)^2)
FewithinVarTorrivio <- sum((TorrivioSubset$Fe-round_mean_TorrivioFe)^2)
FewithinVarB6 <- sum((B6Subset$Fe-round_mean_B6Fe)^2)
FewithinVarUpperGallup <- sum((UpperGallupSubset$Fe-round_mean_UpperGallupFe)^2)

CewithinVarSST <- sum((SSTSubset$Ce-round_mean_SSTCe)^2) 
CewithinVarShale <- sum((ShaleSubset$Ce-round_mean_ShaleCe)^2)
CewithinVarEBase <- sum((EBaseSubset$Ce-round_mean_EBaseCe)^2) 
CewithinVarFBase1 <- sum((FBase1Subset$Ce-round_mean_FBase1Ce)^2)
CewithinVarFBase2 <- sum((FBase2Subset$Ce-round_mean_FBase2Ce)^2)
CewithinVarFBase3 <- sum((FBase3Subset$Ce-round_mean_FBase3Ce)^2)
CewithinVarTorrivio <- sum((TorrivioSubset$Ce-round_mean_TorrivioCe)^2)
CewithinVarB6 <- sum((B6Subset$Ce-round_mean_B6Ce)^2)
CewithinVarUpperGallup <- sum((UpperGallupSubset$Ce-round_mean_UpperGallupCe)^2)

NiwithinVarSST <- sum((SSTSubset$Ni-round_mean_SSTNi)^2) 
NiwithinVarShale <- sum((ShaleSubset$Ni-round_mean_ShaleNi)^2)
NiwithinVarEBase <- sum((EBaseSubset$Ni-round_mean_EBaseNi)^2) 
NiwithinVarFBase1 <- sum((FBase1Subset$Ni-round_mean_FBase1Ni)^2)
NiwithinVarFBase2 <- sum((FBase2Subset$Ni-round_mean_FBase2Ni)^2)
NiwithinVarFBase3 <- sum((FBase3Subset$Ni-round_mean_FBase3Ni)^2)
NiwithinVarTorrivio <- sum((TorrivioSubset$Ni-round_mean_TorrivioNi)^2)
NiwithinVarB6 <- sum((B6Subset$Ni-round_mean_B6Ni)^2)
NiwithinVarUpperGallup <- sum((UpperGallupSubset$Ni-round_mean_UpperGallupNi)^2)

SwithinVarSST <- sum((SSTSubset$S-round_mean_SSTS)^2) 
SwithinVarShale <- sum((ShaleSubset$S-round_mean_ShaleS)^2)
SwithinVarEBase <- sum((EBaseSubset$S-round_mean_EBaseS)^2) 
SwithinVarFBase1 <- sum((FBase1Subset$S-round_mean_FBase1S)^2)
SwithinVarFBase2 <- sum((FBase2Subset$S-round_mean_FBase2S)^2)
SwithinVarFBase3 <- sum((FBase3Subset$S-round_mean_FBase3S)^2)
SwithinVarTorrivio <- sum((TorrivioSubset$S-round_mean_TorrivioS)^2)
SwithinVarB6 <- sum((B6Subset$S-round_mean_B6S)^2)
SwithinVarUpperGallup <- sum((UpperGallupSubset$S-round_mean_UpperGallupS)^2)

ZrwithinVarSST <- sum((SSTSubset$Zr-round_mean_SSTZr)^2) 
ZrwithinVarShale <- sum((ShaleSubset$Zr-round_mean_ShaleZr)^2)
ZrwithinVarEBase <- sum((EBaseSubset$Zr-round_mean_EBaseZr)^2) 
ZrwithinVarFBase1 <- sum((FBase1Subset$Zr-round_mean_FBase1Zr)^2)
ZrwithinVarFBase2 <- sum((FBase2Subset$Zr-round_mean_FBase2Zr)^2)
ZrwithinVarFBase3 <- sum((FBase3Subset$Zr-round_mean_FBase3Zr)^2)
ZrwithinVarTorrivio <- sum((TorrivioSubset$Zr-round_mean_TorrivioZr)^2)
ZrwithinVarB6 <- sum((B6Subset$Zr-round_mean_B6Zr)^2)
ZrwithinVarUpperGallup <- sum((UpperGallupSubset$Zr-round_mean_UpperGallupZr)^2)

CawithinVarSST <- sum((SSTSubset$Ca-round_mean_SSTCa)^2) 
CawithinVarShale <- sum((ShaleSubset$Ca-round_mean_ShaleCa)^2)
CawithinVarEBase <- sum((EBaseSubset$Ca-round_mean_EBaseCa)^2) 
CawithinVarFBase1 <- sum((FBase1Subset$Ca-round_mean_FBase1Ca)^2)
CawithinVarFBase2 <- sum((FBase2Subset$Ca-round_mean_FBase2Ca)^2)
CawithinVarFBase3 <- sum((FBase3Subset$Ca-round_mean_FBase3Ca)^2)
CawithinVarTorrivio <- sum((TorrivioSubset$Ca-round_mean_TorrivioCa)^2)
CawithinVarB6 <- sum((B6Subset$Ca-round_mean_B6Ca)^2)
CawithinVarUpperGallup <- sum((UpperGallupSubset$Ca-round_mean_UpperGallupCa)^2)

YwithinVarSST <- sum((SSTSubset$Y-round_mean_SSTY)^2) 
YwithinVarShale <- sum((ShaleSubset$Y-round_mean_ShaleY)^2)
YwithinVarEBase <- sum((EBaseSubset$Y-round_mean_EBaseY)^2) 
YwithinVarFBase1 <- sum((FBase1Subset$Y-round_mean_FBase1Y)^2)
YwithinVarFBase2 <- sum((FBase2Subset$Y-round_mean_FBase2Y)^2)
YwithinVarFBase3 <- sum((FBase3Subset$Y-round_mean_FBase3Y)^2)
YwithinVarTorrivio <- sum((TorrivioSubset$Y-round_mean_TorrivioY)^2)
YwithinVarB6 <- sum((B6Subset$Y-round_mean_B6Y)^2)
YwithinVarUpperGallup <- sum((UpperGallupSubset$Y-round_mean_UpperGallupY)^2)

MnwithinVarSST <- sum((SSTSubset$Mn-round_mean_SSTMn)^2) 
MnwithinVarShale <- sum((ShaleSubset$Mn-round_mean_ShaleMn)^2)
MnwithinVarEBase <- sum((EBaseSubset$Mn-round_mean_EBaseMn)^2) 
MnwithinVarFBase1 <- sum((FBase1Subset$Mn-round_mean_FBase1Mn)^2)
MnwithinVarFBase2 <- sum((FBase2Subset$Mn-round_mean_FBase2Mn)^2)
MnwithinVarFBase3 <- sum((FBase3Subset$Mn-round_mean_FBase3Mn)^2)
MnwithinVarTorrivio <- sum((TorrivioSubset$Mn-round_mean_TorrivioMn)^2)
MnwithinVarB6 <- sum((B6Subset$Mn-round_mean_B6Mn)^2)
MnwithinVarUpperGallup <- sum((UpperGallupSubset$Mn-round_mean_UpperGallupMn)^2)

NdwithinVarSST <- sum((SSTSubset$Nd-round_mean_SSTNd)^2) 
NdwithinVarShale <- sum((ShaleSubset$Nd-round_mean_ShaleNd)^2)
NdwithinVarEBase <- sum((EBaseSubset$Nd-round_mean_EBaseNd)^2) 
NdwithinVarFBase1 <- sum((FBase1Subset$Nd-round_mean_FBase1Nd)^2)
NdwithinVarFBase2 <- sum((FBase2Subset$Nd-round_mean_FBase2Nd)^2)
NdwithinVarFBase3 <- sum((FBase3Subset$Nd-round_mean_FBase3Nd)^2)
NdwithinVarTorrivio <- sum((TorrivioSubset$Nd-round_mean_TorrivioNd)^2)
NdwithinVarB6 <- sum((B6Subset$Nd-round_mean_B6Nd)^2)
NdwithinVarUpperGallup <- sum((UpperGallupSubset$Nd-round_mean_UpperGallupNd)^2)

ThwithinVarSST <- sum((SSTSubset$Th-round_mean_SSTTh)^2) 
ThwithinVarShale <- sum((ShaleSubset$Th-round_mean_ShaleTh)^2)
ThwithinVarEBase <- sum((EBaseSubset$Th-round_mean_EBaseTh)^2) 
ThwithinVarFBase1 <- sum((FBase1Subset$Th-round_mean_FBase1Th)^2)
ThwithinVarFBase2 <- sum((FBase2Subset$Th-round_mean_FBase2Th)^2)
ThwithinVarFBase3 <- sum((FBase3Subset$Th-round_mean_FBase3Th)^2)
ThwithinVarTorrivio <- sum((TorrivioSubset$Th-round_mean_TorrivioTh)^2)
ThwithinVarB6 <- sum((B6Subset$Th-round_mean_B6Th)^2)
ThwithinVarUpperGallup <- sum((UpperGallupSubset$Th-round_mean_UpperGallupTh)^2)

LawithinVarSST <- sum((SSTSubset$La-round_mean_SSTLa)^2) 
LawithinVarShale <- sum((ShaleSubset$La-round_mean_ShaleLa)^2)
LawithinVarEBase <- sum((EBaseSubset$La-round_mean_EBaseLa)^2) 
LawithinVarFBase1 <- sum((FBase1Subset$La-round_mean_FBase1La)^2)
LawithinVarFBase2 <- sum((FBase2Subset$La-round_mean_FBase2La)^2)
LawithinVarFBase3 <- sum((FBase3Subset$La-round_mean_FBase3La)^2)
LawithinVarTorrivio <- sum((TorrivioSubset$La-round_mean_TorrivioLa)^2)
LawithinVarB6 <- sum((B6Subset$La-round_mean_B6La)^2)
LawithinVarUpperGallup <- sum((UpperGallupSubset$La-round_mean_UpperGallupLa)^2)

#Calculate the variance for each element for each sample grouping
within_group_varFe <- ((FewithinVarB6+FewithinVarUpperGallup+FewithinVarFBase1+FewithinVarEBase+FewithinVarFBase2+FewithinVarFBase3+FewithinVarShale
                        +FewithinVarSST+FewithinVarTorrivio)/((31+27+25+28+29+26+26+25+30)-9))

within_group_varTh <- ((ThwithinVarB6+ThwithinVarUpperGallup+ThwithinVarFBase1+ThwithinVarEBase+ThwithinVarFBase2+ThwithinVarFBase3+ThwithinVarShale
                        +ThwithinVarSST+ThwithinVarTorrivio)/((31+27+25+28+29+26+26+25+30)-9))

within_group_varNi <- ((NiwithinVarB6+NiwithinVarUpperGallup+NiwithinVarFBase1+NiwithinVarEBase+NiwithinVarFBase2+NiwithinVarFBase3+NiwithinVarShale
                        +NiwithinVarSST+NiwithinVarTorrivio)/((31+27+25+28+29+26+26+25+30)-9))

within_group_varNd <- ((NdwithinVarB6+NdwithinVarUpperGallup+NdwithinVarFBase1+NdwithinVarEBase+NdwithinVarFBase2+NdwithinVarFBase3+NdwithinVarShale
                        +NdwithinVarSST+NdwithinVarTorrivio)/((31+27+25+28+29+26+26+25+30)-9))

within_group_varS <- ((SwithinVarB6+SwithinVarUpperGallup+SwithinVarFBase1+SwithinVarEBase+SwithinVarFBase2+SwithinVarFBase3+SwithinVarShale
                       +SwithinVarSST+SwithinVarTorrivio)/((31+27+25+28+29+26+26+25+30)-9))

within_group_varMn <- ((MnwithinVarB6+MnwithinVarUpperGallup+MnwithinVarFBase1+MnwithinVarEBase+MnwithinVarFBase2+MnwithinVarFBase3+MnwithinVarShale
                        +MnwithinVarSST+MnwithinVarTorrivio)/((31+27+25+28+29+26+26+25+30)-9))

within_group_varCe <- ((CewithinVarB6+CewithinVarUpperGallup+CewithinVarFBase1+CewithinVarEBase+CewithinVarFBase2+CewithinVarFBase3+CewithinVarShale
                        +CewithinVarSST+CewithinVarTorrivio)/((31+27+25+28+29+26+26+25+30)-9))

within_group_varZr <- ((ZrwithinVarB6+ZrwithinVarUpperGallup+ZrwithinVarFBase1+ZrwithinVarEBase+ZrwithinVarFBase2+ZrwithinVarFBase3+ZrwithinVarShale
                        +ZrwithinVarSST+ZrwithinVarTorrivio)/((31+27+25+28+29+26+26+25+30)-9))

within_group_varCa <- ((CawithinVarB6+CawithinVarUpperGallup+CawithinVarFBase1+CawithinVarEBase+CawithinVarFBase2+CawithinVarFBase3+CawithinVarShale
                        +CawithinVarSST+CawithinVarTorrivio)/((31+27+25+28+29+26+26+25+30)-9))

within_group_varY <- ((YwithinVarB6+YwithinVarUpperGallup+YwithinVarFBase1+YwithinVarEBase+YwithinVarFBase2+YwithinVarFBase3+YwithinVarShale
                       +YwithinVarSST+YwithinVarTorrivio)/((31+27+25+28+29+26+26+25+30)-9))

within_group_varLa <- ((LawithinVarB6+LawithinVarUpperGallup+LawithinVarFBase1+LawithinVarEBase+LawithinVarFBase2+LawithinVarFBase3+LawithinVarShale
                        +LawithinVarSST+LawithinVarTorrivio)/((31+27+25+28+29+26+26+25+30)-9))

#Calculate F-Statistic
F_Fe <- FeVar/within_group_varFe

F_Ni <- NiVar/within_group_varNi

F_Th <- ThVar/within_group_varTh

F_Nd <- NdVar/within_group_varNd

F_Y <- YVar/within_group_varY

F_Ce <- CeVar/within_group_varCe

F_Ca <- CaVar/within_group_varCa

F_Zr <- ZrVar/within_group_varZr

F_S <- SVar/within_group_varS

F_Mn <- MnVar/within_group_varMn

F_La <- LaVar/within_group_varLa

#Calculate first degree of freedom
df1 <- 9-1

#Calculate second degree of freedom
df2 <- (31+27+25+28+29+26+26+25+30)-9
  
#Calculate p-values
  pFe <- pf(F_Fe,df1,df2, lower.tail = FALSE)
  pMn <- pf(F_Mn,df1,df2, lower.tail = FALSE)
  pS <- pf(F_S,df1,df2, lower.tail = FALSE)
  pNi <- pf(F_Ni,df1,df2, lower.tail = FALSE)
  pCa <- pf(F_Ca,df1,df2, lower.tail = FALSE)
  pCe <- pf(F_Ce,df1,df2, lower.tail = FALSE)
  pTh <- pf(F_Th,df1,df2, lower.tail = FALSE)
  pZr <- pf(F_Zr,df1,df2, lower.tail = FALSE)
  pY <- pf(F_Y,df1,df2, lower.tail = FALSE)
  pNd <- pf(F_Nd,df1,df2, lower.tail = FALSE)
  pLa <- pf(F_La,df1,df2, lower.tail = FALSE)
  
  
  library(readxl)
  Bentonite2 <- read_excel("C:/Bentonite Scan Moly Tube Test 2.xlsx")
  
  BentonitePCA2 <- Bentonite2[,1:57]
  install.packages("tidyverse")
  
  Bentonite2.pr <- prcomp(BentonitePCA2, center = TRUE, scale = TRUE)
  summary(Bentonite2.pr)
  
  screeplot(Bentonite2.pr, type = "l", npcs = 58, main = "Screeplot of the 58 PCs")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(Bentonite.pr$sdev^2 / sum(Bentonite.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
  abline(v = 6, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC6"),
         col=c("blue"), lty=5, cex=0.6)
  
  #Calculate the rounded mean for Fe, Th, Ni, Mn, Nd, Zr, Y, S, Ca, Ce, and La
  round_mean_Fe2 <- round(mean(BentonitePCA2$Fe), digits =2)
  round_mean_Th2 <- round(mean(BentonitePCA2$Th), digits =2)
  round_mean_Ni2 <- round(mean(BentonitePCA2$Ni), digits =2)
  round_mean_Mn2 <- round(mean(BentonitePCA2$Mn), digits =2)
  round_mean_Nd2 <- round(mean(BentonitePCA2$Nd), digits =2)
  round_mean_Zr2 <- round(mean(BentonitePCA2$Zr), digits =2)
  round_mean_Y2 <- round(mean(BentonitePCA2$Y), digits =2)
  round_mean_S2 <- round(mean(BentonitePCA2$S), digits =2)
  round_mean_Ca2 <- round(mean(BentonitePCA2$Ca), digits =2)
  round_mean_Cr2 <- round(mean(BentonitePCA2$Cr), digits =2)
  round_mean_Ce2 <- round(mean(BentonitePCA2$Ce), digits =2)
  round_mean_La2 <- round(mean(BentonitePCA2$La), digits =2)
  
  #Subset Juana Lopez Data according to bentonite sample groups
  JL_4_A <- BentonitePCA2[1:26,1:57]
  
  JL_4_B <- BentonitePCA2[27:55,1:57]
  
  JL_1 <- BentonitePCA2[56:84,1:57]
  
  JL_2 <- BentonitePCA2[85:109,1:57]
  
  JL_3 <- BentonitePCA2[110:137,1:57]
  
  JL_4_C <- BentonitePCA2[138:165,1:57]
  
  JL_4_D <- BentonitePCA2[166:194,1:57]
  
  
  FeBentonite <- read_excel("C:/Fe Table.xlsx")
  
  #Calculate correlation matrices 
  pairs.panels(FeBentonite[1:7])
  
  CaBentonite <- read_excel("C:/Ca Table.xlsx")
  
  #Calculate correlation matrices 
  pairs.panels(CaBentonite[1:7])
  
  KBentonite <- read_excel("C:/K Table.xlsx")
  
  #Calculate correlation matrices 
  pairs.panels(KBentonite[1:7])
  
  ZrBentonite <- read_excel("C:/Zr Table.xlsx")
  
  #Calculate correlation matrices 
  pairs.panels(ZrBentonite[1:7])
  
  CrBentonite <- read_excel("C:/Cr Table.xlsx")
  
  #Calculate correlation matrices 
  pairs.panels(CrBentonite[1:7])
  
  SBentonite <- read_excel("C:/S Table.xlsx")
  
  #Calculate correlation matrices 
  pairs.panels(SBentonite[1:7])
  
  MnBentonite <- read_excel("C:/Mn Table.xlsx")
  
  #Calculate correlation matrices 
  pairs.panels(MnBentonite[1:7])
  
  TiBentonite <- read_excel("C:/Ti Table.xlsx")
  
  #Calculate correlation matrices 
  pairs.panels(TiBentonite[1:7])
  
  
  #Calculate correlation matrices for all elements for each sample group
  library(psych)
  pairs.panels(JL_1[2:10])
  pairs.panels(JL_1[11:21])
  pairs.panels(JL_1[22:32])
  pairs.panels(JL_1[33:43])
  pairs.panels(JL_1[44:54])
  pairs.panels(JL_1[55:57])
  
  pairs.panels(JL_2[2:10])
  pairs.panels(JL_2[11:21])
  pairs.panels(JL_2[22:32])
  pairs.panels(JL_2[33:43])
  pairs.panels(JL_2[44:54])
  pairs.panels(JL_2[55:57])
  
  pairs.panels(JL_3[2:10])
  pairs.panels(JL_3[11:21])
  pairs.panels(JL_3[22:32])
  pairs.panels(JL_3[33:43])
  pairs.panels(JL_3[44:54])
  pairs.panels(JL_3[55:57])
  
  pairs.panels(JL_4_A[2:10])
  pairs.panels(JL_4_A[11:21])
  pairs.panels(JL_4_A[22:32])
  pairs.panels(JL_4_A[33:43])
  pairs.panels(JL_4_A[44:54])
  pairs.panels(JL_4_A[55:57])
  
  pairs.panels(JL_4_B[2:10])
  pairs.panels(JL_4_B[11:21])
  pairs.panels(JL_4_B[22:32])
  pairs.panels(JL_4_B[33:43])
  pairs.panels(JL_4_B[44:54])
  pairs.panels(JL_4_B[55:57])
  
  pairs.panels(JL_4_C[2:10])
  pairs.panels(JL_4_C[11:21])
  pairs.panels(JL_4_C[22:32])
  pairs.panels(JL_4_C[33:43])
  pairs.panels(JL_4_C[44:54])
  pairs.panels(JL_4_C[55:57])
  
  pairs.panels(JL_4_D[2:10])
  pairs.panels(JL_4_D[11:21])
  pairs.panels(JL_4_D[22:32])
  pairs.panels(JL_4_D[33:43])
  pairs.panels(JL_4_D[44:54])
  pairs.panels(JL_4_D[55:57])
  
  JL_1Subset <- JL_1[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Cr", "Ce","La")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_1Subset[1:11])
  
  JL_1Terra <- JL_1[,c("Al","K","Si","Ca","Ti","Cl","S","Zr")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_1Terra[1:8])
  
  JL_2Subset <- JL_2[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Cr", "Ce","La")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_2Subset[1:11])
  
  JL_2Terra <- JL_2[,c("Al","K","Si","Ca","Ti","Cl","S","Zr")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_2Terra[1:8])
  
  JL_3Subset <- JL_3[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Cr", "Ce","La")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_3Subset[1:11])
  
  JL_3Terra <- JL_3[,c("Al","K","Si","Ca","Ti","Cl","S","Zr")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_3Terra[1:8])
  
  JL_4_ASubset <- JL_4_A[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Cr", "Ce","La")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_4_ASubset[1:11])
  
  JL_4_ATerra <- JL_4_A[,c("Al","K","Si","Ca","Ti","Cl","S","Zr")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_4_ATerra[1:8])
  
  JL_4_BSubset <- JL_4_B[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Cr", "Ce","La")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_4_BSubset[1:11])
  
  JL_4_BTerra <- JL_4_B[,c("Al","K","Si","Ca","Ti","Cl","S","Zr")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_4_BTerra[1:8])
  
  JL_4_CSubset <- JL_4_C[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Cr", "Ce","La")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_4_CSubset[1:11])
  
  JL_4_CTerra <- JL_4_C[,c("Al","K","Si","Ca","Ti","Cl","S","Zr")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_4_CTerra[1:8])
  
  JL_4_DSubset <- JL_4_D[,c("Fe","Th","Ni", "Mn", "Nd", "Zr", "Y", "S", "Ca", "Cr", "Ce","La")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_4_DSubset[1:11])
  
  
  JL_4_DTerra <- JL_4_D[,c("Al","K","Si","Ca","Ti","Cl","S","Zr")]
  
  #Calculate correlation matrices 
  pairs.panels(JL_4_DTerra[1:8])
  
  #Calculate the rounded mean for Fe, Mn, S, Zr, Ce, Ca, Nd, Hf, Th,La, Ni, and Y for JL-1, JL-2, JL-3, JL-4A, JL-4B, JL-4C, and JL-4D
  round_mean_JL_1Fe <- round(mean(JL_1Subset$Fe), digits =2)
  round_mean_JL_1Th <- round(mean(JL_1Subset$Th), digits =2)
  round_mean_JL_1Ni <- round(mean(JL_1Subset$Ni), digits =2)
  round_mean_JL_1Mn <- round(mean(JL_1Subset$Mn), digits =2)
  round_mean_JL_1Nd <- round(mean(JL_1Subset$Nd), digits =2)
  round_mean_JL_1Zr <- round(mean(JL_1Subset$Zr), digits =2)
  round_mean_JL_1Y <- round(mean(JL_1Subset$Y), digits =2)
  round_mean_JL_1S <- round(mean(JL_1Subset$S), digits =2)
  round_mean_JL_1Ca <- round(mean(JL_1Subset$Ca), digits =2)
  round_mean_JL_1Cr <- round(mean(JL_1Subset$Cr), digits =2)
  round_mean_JL_1Ce <- round(mean(JL_1Subset$Ce), digits =2)
  round_mean_JL_1La <- round(mean(JL_1Subset$La), digits =2)
  round_mean_JL_1U <- round(mean(JL_1Subset$U), digits =2)
  round_mean_JL_1Hf <- round(mean(JL_1Subset$Hf), digits =2)
  
  round_mean_JL_2Fe <- round(mean(JL_2Subset$Fe), digits =2)
  round_mean_JL_2Th <- round(mean(JL_2Subset$Th), digits =2)
  round_mean_JL_2Ni <- round(mean(JL_2Subset$Ni), digits =2)
  round_mean_JL_2Mn <- round(mean(JL_2Subset$Mn), digits =2)
  round_mean_JL_2Nd <- round(mean(JL_2Subset$Nd), digits =2)
  round_mean_JL_2Zr <- round(mean(JL_2Subset$Zr), digits =2)
  round_mean_JL_2Y <- round(mean(JL_2Subset$Y), digits =2)
  round_mean_JL_2S <- round(mean(JL_2Subset$S), digits =2)
  round_mean_JL_2Ca <- round(mean(JL_2Subset$Ca), digits =2)
  round_mean_JL_2Cr <- round(mean(JL_2Subset$Cr), digits =2)
  round_mean_JL_2Ce <- round(mean(JL_2Subset$Ce), digits =2)
  round_mean_JL_2La <- round(mean(JL_2Subset$La), digits =2)
  
  round_mean_JL_3Fe <- round(mean(JL_3Subset$Fe), digits =2)
  round_mean_JL_3Th <- round(mean(JL_3Subset$Th), digits =2)
  round_mean_JL_3Ni <- round(mean(JL_3Subset$Ni), digits =2)
  round_mean_JL_3Mn <- round(mean(JL_3Subset$Mn), digits =2)
  round_mean_JL_3Nd <- round(mean(JL_3Subset$Nd), digits =2)
  round_mean_JL_3Zr <- round(mean(JL_3Subset$Zr), digits =2)
  round_mean_JL_3Y <- round(mean(JL_3Subset$Y), digits =2)
  round_mean_JL_3S <- round(mean(JL_3Subset$S), digits =2)
  round_mean_JL_3Ca <- round(mean(JL_3Subset$Ca), digits =2)
  round_mean_JL_3Cr <- round(mean(JL_3Subset$Cr), digits =2)
  round_mean_JL_3Ce <- round(mean(JL_3Subset$Ce), digits =2)
  round_mean_JL_3La <- round(mean(JL_3Subset$La), digits =2)
  
  round_mean_JL_4_AFe <- round(mean(JL_4_ASubset$Fe), digits =2)
  round_mean_JL_4_ATh <- round(mean(JL_4_ASubset$Th), digits =2)
  round_mean_JL_4_ANi <- round(mean(JL_4_ASubset$Ni), digits =2)
  round_mean_JL_4_AMn <- round(mean(JL_4_ASubset$Mn), digits =2)
  round_mean_JL_4_ANd <- round(mean(JL_4_ASubset$Nd), digits =2)
  round_mean_JL_4_AZr <- round(mean(JL_4_ASubset$Zr), digits =2)
  round_mean_JL_4_AY <- round(mean(JL_4_ASubset$Y), digits =2)
  round_mean_JL_4_AS <- round(mean(JL_4_ASubset$S), digits =2)
  round_mean_JL_4_ACa <- round(mean(JL_4_ASubset$Ca), digits =2)
  round_mean_JL_4_ACr <- round(mean(JL_4_ASubset$Cr), digits =2)
  round_mean_JL_4_ACe <- round(mean(JL_4_ASubset$Ce), digits =2)
  round_mean_JL_4_ALa <- round(mean(JL_4_ASubset$La), digits =2)
  
  round_mean_JL_4_BFe <- round(mean(JL_4_BSubset$Fe), digits =2)
  round_mean_JL_4_BTh <- round(mean(JL_4_BSubset$Th), digits =2)
  round_mean_JL_4_BNi <- round(mean(JL_4_BSubset$Ni), digits =2)
  round_mean_JL_4_BMn <- round(mean(JL_4_BSubset$Mn), digits =2)
  round_mean_JL_4_BNd <- round(mean(JL_4_BSubset$Nd), digits =2)
  round_mean_JL_4_BZr <- round(mean(JL_4_BSubset$Zr), digits =2)
  round_mean_JL_4_BY <- round(mean(JL_4_BSubset$Y), digits =2)
  round_mean_JL_4_BS <- round(mean(JL_4_BSubset$S), digits =2)
  round_mean_JL_4_BCa <- round(mean(JL_4_BSubset$Ca), digits =2)
  round_mean_JL_4_BCr <- round(mean(JL_4_BSubset$Cr), digits =2)
  round_mean_JL_4_BCe <- round(mean(JL_4_BSubset$Ce), digits =2)
  round_mean_JL_4_BLa <- round(mean(JL_4_BSubset$La), digits =2)
  
  round_mean_JL_4_CFe <- round(mean(JL_4_CSubset$Fe), digits =2)
  round_mean_JL_4_CTh <- round(mean(JL_4_CSubset$Th), digits =2)
  round_mean_JL_4_CNi <- round(mean(JL_4_CSubset$Ni), digits =2)
  round_mean_JL_4_CMn <- round(mean(JL_4_CSubset$Mn), digits =2)
  round_mean_JL_4_CNd <- round(mean(JL_4_CSubset$Nd), digits =2)
  round_mean_JL_4_CZr <- round(mean(JL_4_CSubset$Zr), digits =2)
  round_mean_JL_4_CY <- round(mean(JL_4_CSubset$Y), digits =2)
  round_mean_JL_4_CS <- round(mean(JL_4_CSubset$S), digits =2)
  round_mean_JL_4_CCa <- round(mean(JL_4_CSubset$Ca), digits =2)
  round_mean_JL_4_CCr <- round(mean(JL_4_CSubset$Cr), digits =2)
  round_mean_JL_4_CCe <- round(mean(JL_4_CSubset$Ce), digits =2)
  round_mean_JL_4_CLa <- round(mean(JL_4_CSubset$La), digits =2)
  
  round_mean_JL_4_DFe <- round(mean(JL_4_DSubset$Fe), digits =2)
  round_mean_JL_4_DTh <- round(mean(JL_4_DSubset$Th), digits =2)
  round_mean_JL_4_DNi <- round(mean(JL_4_DSubset$Ni), digits =2)
  round_mean_JL_4_DMn <- round(mean(JL_4_DSubset$Mn), digits =2)
  round_mean_JL_4_DNd <- round(mean(JL_4_DSubset$Nd), digits =2)
  round_mean_JL_4_DZr <- round(mean(JL_4_DSubset$Zr), digits =2)
  round_mean_JL_4_DY <- round(mean(JL_4_DSubset$Y), digits =2)
  round_mean_JL_4_DS <- round(mean(JL_4_DSubset$S), digits =2)
  round_mean_JL_4_DCa <- round(mean(JL_4_DSubset$Ca), digits =2)
  round_mean_JL_4_DCr <- round(mean(JL_4_DSubset$Cr), digits =2)
  round_mean_JL_4_DCe <- round(mean(JL_4_DSubset$Ce), digits =2)
  round_mean_JL_4_DLa <- round(mean(JL_4_DSubset$La), digits =2)
  
  #Calculate the variance for Fe, Mn, S, Zr, Ce, Ca, Nd, Hf, Th,La, Ni, and Y
  MnVar <- (29 * (round_mean_JL_1Mn-round_mean_Mn2)^2 + 25 * (round_mean_JL_2Mn-round_mean_Mn2)^2 + 28* (round_mean_JL_3Mn-round_mean_Mn2)^2
            + 26 * (round_mean_JL_4_AMn-round_mean_Mn2)^2 + 29 * (round_mean_JL_4_BMn-round_mean_Mn2)^2 + 28 * (round_mean_JL_4_CMn-round_mean_Mn2)^2
            + 29 * (round_mean_JL_4_DMn-round_mean_Mn2)^2/11)
  
  AlVar <- (29 * (round_mean_JL_1Al-round_mean_Al2)^2 + 25 * (round_mean_JL_2Al-round_mean_Al2)^2 + 28* (round_mean_JL_3Al-round_mean_Al2)^2
            + 26 * (round_mean_JL_4_AAl-round_mean_Al2)^2 + 29 * (round_mean_JL_4_BAl-round_mean_Al2)^2 + 28 * (round_mean_JL_4_CAl-round_mean_Al2)^2
            + 29 * (round_mean_JL_4_DAl-round_mean_Al2)^2/11)
  
  SVar <- (29 * (round_mean_JL_1S-round_mean_S2)^2 + 25 * (round_mean_JL_2S-round_mean_S2)^2 + 28* (round_mean_JL_3S-round_mean_S2)^2
           + 26 * (round_mean_JL_4_AS-round_mean_S2)^2 + 29 * (round_mean_JL_4_BS-round_mean_S2)^2 + 28 * (round_mean_JL_4_CS-round_mean_S2)^2
           + 29 * (round_mean_JL_4_DS-round_mean_S2)^2/11)
  
  TiVar <- (29 * (round_mean_JL_1Ti-round_mean_Ti2)^2 + 25 * (round_mean_JL_2Ti-round_mean_Ti2)^2 + 28* (round_mean_JL_3Ti-round_mean_Ti2)^2
            + 26 * (round_mean_JL_4_ATi-round_mean_Ti2)^2 + 29 * (round_mean_JL_4_BTi-round_mean_Ti2)^2 + 28 * (round_mean_JL_4_CTi-round_mean_Ti2)^2
            + 29 * (round_mean_JL_4_DTi-round_mean_Ti2)^2/11)
  
  CeVar <- (29 * (round_mean_JL_1Ce-round_mean_Ce2)^2 + 25 * (round_mean_JL_2Ce-round_mean_Ce2)^2 + 28* (round_mean_JL_3Ce-round_mean_Ce2)^2
            + 26 * (round_mean_JL_4_ACe-round_mean_Ce2)^2 + 29 * (round_mean_JL_4_BCe-round_mean_Ce2)^2 + 28 * (round_mean_JL_4_CCe-round_mean_Ce2)^2
            + 29 * (round_mean_JL_4_DCe-round_mean_Ce2)^2/11)
  
  NiVar <- (29 * (round_mean_JL_1Ni-round_mean_Ni2)^2 + 25 * (round_mean_JL_2Ni-round_mean_Ni2)^2 + 28* (round_mean_JL_3Ni-round_mean_Ni2)^2
            + 26 * (round_mean_JL_4_ANi-round_mean_Ni2)^2 + 29 * (round_mean_JL_4_BNi-round_mean_Ni2)^2 + 28 * (round_mean_JL_4_CNi-round_mean_Ni2)^2
            + 29 * (round_mean_JL_4_DNi-round_mean_Ni2)^2/11)
  
  NdVar <- (29 * (round_mean_JL_1Nd-round_mean_Nd2)^2 + 25 * (round_mean_JL_2Nd-round_mean_Nd2)^2 + 28* (round_mean_JL_3Nd-round_mean_Nd2)^2
            + 26 * (round_mean_JL_4_ANd-round_mean_Nd2)^2 + 29 * (round_mean_JL_4_BNd-round_mean_Nd2)^2 + 28 * (round_mean_JL_4_CNd-round_mean_Nd2)^2
            + 29 * (round_mean_JL_4_DNd-round_mean_Nd2)^2/11)
  
  HfVar <- (29 * (round_mean_JL_1Hf-round_mean_Hf2)^2 + 25 * (round_mean_JL_2Hf-round_mean_Hf2)^2 + 28* (round_mean_JL_3Hf-round_mean_Hf2)^2
            + 26 * (round_mean_JL_4_AHf-round_mean_Hf2)^2 + 29 * (round_mean_JL_4_BHf-round_mean_Hf2)^2 + 28 * (round_mean_JL_4_CHf-round_mean_Hf2)^2
            + 29 * (round_mean_JL_4_DHf-round_mean_Hf2)^2/11)
  
  UVar <- (29 * (round_mean_JL_1U-round_mean_U2)^2 + 25 * (round_mean_JL_2U-round_mean_U2)^2 + 28* (round_mean_JL_3U-round_mean_U2)^2
           + 26 * (round_mean_JL_4_AU-round_mean_U2)^2 + 29 * (round_mean_JL_4_BU-round_mean_U2)^2 + 28 * (round_mean_JL_4_CU-round_mean_U2)^2
           + 29 * (round_mean_JL_4_DU-round_mean_U2)^2/11)
  
  ThVar <- (29 * (round_mean_JL_1Th-round_mean_Th2)^2 + 25 * (round_mean_JL_2Th-round_mean_Th2)^2 + 28* (round_mean_JL_3Th-round_mean_Th2)^2
            + 26 * (round_mean_JL_4_ATh-round_mean_Th2)^2 + 29 * (round_mean_JL_4_BTh-round_mean_Th2)^2 + 28 * (round_mean_JL_4_CTh-round_mean_Th2)^2
            + 29 * (round_mean_JL_4_DTh-round_mean_Th2)^2/11)
  
  CrVar <- (29 * (round_mean_JL_1Cr-round_mean_Cr2)^2 + 25 * (round_mean_JL_2Cr-round_mean_Cr2)^2 + 28* (round_mean_JL_3Cr-round_mean_Cr2)^2
            + 26 * (round_mean_JL_4_ACr-round_mean_Cr2)^2 + 29 * (round_mean_JL_4_BCr-round_mean_Cr2)^2 + 28 * (round_mean_JL_4_CCr-round_mean_Cr2)^2
            + 29 * (round_mean_JL_4_DCr-round_mean_Cr2)^2/11)
  
  CaVar <- (29 * (round_mean_JL_1Ca-round_mean_Ca2)^2 + 25 * (round_mean_JL_2Ca-round_mean_Ca2)^2 + 28* (round_mean_JL_3Ca-round_mean_Ca2)^2
            + 26 * (round_mean_JL_4_ACa-round_mean_Ca2)^2 + 29 * (round_mean_JL_4_BCa-round_mean_Ca2)^2 + 28 * (round_mean_JL_4_CCa-round_mean_Ca2)^2
            + 29 * (round_mean_JL_4_DCa-round_mean_Ca2)^2/11)
  
  KVar <- (29 * (round_mean_JL_1K-round_mean_K2)^2 + 25 * (round_mean_JL_2K-round_mean_K2)^2 + 28* (round_mean_JL_3K-round_mean_K2)^2
           + 26 * (round_mean_JL_4_AK-round_mean_K2)^2 + 29 * (round_mean_JL_4_BK-round_mean_K2)^2 + 28 * (round_mean_JL_4_CK-round_mean_K2)^2
           + 29 * (round_mean_JL_4_DK-round_mean_K2)^2/11)
  
  ZrVar <- (29 * (round_mean_JL_1Zr-round_mean_Zr2)^2 + 25 * (round_mean_JL_2Zr-round_mean_Zr2)^2 + 28* (round_mean_JL_3Zr-round_mean_Zr2)^2
            + 26 * (round_mean_JL_4_AZr-round_mean_Zr2)^2 + 29 * (round_mean_JL_4_BZr-round_mean_Zr2)^2 + 28 * (round_mean_JL_4_CZr-round_mean_Zr2)^2
            + 29 * (round_mean_JL_4_DZr-round_mean_Zr2)^2/11)
  
  FeVar <- (29 * (round_mean_JL_1Fe-round_mean_Fe2)^2 + 25 * (round_mean_JL_2Fe-round_mean_Fe2)^2 + 28* (round_mean_JL_3Fe-round_mean_Fe2)^2
            + 26 * (round_mean_JL_4_AFe-round_mean_Fe2)^2 + 29 * (round_mean_JL_4_BFe-round_mean_Fe2)^2 + 28 * (round_mean_JL_4_CFe-round_mean_Fe2)^2
            + 29 * (round_mean_JL_4_DFe-round_mean_Fe2)^2/11)
  
  LaVar <- (29 * (round_mean_JL_1La-round_mean_La2)^2 + 25 * (round_mean_JL_2La-round_mean_La2)^2 + 28* (round_mean_JL_3La-round_mean_La2)^2
            + 26 * (round_mean_JL_4_ALa-round_mean_La2)^2 + 29 * (round_mean_JL_4_BLa-round_mean_La2)^2 + 28 * (round_mean_JL_4_CLa-round_mean_La2)^2
            + 29 * (round_mean_JL_4_DLa-round_mean_La2)^2/11)
  
  YVar <- (29 * (round_mean_JL_1Y-round_mean_Y2)^2 + 25 * (round_mean_JL_2Y-round_mean_Y2)^2 + 28* (round_mean_JL_3Y-round_mean_Y2)^2
           + 26 * (round_mean_JL_4_AY-round_mean_Y2)^2 + 29 * (round_mean_JL_4_BY-round_mean_Y2)^2 + 28 * (round_mean_JL_4_CY-round_mean_Y2)^2
           + 29 * (round_mean_JL_4_DY-round_mean_Y2)^2/11)
  
  #Calculate the within group variance for Fe, Mn, S, Zr, Ce, Ca, Nd, Hf, Th,La, Ni, and Y
  FewithinVarJL_1 <- sum((JL_1Subset$Fe-round_mean_JL_1Fe)^2) 
  FewithinVarJL_2 <- sum((JL_2Subset$Fe-round_mean_JL_2Fe)^2)
  FewithinVarJL_3 <- sum((JL_3Subset$Fe-round_mean_JL_3Fe)^2) 
  FewithinVarJL_4_A <- sum((JL_4_ASubset$Fe-round_mean_JL_4_AFe)^2)
  FewithinVarJL_4_B <- sum((JL_4_BSubset$Fe-round_mean_JL_4_BFe)^2)
  FewithinVarJL_4_C <- sum((JL_4_CSubset$Fe-round_mean_JL_4_CFe)^2)
  FewithinVarJL_4_D <- sum((JL_4_DSubset$Fe-round_mean_JL_4_DFe)^2)
  
  MnwithinVarJL_1 <- sum((JL_1Subset$Mn-round_mean_JL_1Mn)^2) 
  MnwithinVarJL_2 <- sum((JL_2Subset$Mn-round_mean_JL_2Mn)^2)
  MnwithinVarJL_3 <- sum((JL_3Subset$Mn-round_mean_JL_3Mn)^2) 
  MnwithinVarJL_4_A <- sum((JL_4_ASubset$Mn-round_mean_JL_4_AMn)^2)
  MnwithinVarJL_4_B <- sum((JL_4_BSubset$Mn-round_mean_JL_4_BMn)^2)
  MnwithinVarJL_4_C <- sum((JL_4_CSubset$Mn-round_mean_JL_4_CMn)^2)
  MnwithinVarJL_4_D <- sum((JL_4_DSubset$Mn-round_mean_JL_4_DMn)^2)
  
  SwithinVarJL_1 <- sum((JL_1Subset$S-round_mean_JL_1S)^2) 
  SwithinVarJL_2 <- sum((JL_2Subset$S-round_mean_JL_2S)^2)
  SwithinVarJL_3 <- sum((JL_3Subset$S-round_mean_JL_3S)^2) 
  SwithinVarJL_4_A <- sum((JL_4_ASubset$S-round_mean_JL_4_AS)^2)
  SwithinVarJL_4_B <- sum((JL_4_BSubset$S-round_mean_JL_4_BS)^2)
  SwithinVarJL_4_C <- sum((JL_4_CSubset$S-round_mean_JL_4_CS)^2)
  SwithinVarJL_4_D <- sum((JL_4_DSubset$S-round_mean_JL_4_DS)^2)
  
  ZrwithinVarJL_1 <- sum((JL_1Subset$Zr-round_mean_JL_1Zr)^2) 
  ZrwithinVarJL_2 <- sum((JL_2Subset$Zr-round_mean_JL_2Zr)^2)
  ZrwithinVarJL_3 <- sum((JL_3Subset$Zr-round_mean_JL_3Zr)^2) 
  ZrwithinVarJL_4_A <- sum((JL_4_ASubset$Zr-round_mean_JL_4_AZr)^2)
  ZrwithinVarJL_4_B <- sum((JL_4_BSubset$Zr-round_mean_JL_4_BZr)^2)
  ZrwithinVarJL_4_C <- sum((JL_4_CSubset$Zr-round_mean_JL_4_CZr)^2)
  ZrwithinVarJL_4_D <- sum((JL_4_DSubset$Zr-round_mean_JL_4_DZr)^2)
  
  CewithinVarJL_1 <- sum((JL_1Subset$Ce-round_mean_JL_1Ce)^2) 
  CewithinVarJL_2 <- sum((JL_2Subset$Ce-round_mean_JL_2Ce)^2)
  CewithinVarJL_3 <- sum((JL_3Subset$Ce-round_mean_JL_3Ce)^2) 
  CewithinVarJL_4_A <- sum((JL_4_ASubset$Ce-round_mean_JL_4_ACe)^2)
  CewithinVarJL_4_B <- sum((JL_4_BSubset$Ce-round_mean_JL_4_BCe)^2)
  CewithinVarJL_4_C <- sum((JL_4_CSubset$Ce-round_mean_JL_4_CCe)^2)
  CewithinVarJL_4_D <- sum((JL_4_DSubset$Ce-round_mean_JL_4_DCe)^2)
  
  CawithinVarJL_1 <- sum((JL_1Subset$Ca-round_mean_JL_1Ca)^2) 
  CawithinVarJL_2 <- sum((JL_2Subset$Ca-round_mean_JL_2Ca)^2)
  CawithinVarJL_3 <- sum((JL_3Subset$Ca-round_mean_JL_3Ca)^2) 
  CawithinVarJL_4_A <- sum((JL_4_ASubset$Ca-round_mean_JL_4_ACa)^2)
  CawithinVarJL_4_B <- sum((JL_4_BSubset$Ca-round_mean_JL_4_BCa)^2)
  CawithinVarJL_4_C <- sum((JL_4_CSubset$Ca-round_mean_JL_4_CCa)^2)
  CawithinVarJL_4_D <- sum((JL_4_DSubset$Ca-round_mean_JL_4_DCa)^2)
  
  NdwithinVarJL_1 <- sum((JL_1Subset$Nd-round_mean_JL_1Nd)^2) 
  NdwithinVarJL_2 <- sum((JL_2Subset$Nd-round_mean_JL_2Nd)^2)
  NdwithinVarJL_3 <- sum((JL_3Subset$Nd-round_mean_JL_3Nd)^2) 
  NdwithinVarJL_4_A <- sum((JL_4_ASubset$Nd-round_mean_JL_4_ANd)^2)
  NdwithinVarJL_4_B <- sum((JL_4_BSubset$Nd-round_mean_JL_4_BNd)^2)
  NdwithinVarJL_4_C <- sum((JL_4_CSubset$Nd-round_mean_JL_4_CNd)^2)
  NdwithinVarJL_4_D <- sum((JL_4_DSubset$Nd-round_mean_JL_4_DNd)^2)
  
  HfwithinVarJL_1 <- sum((JL_1Subset$Hf-round_mean_JL_1Hf)^2) 
  HfwithinVarJL_2 <- sum((JL_2Subset$Hf-round_mean_JL_2Hf)^2)
  HfwithinVarJL_3 <- sum((JL_3Subset$Hf-round_mean_JL_3Hf)^2) 
  HfwithinVarJL_4_A <- sum((JL_4_ASubset$Hf-round_mean_JL_4_AHf)^2)
  HfwithinVarJL_4_B <- sum((JL_4_BSubset$Hf-round_mean_JL_4_BHf)^2)
  HfwithinVarJL_4_C <- sum((JL_4_CSubset$Hf-round_mean_JL_4_CHf)^2)
  HfwithinVarJL_4_D <- sum((JL_4_DSubset$Hf-round_mean_JL_4_DHf)^2)
  
  ThwithinVarJL_1 <- sum((JL_1Subset$Th-round_mean_JL_1Th)^2) 
  ThwithinVarJL_2 <- sum((JL_2Subset$Th-round_mean_JL_2Th)^2)
  ThwithinVarJL_3 <- sum((JL_3Subset$Th-round_mean_JL_3Th)^2) 
  ThwithinVarJL_4_A <- sum((JL_4_ASubset$Th-round_mean_JL_4_ATh)^2)
  ThwithinVarJL_4_B <- sum((JL_4_BSubset$Th-round_mean_JL_4_BTh)^2)
  ThwithinVarJL_4_C <- sum((JL_4_CSubset$Th-round_mean_JL_4_CTh)^2)
  ThwithinVarJL_4_D <- sum((JL_4_DSubset$Th-round_mean_JL_4_DTh)^2)
  
  LawithinVarJL_1 <- sum((JL_1Subset$La-round_mean_JL_1La)^2) 
  LawithinVarJL_2 <- sum((JL_2Subset$La-round_mean_JL_2La)^2)
  LawithinVarJL_3 <- sum((JL_3Subset$La-round_mean_JL_3La)^2) 
  LawithinVarJL_4_A <- sum((JL_4_ASubset$La-round_mean_JL_4_ALa)^2)
  LawithinVarJL_4_B <- sum((JL_4_BSubset$La-round_mean_JL_4_BLa)^2)
  LawithinVarJL_4_C <- sum((JL_4_CSubset$La-round_mean_JL_4_CLa)^2)
  LawithinVarJL_4_D <- sum((JL_4_DSubset$La-round_mean_JL_4_DLa)^2)
  
  NiwithinVarJL_1 <- sum((JL_1Subset$Ni-round_mean_JL_1Ni)^2) 
  NiwithinVarJL_2 <- sum((JL_2Subset$Ni-round_mean_JL_2Ni)^2)
  NiwithinVarJL_3 <- sum((JL_3Subset$Ni-round_mean_JL_3Ni)^2) 
  NiwithinVarJL_4_A <- sum((JL_4_ASubset$Ni-round_mean_JL_4_ANi)^2)
  NiwithinVarJL_4_B <- sum((JL_4_BSubset$Ni-round_mean_JL_4_BNi)^2)
  NiwithinVarJL_4_C <- sum((JL_4_CSubset$Ni-round_mean_JL_4_CNi)^2)
  NiwithinVarJL_4_D <- sum((JL_4_DSubset$Ni-round_mean_JL_4_DNi)^2)
  
  YwithinVarJL_1 <- sum((JL_1Subset$Y-round_mean_JL_1Y)^2) 
  YwithinVarJL_2 <- sum((JL_2Subset$Y-round_mean_JL_2Y)^2)
  YwithinVarJL_3 <- sum((JL_3Subset$Y-round_mean_JL_3Y)^2) 
  YwithinVarJL_4_A <- sum((JL_4_ASubset$Y-round_mean_JL_4_AY)^2)
  YwithinVarJL_4_B <- sum((JL_4_BSubset$Y-round_mean_JL_4_BY)^2)
  YwithinVarJL_4_C <- sum((JL_4_CSubset$Y-round_mean_JL_4_CY)^2)
  YwithinVarJL_4_D <- sum((JL_4_DSubset$Y-round_mean_JL_4_DY)^2)
  
  #Calculate the between group variance for Fe, Mn, S, Zr, Ce, Ca, Nd, Hf, Th,La, Ni, and Y
  within_group_varFe <- ((FewithinVarJL_1+FewithinVarJL_2+FewithinVarJL_3+FewithinVarJL_4_A+FewithinVarJL_4_B
                          +FewithinVarJL_4_C+FewithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varY <- ((YwithinVarJL_1+YwithinVarJL_2+YwithinVarJL_3+YwithinVarJL_4_A+YwithinVarJL_4_B
                         +YwithinVarJL_4_C+YwithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varTh <- ((ThwithinVarJL_1+ThwithinVarJL_2+ThwithinVarJL_3+ThwithinVarJL_4_A+ThwithinVarJL_4_B
                          +ThwithinVarJL_4_C+ThwithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varU <- ((UwithinVarJL_1+UwithinVarJL_2+UwithinVarJL_3+UwithinVarJL_4_A+UwithinVarJL_4_B
                         +UwithinVarJL_4_C+UwithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varCe <- ((CewithinVarJL_1+CewithinVarJL_2+CewithinVarJL_3+CewithinVarJL_4_A+CewithinVarJL_4_B
                          +CewithinVarJL_4_C+CewithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varCa <- ((CawithinVarJL_1+CawithinVarJL_2+CawithinVarJL_3+CawithinVarJL_4_A+CawithinVarJL_4_B
                          +CawithinVarJL_4_C+CawithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varCr <- ((CrwithinVarJL_1+CrwithinVarJL_2+CrwithinVarJL_3+CrwithinVarJL_4_A+CrwithinVarJL_4_B
                          +CrwithinVarJL_4_C+CrwithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varZr <- ((ZrwithinVarJL_1+ZrwithinVarJL_2+ZrwithinVarJL_3+ZrwithinVarJL_4_A+ZrwithinVarJL_4_B
                          +ZrwithinVarJL_4_C+ZrwithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varLa <- ((LawithinVarJL_1+LawithinVarJL_2+LawithinVarJL_3+LawithinVarJL_4_A+LawithinVarJL_4_B
                          +LawithinVarJL_4_C+LawithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varNd <- ((NdwithinVarJL_1+NdwithinVarJL_2+NdwithinVarJL_3+NdwithinVarJL_4_A+NdwithinVarJL_4_B
                          +NdwithinVarJL_4_C+NdwithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varS <- ((SwithinVarJL_1+SwithinVarJL_2+SwithinVarJL_3+SwithinVarJL_4_A+SwithinVarJL_4_B
                         +SwithinVarJL_4_C+SwithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varLa <- ((LawithinVarJL_1+LawithinVarJL_2+LawithinVarJL_3+LawithinVarJL_4_A+LawithinVarJL_4_B
                          +LawithinVarJL_4_C+LawithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varMn <- ((MnwithinVarJL_1+MnwithinVarJL_2+MnwithinVarJL_3+MnwithinVarJL_4_A+MnwithinVarJL_4_B
                          +MnwithinVarJL_4_C+MnwithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  within_group_varNi <- ((NiwithinVarJL_1+NiwithinVarJL_2+NiwithinVarJL_3+NiwithinVarJL_4_A+NiwithinVarJL_4_B
                          +NiwithinVarJL_4_C+NiwithinVarJL_4_D)/((29+25+28+26+29+28+29)-7))
  
  #Calculate the F Statistic for Fe, Mn, S, Zr, Ce, Ca, Nd, Hf, Th,La, Ni, and Y
  F_Fe <- FeVar/within_group_varFe
  
  F_Ni <- NiVar/within_group_varNi
  
  F_Th <- ThVar/within_group_varTh
  
  F_Nd <- NdVar/within_group_varNd
  
  F_Y <- YVar/within_group_varY
  
  F_Cr <- CrVar/within_group_varCr
  
  F_Ce <- CeVar/within_group_varCe
  
  F_Ca <- CaVar/within_group_varCa
  
  F_Zr <- ZrVar/within_group_varZr
  
  F_S <- SVar/within_group_varS
  
  F_Mn <- MnVar/within_group_varMn
  
  F_La <- LaVar/within_group_varLa
  
  #Calculate the first degree of freedom
  df1 <- 7-1
  
  #Calculate the second degree of freedom
  df2 <- (29+25+28+26+29+28+29)-7
  
  #Calculate the p-value
  pFe <- pf(F_Fe,df1,df2, lower.tail = FALSE)
  pMn <- pf(F_Mn,df1,df2, lower.tail = FALSE)
  pS <- pf(F_S,df1,df2, lower.tail = FALSE)
  pNi <- pf(F_Ni,df1,df2, lower.tail = FALSE)
  pCa <- pf(F_Ca,df1,df2, lower.tail = FALSE)
  pCe <- pf(F_Ce,df1,df2, lower.tail = FALSE)
  pTh <- pf(F_Th,df1,df2, lower.tail = FALSE)
  pZr <- pf(F_Zr,df1,df2, lower.tail = FALSE)
  pY <- pf(F_Y,df1,df2, lower.tail = FALSE)
  pNd <- pf(F_Nd,df1,df2, lower.tail = FALSE)
  pLa <- pf(F_La,df1,df2, lower.tail = FALSE)
  
  
  ##Juana Lopez F Statistic Calculations
  
  #Calculate the variance for Fe, Mn, S, Zr, Ce, Ca, Nd, Hf, Th,La, Ni, and Y
  MnVar <- (26 * (round_mean_JL_4_AMn-round_mean_Mn2)^2 + 29 * (round_mean_JL_4_BMn-round_mean_Mn2)^2 + 28 * (round_mean_JL_4_CMn-round_mean_Mn2)^2
            + 29 * (round_mean_JL_4_DMn-round_mean_Mn2)^2/11)
  
  AlVar <- (26 * (round_mean_JL_4_AAl-round_mean_Al2)^2 + 29 * (round_mean_JL_4_BAl-round_mean_Al2)^2 + 28 * (round_mean_JL_4_CAl-round_mean_Al2)^2
            + 29 * (round_mean_JL_4_DAl-round_mean_Al2)^2/11)
  
  SVar <- (26 * (round_mean_JL_4_AS-round_mean_S2)^2 + 29 * (round_mean_JL_4_BS-round_mean_S2)^2 + 28 * (round_mean_JL_4_CS-round_mean_S2)^2
           + 29 * (round_mean_JL_4_DS-round_mean_S2)^2/11)
  
  TiVar <- (26 * (round_mean_JL_4_ATi-round_mean_Ti2)^2 + 29 * (round_mean_JL_4_BTi-round_mean_Ti2)^2 + 28 * (round_mean_JL_4_CTi-round_mean_Ti2)^2
            + 29 * (round_mean_JL_4_DTi-round_mean_Ti2)^2/11)
  
  CeVar <- (26 * (round_mean_JL_4_ACe-round_mean_Ce2)^2 + 29 * (round_mean_JL_4_BCe-round_mean_Ce2)^2 + 28 * (round_mean_JL_4_CCe-round_mean_Ce2)^2
            + 29 * (round_mean_JL_4_DCe-round_mean_Ce2)^2/11)
  
  NiVar <- (26 * (round_mean_JL_4_ANi-round_mean_Ni2)^2 + 29 * (round_mean_JL_4_BNi-round_mean_Ni2)^2 + 28 * (round_mean_JL_4_CNi-round_mean_Ni2)^2
            + 29 * (round_mean_JL_4_DNi-round_mean_Ni2)^2/11)
  
  NdVar <- (26 * (round_mean_JL_4_ANd-round_mean_Nd2)^2 + 29 * (round_mean_JL_4_BNd-round_mean_Nd2)^2 + 28 * (round_mean_JL_4_CNd-round_mean_Nd2)^2
            + 29 * (round_mean_JL_4_DNd-round_mean_Nd2)^2/11)
  
  HfVar <- (26 * (round_mean_JL_4_AHf-round_mean_Hf2)^2 + 29 * (round_mean_JL_4_BHf-round_mean_Hf2)^2 + 28 * (round_mean_JL_4_CHf-round_mean_Hf2)^2
            + 29 * (round_mean_JL_4_DHf-round_mean_Hf2)^2/11)
  
  UVar <- (26 * (round_mean_JL_4_AU-round_mean_U2)^2 + 29 * (round_mean_JL_4_BU-round_mean_U2)^2 + 28 * (round_mean_JL_4_CU-round_mean_U2)^2
           + 29 * (round_mean_JL_4_DU-round_mean_U2)^2/11)
  
  ThVar <- (26 * (round_mean_JL_4_ATh-round_mean_Th2)^2 + 29 * (round_mean_JL_4_BTh-round_mean_Th2)^2 + 28 * (round_mean_JL_4_CTh-round_mean_Th2)^2
            + 29 * (round_mean_JL_4_DTh-round_mean_Th2)^2/11)

  CaVar <- (26 * (round_mean_JL_4_ACa-round_mean_Ca2)^2 + 29 * (round_mean_JL_4_BCa-round_mean_Ca2)^2 + 28 * (round_mean_JL_4_CCa-round_mean_Ca2)^2
            + 29 * (round_mean_JL_4_DCa-round_mean_Ca2)^2/11)
  
  KVar <- (26 * (round_mean_JL_4_AK-round_mean_K2)^2 + 29 * (round_mean_JL_4_BK-round_mean_K2)^2 + 28 * (round_mean_JL_4_CK-round_mean_K2)^2
           + 29 * (round_mean_JL_4_DK-round_mean_K2)^2/11)
  
  ZrVar <- (26 * (round_mean_JL_4_AZr-round_mean_Zr2)^2 + 29 * (round_mean_JL_4_BZr-round_mean_Zr2)^2 + 28 * (round_mean_JL_4_CZr-round_mean_Zr2)^2
            + 29 * (round_mean_JL_4_DZr-round_mean_Zr2)^2/11)
  
  FeVar <- (26 * (round_mean_JL_4_AFe-round_mean_Fe2)^2 + 29 * (round_mean_JL_4_BFe-round_mean_Fe2)^2 + 28 * (round_mean_JL_4_CFe-round_mean_Fe2)^2
            + 29 * (round_mean_JL_4_DFe-round_mean_Fe2)^2/11)
  
  LaVar <- (26 * (round_mean_JL_4_ALa-round_mean_La2)^2 + 29 * (round_mean_JL_4_BLa-round_mean_La2)^2 + 28 * (round_mean_JL_4_CLa-round_mean_La2)^2
            + 29 * (round_mean_JL_4_DLa-round_mean_La2)^2/11)
  
  YVar <- (26 * (round_mean_JL_4_AY-round_mean_Y2)^2 + 29 * (round_mean_JL_4_BY-round_mean_Y2)^2 + 28 * (round_mean_JL_4_CY-round_mean_Y2)^2
           + 29 * (round_mean_JL_4_DY-round_mean_Y2)^2/11)
  
  
  #Calculate the within group variance for Fe, Mn, S, Zr, Ce, Ca, Nd, Hf, Th,La, Ni, and Y
  FewithinVarJL_1 <- sum((JL_1Subset$Fe-round_mean_JL_1Fe)^2) 
  FewithinVarJL_2 <- sum((JL_2Subset$Fe-round_mean_JL_2Fe)^2)
  FewithinVarJL_3 <- sum((JL_3Subset$Fe-round_mean_JL_3Fe)^2) 
  FewithinVarJL_4_A <- sum((JL_4_ASubset$Fe-round_mean_JL_4_AFe)^2)
  FewithinVarJL_4_B <- sum((JL_4_BSubset$Fe-round_mean_JL_4_BFe)^2)
  FewithinVarJL_4_C <- sum((JL_4_CSubset$Fe-round_mean_JL_4_CFe)^2)
  FewithinVarJL_4_D <- sum((JL_4_DSubset$Fe-round_mean_JL_4_DFe)^2)
  
  MnwithinVarJL_1 <- sum((JL_1Subset$Mn-round_mean_JL_1Mn)^2) 
  MnwithinVarJL_2 <- sum((JL_2Subset$Mn-round_mean_JL_2Mn)^2)
  MnwithinVarJL_3 <- sum((JL_3Subset$Mn-round_mean_JL_3Mn)^2) 
  MnwithinVarJL_4_A <- sum((JL_4_ASubset$Mn-round_mean_JL_4_AMn)^2)
  MnwithinVarJL_4_B <- sum((JL_4_BSubset$Mn-round_mean_JL_4_BMn)^2)
  MnwithinVarJL_4_C <- sum((JL_4_CSubset$Mn-round_mean_JL_4_CMn)^2)
  MnwithinVarJL_4_D <- sum((JL_4_DSubset$Mn-round_mean_JL_4_DMn)^2)
  
  SwithinVarJL_1 <- sum((JL_1Subset$S-round_mean_JL_1S)^2) 
  SwithinVarJL_2 <- sum((JL_2Subset$S-round_mean_JL_2S)^2)
  SwithinVarJL_3 <- sum((JL_3Subset$S-round_mean_JL_3S)^2) 
  SwithinVarJL_4_A <- sum((JL_4_ASubset$S-round_mean_JL_4_AS)^2)
  SwithinVarJL_4_B <- sum((JL_4_BSubset$S-round_mean_JL_4_BS)^2)
  SwithinVarJL_4_C <- sum((JL_4_CSubset$S-round_mean_JL_4_CS)^2)
  SwithinVarJL_4_D <- sum((JL_4_DSubset$S-round_mean_JL_4_DS)^2)
  
  ZrwithinVarJL_1 <- sum((JL_1Subset$Zr-round_mean_JL_1Zr)^2) 
  ZrwithinVarJL_2 <- sum((JL_2Subset$Zr-round_mean_JL_2Zr)^2)
  ZrwithinVarJL_3 <- sum((JL_3Subset$Zr-round_mean_JL_3Zr)^2) 
  ZrwithinVarJL_4_A <- sum((JL_4_ASubset$Zr-round_mean_JL_4_AZr)^2)
  ZrwithinVarJL_4_B <- sum((JL_4_BSubset$Zr-round_mean_JL_4_BZr)^2)
  ZrwithinVarJL_4_C <- sum((JL_4_CSubset$Zr-round_mean_JL_4_CZr)^2)
  ZrwithinVarJL_4_D <- sum((JL_4_DSubset$Zr-round_mean_JL_4_DZr)^2)
  
  CewithinVarJL_1 <- sum((JL_1Subset$Ce-round_mean_JL_1Ce)^2) 
  CewithinVarJL_2 <- sum((JL_2Subset$Ce-round_mean_JL_2Ce)^2)
  CewithinVarJL_3 <- sum((JL_3Subset$Ce-round_mean_JL_3Ce)^2) 
  CewithinVarJL_4_A <- sum((JL_4_ASubset$Ce-round_mean_JL_4_ACe)^2)
  CewithinVarJL_4_B <- sum((JL_4_BSubset$Ce-round_mean_JL_4_BCe)^2)
  CewithinVarJL_4_C <- sum((JL_4_CSubset$Ce-round_mean_JL_4_CCe)^2)
  CewithinVarJL_4_D <- sum((JL_4_DSubset$Ce-round_mean_JL_4_DCe)^2)
  
  CawithinVarJL_1 <- sum((JL_1Subset$Ca-round_mean_JL_1Ca)^2) 
  CawithinVarJL_2 <- sum((JL_2Subset$Ca-round_mean_JL_2Ca)^2)
  CawithinVarJL_3 <- sum((JL_3Subset$Ca-round_mean_JL_3Ca)^2) 
  CawithinVarJL_4_A <- sum((JL_4_ASubset$Ca-round_mean_JL_4_ACa)^2)
  CawithinVarJL_4_B <- sum((JL_4_BSubset$Ca-round_mean_JL_4_BCa)^2)
  CawithinVarJL_4_C <- sum((JL_4_CSubset$Ca-round_mean_JL_4_CCa)^2)
  CawithinVarJL_4_D <- sum((JL_4_DSubset$Ca-round_mean_JL_4_DCa)^2)
  
  NdwithinVarJL_1 <- sum((JL_1Subset$Nd-round_mean_JL_1Nd)^2) 
  NdwithinVarJL_2 <- sum((JL_2Subset$Nd-round_mean_JL_2Nd)^2)
  NdwithinVarJL_3 <- sum((JL_3Subset$Nd-round_mean_JL_3Nd)^2) 
  NdwithinVarJL_4_A <- sum((JL_4_ASubset$Nd-round_mean_JL_4_ANd)^2)
  NdwithinVarJL_4_B <- sum((JL_4_BSubset$Nd-round_mean_JL_4_BNd)^2)
  NdwithinVarJL_4_C <- sum((JL_4_CSubset$Nd-round_mean_JL_4_CNd)^2)
  NdwithinVarJL_4_D <- sum((JL_4_DSubset$Nd-round_mean_JL_4_DNd)^2)
  
  HfwithinVarJL_1 <- sum((JL_1Subset$Hf-round_mean_JL_1Hf)^2) 
  HfwithinVarJL_2 <- sum((JL_2Subset$Hf-round_mean_JL_2Hf)^2)
  HfwithinVarJL_3 <- sum((JL_3Subset$Hf-round_mean_JL_3Hf)^2) 
  HfwithinVarJL_4_A <- sum((JL_4_ASubset$Hf-round_mean_JL_4_AHf)^2)
  HfwithinVarJL_4_B <- sum((JL_4_BSubset$Hf-round_mean_JL_4_BHf)^2)
  HfwithinVarJL_4_C <- sum((JL_4_CSubset$Hf-round_mean_JL_4_CHf)^2)
  HfwithinVarJL_4_D <- sum((JL_4_DSubset$Hf-round_mean_JL_4_DHf)^2)
  
  ThwithinVarJL_1 <- sum((JL_1Subset$Th-round_mean_JL_1Th)^2) 
  ThwithinVarJL_2 <- sum((JL_2Subset$Th-round_mean_JL_2Th)^2)
  ThwithinVarJL_3 <- sum((JL_3Subset$Th-round_mean_JL_3Th)^2) 
  ThwithinVarJL_4_A <- sum((JL_4_ASubset$Th-round_mean_JL_4_ATh)^2)
  ThwithinVarJL_4_B <- sum((JL_4_BSubset$Th-round_mean_JL_4_BTh)^2)
  ThwithinVarJL_4_C <- sum((JL_4_CSubset$Th-round_mean_JL_4_CTh)^2)
  ThwithinVarJL_4_D <- sum((JL_4_DSubset$Th-round_mean_JL_4_DTh)^2)
  
  UwithinVarJL_1 <- sum((JL_1Subset$U-round_mean_JL_1U)^2) 
  UwithinVarJL_2 <- sum((JL_2Subset$U-round_mean_JL_2U)^2)
  UwithinVarJL_3 <- sum((JL_3Subset$U-round_mean_JL_3U)^2) 
  UwithinVarJL_4_A <- sum((JL_4_ASubset$U-round_mean_JL_4_AU)^2)
  UwithinVarJL_4_B <- sum((JL_4_BSubset$U-round_mean_JL_4_BU)^2)
  UwithinVarJL_4_C <- sum((JL_4_CSubset$U-round_mean_JL_4_CU)^2)
  UwithinVarJL_4_D <- sum((JL_4_DSubset$U-round_mean_JL_4_DU)^2)
  
  LawithinVarJL_1 <- sum((JL_1Subset$La-round_mean_JL_1La)^2) 
  LawithinVarJL_2 <- sum((JL_2Subset$La-round_mean_JL_2La)^2)
  LawithinVarJL_3 <- sum((JL_3Subset$La-round_mean_JL_3La)^2) 
  LawithinVarJL_4_A <- sum((JL_4_ASubset$La-round_mean_JL_4_ALa)^2)
  LawithinVarJL_4_B <- sum((JL_4_BSubset$La-round_mean_JL_4_BLa)^2)
  LawithinVarJL_4_C <- sum((JL_4_CSubset$La-round_mean_JL_4_CLa)^2)
  LawithinVarJL_4_D <- sum((JL_4_DSubset$La-round_mean_JL_4_DLa)^2)
  
  NiwithinVarJL_1 <- sum((JL_1Subset$Ni-round_mean_JL_1Ni)^2) 
  NiwithinVarJL_2 <- sum((JL_2Subset$Ni-round_mean_JL_2Ni)^2)
  NiwithinVarJL_3 <- sum((JL_3Subset$Ni-round_mean_JL_3Ni)^2) 
  NiwithinVarJL_4_A <- sum((JL_4_ASubset$Ni-round_mean_JL_4_ANi)^2)
  NiwithinVarJL_4_B <- sum((JL_4_BSubset$Ni-round_mean_JL_4_BNi)^2)
  NiwithinVarJL_4_C <- sum((JL_4_CSubset$Ni-round_mean_JL_4_CNi)^2)
  NiwithinVarJL_4_D <- sum((JL_4_DSubset$Ni-round_mean_JL_4_DNi)^2)
  
  YwithinVarJL_1 <- sum((JL_1Subset$Y-round_mean_JL_1Y)^2) 
  YwithinVarJL_2 <- sum((JL_2Subset$Y-round_mean_JL_2Y)^2)
  YwithinVarJL_3 <- sum((JL_3Subset$Y-round_mean_JL_3Y)^2) 
  YwithinVarJL_4_A <- sum((JL_4_ASubset$Y-round_mean_JL_4_AY)^2)
  YwithinVarJL_4_B <- sum((JL_4_BSubset$Y-round_mean_JL_4_BY)^2)
  YwithinVarJL_4_C <- sum((JL_4_CSubset$Y-round_mean_JL_4_CY)^2)
  YwithinVarJL_4_D <- sum((JL_4_DSubset$Y-round_mean_JL_4_DY)^2)
  
  #Calculate the between group variance for Fe, Mn, S, Zr, Ce, Ca, Nd, Hf, Th,La, Ni, and Y
  within_group_varFe <- ((FewithinVarJL_4_A+FewithinVarJL_4_B
                          +FewithinVarJL_4_C+FewithinVarJL_4_D)/((26+29+28+29)-4))
  
  within_group_varY <- ((YwithinVarJL_4_A+YwithinVarJL_4_B
                         +YwithinVarJL_4_C+YwithinVarJL_4_D)/((26+29+28+29)-4))
  
  within_group_varTh <- ((ThwithinVarJL_4_A+ThwithinVarJL_4_B
                          +ThwithinVarJL_4_C+ThwithinVarJL_4_D)/((26+29+28+29)-4))
  
  within_group_varU <- ((UwithinVarJL_4_A+UwithinVarJL_4_B
                         +UwithinVarJL_4_C+UwithinVarJL_4_D)/((26+29+28+29)-4))
  
  within_group_varCe <- ((CewithinVarJL_4_A+CewithinVarJL_4_B
                          +CewithinVarJL_4_C+CewithinVarJL_4_D)/((26+29+28+29)-4))
  
  within_group_varCa <- ((CawithinVarJL_4_A+CawithinVarJL_4_B
                          +CawithinVarJL_4_C+CawithinVarJL_4_D)/((26+29+28+29)-4))
  
  within_group_varZr <- ((ZrwithinVarJL_4_A+ZrwithinVarJL_4_B
                          +ZrwithinVarJL_4_C+ZrwithinVarJL_4_D)/((26+29+28+29)-4))
  
  within_group_varLa <- ((LawithinVarJL_4_A+LawithinVarJL_4_B
                          +LawithinVarJL_4_C+LawithinVarJL_4_D)/((26+29+28+29)-4))
  
  within_group_varNd <- ((NdwithinVarJL_4_A+NdwithinVarJL_4_B
                          +NdwithinVarJL_4_C+NdwithinVarJL_4_D)/((26+29+28+29)-4))
  
  within_group_varS <- ((SwithinVarJL_4_A+SwithinVarJL_4_B
                         +SwithinVarJL_4_C+SwithinVarJL_4_D)/((26+29+28+29)-4))
  
  within_group_varLa <- ((LawithinVarJL_4_A+LawithinVarJL_4_B
                          +LawithinVarJL_4_C+LawithinVarJL_4_D)/((26+29+28+29)-4))
  
  within_group_varMn <- ((MnwithinVarJL_4_A+MnwithinVarJL_4_B
                          +MnwithinVarJL_4_C+MnwithinVarJL_4_D)/((26+29+28+29)-4))
  
  within_group_varNi <- ((NiwithinVarJL_4_A+NiwithinVarJL_4_B
                          +NiwithinVarJL_4_C+NiwithinVarJL_4_D)/((26+29+28+29)-4))
  
  #Calculate the F Statistic
  F_Fe <- FeVar/within_group_varFe
  
  F_Ni <- NiVar/within_group_varNi
  
  F_Th <- ThVar/within_group_varTh
  
  F_Nd <- NdVar/within_group_varNd
  
  F_Y <- YVar/within_group_varY
  
  F_Ce <- CeVar/within_group_varCe
  
  F_Ca <- CaVar/within_group_varCa
  
  F_Zr <- ZrVar/within_group_varZr
  
  F_S <- SVar/within_group_varS
  
  F_Mn <- MnVar/within_group_varMn
  
  F_La <- LaVar/within_group_varLa
  
  #Calculate first degree of freedom
  df1 <- 4-1
  
  #Calculate the second degree of freedom
  df2 <- (26+29+28+29)-4
  
  #Calculate the p-value
  pFe <- pf(F_Fe,df1,df2, lower.tail = FALSE)
  pMn <- pf(F_Mn,df1,df2, lower.tail = FALSE)
  pS <- pf(F_S,df1,df2, lower.tail = FALSE)
  pNi <- pf(F_Ni,df1,df2, lower.tail = FALSE)
  pCa <- pf(F_Ca,df1,df2, lower.tail = FALSE)
  pCe <- pf(F_Ce,df1,df2, lower.tail = FALSE)
  pTh <- pf(F_Th,df1,df2, lower.tail = FALSE)
  pZr <- pf(F_Zr,df1,df2, lower.tail = FALSE)
  pY <- pf(F_Y,df1,df2, lower.tail = FALSE)
  pNd <- pf(F_Nd,df1,df2, lower.tail = FALSE)
  pLa <- pf(F_La,df1,df2, lower.tail = FALSE)
  
  
  ##LS4 F Statistic Calculations
  
  #Calculate the variance for Fe, Mn, S, Zr, Ce, Ca, Nd, Hf, Th,La, Ni, and Y
  FeVar <- (28* (round_mean_FBase3Fe-round_mean_Fe)^2
            + 29 * (round_mean_FBase2Fe-round_mean_Fe)^2 
            + 26 * (round_mean_FBase1Fe-round_mean_Fe)^2)/11
  
  ThVar <- (28* (round_mean_FBase3Th-round_mean_Th)^2
            + 29 * (round_mean_FBase2Th-round_mean_Th)^2 
            + 26 * (round_mean_FBase1Th-round_mean_Th)^2)/11
  
  NiVar <- (28* (round_mean_FBase3Ni-round_mean_Ni)^2
            + 29 * (round_mean_FBase2Ni-round_mean_Ni)^2 
            + 26 * (round_mean_FBase1Ni-round_mean_Ni)^2)/11
  
  MnVar <- (28* (round_mean_FBase3Mn-round_mean_Mn)^2
            + 29 * (round_mean_FBase2Mn-round_mean_Mn)^2 
            + 26 * (round_mean_FBase1Mn-round_mean_Mn)^2)/11
  
  NdVar <- (28* (round_mean_FBase3Nd-round_mean_Nd)^2
            + 29 * (round_mean_FBase2Nd-round_mean_Nd)^2 
            + 26 * (round_mean_FBase1Nd-round_mean_Nd)^2)/11
  
  ZrVar <- (28* (round_mean_FBase3Zr-round_mean_Zr)^2
            + 29 * (round_mean_FBase2Zr-round_mean_Zr)^2 
            + 26 * (round_mean_FBase1Zr-round_mean_Zr)^2)/11
  
  YVar <- (28* (round_mean_FBase3Y-round_mean_Y)^2
           + 29 * (round_mean_FBase2Y-round_mean_Y)^2 
           + 26 * (round_mean_FBase1Y-round_mean_Y)^2)/11
  
  CeVar <- (28* (round_mean_FBase3Ce-round_mean_Ce)^2
            + 29 * (round_mean_FBase2Ce-round_mean_Ce)^2 
            + 26 * (round_mean_FBase1Ce-round_mean_Ce)^2)/11
  
  SVar <- (28* (round_mean_FBase3S-round_mean_S)^2
           + 29 * (round_mean_FBase2S-round_mean_S)^2 
           + 26 * (round_mean_FBase1S-round_mean_S)^2)/11
  
  CaVar <- (28* (round_mean_FBase3Ca-round_mean_Ca)^2
            + 29 * (round_mean_FBase2Ca-round_mean_Ca)^2 
            + 26 * (round_mean_FBase1Ca-round_mean_Ca)^2)/11
  
  LaVar <- (28* (round_mean_FBase3La-round_mean_La)^2
            + 29 * (round_mean_FBase2La-round_mean_La)^2 
            + 26 * (round_mean_FBase1La-round_mean_La)^2)/11
  
  
  ##Calculate the within group variance for Fe, Mn, S, Zr, Ce, Ca, Nd, Hf, Th,La, Ni, and Y
  FewithinVarSST <- sum((SSTSubset$Fe-round_mean_SSTFe)^2) 
  FewithinVarShale <- sum((ShaleSubset$Fe-round_mean_ShaleFe)^2)
  FewithinVarEBase <- sum((EBaseSubset$Fe-round_mean_EBaseFe)^2) 
  FewithinVarFBase1 <- sum((FBase1Subset$Fe-round_mean_FBase1Fe)^2)
  FewithinVarFBase2 <- sum((FBase2Subset$Fe-round_mean_FBase2Fe)^2)
  FewithinVarFBase3 <- sum((FBase3Subset$Fe-round_mean_FBase3Fe)^2)
  FewithinVarTorrivio <- sum((TorrivioSubset$Fe-round_mean_TorrivioFe)^2)
  FewithinVarB6 <- sum((B6Subset$Fe-round_mean_B6Fe)^2)
  FewithinVarUpperGallup <- sum((UpperGallupSubset$Fe-round_mean_UpperGallupFe)^2)
 
  CewithinVarSST <- sum((SSTSubset$Ce-round_mean_SSTCe)^2) 
  CewithinVarShale <- sum((ShaleSubset$Ce-round_mean_ShaleCe)^2)
  CewithinVarEBase <- sum((EBaseSubset$Ce-round_mean_EBaseCe)^2) 
  CewithinVarFBase1 <- sum((FBase1Subset$Ce-round_mean_FBase1Ce)^2)
  CewithinVarFBase2 <- sum((FBase2Subset$Ce-round_mean_FBase2Ce)^2)
  CewithinVarFBase3 <- sum((FBase3Subset$Ce-round_mean_FBase3Ce)^2)
  CewithinVarTorrivio <- sum((TorrivioSubset$Ce-round_mean_TorrivioCe)^2)
  CewithinVarB6 <- sum((B6Subset$Ce-round_mean_B6Ce)^2)
  CewithinVarUpperGallup <- sum((UpperGallupSubset$Ce-round_mean_UpperGallupCe)^2)
  
  NiwithinVarSST <- sum((SSTSubset$Ni-round_mean_SSTNi)^2) 
  NiwithinVarShale <- sum((ShaleSubset$Ni-round_mean_ShaleNi)^2)
  NiwithinVarEBase <- sum((EBaseSubset$Ni-round_mean_EBaseNi)^2) 
  NiwithinVarFBase1 <- sum((FBase1Subset$Ni-round_mean_FBase1Ni)^2)
  NiwithinVarFBase2 <- sum((FBase2Subset$Ni-round_mean_FBase2Ni)^2)
  NiwithinVarFBase3 <- sum((FBase3Subset$Ni-round_mean_FBase3Ni)^2)
  NiwithinVarTorrivio <- sum((TorrivioSubset$Ni-round_mean_TorrivioNi)^2)
  NiwithinVarB6 <- sum((B6Subset$Ni-round_mean_B6Ni)^2)
  NiwithinVarUpperGallup <- sum((UpperGallupSubset$Ni-round_mean_UpperGallupNi)^2)
  
  SwithinVarSST <- sum((SSTSubset$S-round_mean_SSTS)^2) 
  SwithinVarShale <- sum((ShaleSubset$S-round_mean_ShaleS)^2)
  SwithinVarEBase <- sum((EBaseSubset$S-round_mean_EBaseS)^2) 
  SwithinVarFBase1 <- sum((FBase1Subset$S-round_mean_FBase1S)^2)
  SwithinVarFBase2 <- sum((FBase2Subset$S-round_mean_FBase2S)^2)
  SwithinVarFBase3 <- sum((FBase3Subset$S-round_mean_FBase3S)^2)
  SwithinVarTorrivio <- sum((TorrivioSubset$S-round_mean_TorrivioS)^2)
  SwithinVarB6 <- sum((B6Subset$S-round_mean_B6S)^2)
  SwithinVarUpperGallup <- sum((UpperGallupSubset$S-round_mean_UpperGallupS)^2)
  
  ZrwithinVarSST <- sum((SSTSubset$Zr-round_mean_SSTZr)^2) 
  ZrwithinVarShale <- sum((ShaleSubset$Zr-round_mean_ShaleZr)^2)
  ZrwithinVarEBase <- sum((EBaseSubset$Zr-round_mean_EBaseZr)^2) 
  ZrwithinVarFBase1 <- sum((FBase1Subset$Zr-round_mean_FBase1Zr)^2)
  ZrwithinVarFBase2 <- sum((FBase2Subset$Zr-round_mean_FBase2Zr)^2)
  ZrwithinVarFBase3 <- sum((FBase3Subset$Zr-round_mean_FBase3Zr)^2)
  ZrwithinVarTorrivio <- sum((TorrivioSubset$Zr-round_mean_TorrivioZr)^2)
  ZrwithinVarB6 <- sum((B6Subset$Zr-round_mean_B6Zr)^2)
  ZrwithinVarUpperGallup <- sum((UpperGallupSubset$Zr-round_mean_UpperGallupZr)^2)
  
  CawithinVarSST <- sum((SSTSubset$Ca-round_mean_SSTCa)^2) 
  CawithinVarShale <- sum((ShaleSubset$Ca-round_mean_ShaleCa)^2)
  CawithinVarEBase <- sum((EBaseSubset$Ca-round_mean_EBaseCa)^2) 
  CawithinVarFBase1 <- sum((FBase1Subset$Ca-round_mean_FBase1Ca)^2)
  CawithinVarFBase2 <- sum((FBase2Subset$Ca-round_mean_FBase2Ca)^2)
  CawithinVarFBase3 <- sum((FBase3Subset$Ca-round_mean_FBase3Ca)^2)
  CawithinVarTorrivio <- sum((TorrivioSubset$Ca-round_mean_TorrivioCa)^2)
  CawithinVarB6 <- sum((B6Subset$Ca-round_mean_B6Ca)^2)
  CawithinVarUpperGallup <- sum((UpperGallupSubset$Ca-round_mean_UpperGallupCa)^2)
  
  YwithinVarSST <- sum((SSTSubset$Y-round_mean_SSTY)^2) 
  YwithinVarShale <- sum((ShaleSubset$Y-round_mean_ShaleY)^2)
  YwithinVarEBase <- sum((EBaseSubset$Y-round_mean_EBaseY)^2) 
  YwithinVarFBase1 <- sum((FBase1Subset$Y-round_mean_FBase1Y)^2)
  YwithinVarFBase2 <- sum((FBase2Subset$Y-round_mean_FBase2Y)^2)
  YwithinVarFBase3 <- sum((FBase3Subset$Y-round_mean_FBase3Y)^2)
  YwithinVarTorrivio <- sum((TorrivioSubset$Y-round_mean_TorrivioY)^2)
  YwithinVarB6 <- sum((B6Subset$Y-round_mean_B6Y)^2)
  YwithinVarUpperGallup <- sum((UpperGallupSubset$Y-round_mean_UpperGallupY)^2)
  
  MnwithinVarSST <- sum((SSTSubset$Mn-round_mean_SSTMn)^2) 
  MnwithinVarShale <- sum((ShaleSubset$Mn-round_mean_ShaleMn)^2)
  MnwithinVarEBase <- sum((EBaseSubset$Mn-round_mean_EBaseMn)^2) 
  MnwithinVarFBase1 <- sum((FBase1Subset$Mn-round_mean_FBase1Mn)^2)
  MnwithinVarFBase2 <- sum((FBase2Subset$Mn-round_mean_FBase2Mn)^2)
  MnwithinVarFBase3 <- sum((FBase3Subset$Mn-round_mean_FBase3Mn)^2)
  MnwithinVarTorrivio <- sum((TorrivioSubset$Mn-round_mean_TorrivioMn)^2)
  MnwithinVarB6 <- sum((B6Subset$Mn-round_mean_B6Mn)^2)
  MnwithinVarUpperGallup <- sum((UpperGallupSubset$Mn-round_mean_UpperGallupMn)^2)
  
  NdwithinVarSST <- sum((SSTSubset$Nd-round_mean_SSTNd)^2) 
  NdwithinVarShale <- sum((ShaleSubset$Nd-round_mean_ShaleNd)^2)
  NdwithinVarEBase <- sum((EBaseSubset$Nd-round_mean_EBaseNd)^2) 
  NdwithinVarFBase1 <- sum((FBase1Subset$Nd-round_mean_FBase1Nd)^2)
  NdwithinVarFBase2 <- sum((FBase2Subset$Nd-round_mean_FBase2Nd)^2)
  NdwithinVarFBase3 <- sum((FBase3Subset$Nd-round_mean_FBase3Nd)^2)
  NdwithinVarTorrivio <- sum((TorrivioSubset$Nd-round_mean_TorrivioNd)^2)
  NdwithinVarB6 <- sum((B6Subset$Nd-round_mean_B6Nd)^2)
  NdwithinVarUpperGallup <- sum((UpperGallupSubset$Nd-round_mean_UpperGallupNd)^2)
  
  ThwithinVarSST <- sum((SSTSubset$Th-round_mean_SSTTh)^2) 
  ThwithinVarShale <- sum((ShaleSubset$Th-round_mean_ShaleTh)^2)
  ThwithinVarEBase <- sum((EBaseSubset$Th-round_mean_EBaseTh)^2) 
  ThwithinVarFBase1 <- sum((FBase1Subset$Th-round_mean_FBase1Th)^2)
  ThwithinVarFBase2 <- sum((FBase2Subset$Th-round_mean_FBase2Th)^2)
  ThwithinVarFBase3 <- sum((FBase3Subset$Th-round_mean_FBase3Th)^2)
  ThwithinVarTorrivio <- sum((TorrivioSubset$Th-round_mean_TorrivioTh)^2)
  ThwithinVarB6 <- sum((B6Subset$Th-round_mean_B6Th)^2)
  ThwithinVarUpperGallup <- sum((UpperGallupSubset$Th-round_mean_UpperGallupTh)^2)
  
  LawithinVarSST <- sum((SSTSubset$La-round_mean_SSTLa)^2) 
  LawithinVarShale <- sum((ShaleSubset$La-round_mean_ShaleLa)^2)
  LawithinVarEBase <- sum((EBaseSubset$La-round_mean_EBaseLa)^2) 
  LawithinVarFBase1 <- sum((FBase1Subset$La-round_mean_FBase1La)^2)
  LawithinVarFBase2 <- sum((FBase2Subset$La-round_mean_FBase2La)^2)
  LawithinVarFBase3 <- sum((FBase3Subset$La-round_mean_FBase3La)^2)
  LawithinVarTorrivio <- sum((TorrivioSubset$La-round_mean_TorrivioLa)^2)
  LawithinVarB6 <- sum((B6Subset$La-round_mean_B6La)^2)
  LawithinVarUpperGallup <- sum((UpperGallupSubset$La-round_mean_UpperGallupLa)^2)
  
  #Calculate the between group variance for Fe, Mn, S, Zr, Ce, Ca, Nd, Hf, Th,La, Ni, and Y
  within_group_varFe <- ((FewithinVarFBase1+FewithinVarFBase2+FewithinVarFBase3)/((26+29+28)-3))
  
  within_group_varTh <- ((ThwithinVarFBase1+ThwithinVarFBase2+ThwithinVarFBase3)/((26+29+28)-3))
  
  within_group_varNi <- ((NiwithinVarFBase1+NiwithinVarFBase2+NiwithinVarFBase3)/((26+29+28)-3))
  
  within_group_varNd <- ((NdwithinVarFBase1+NdwithinVarFBase2+NdwithinVarFBase3)/((26+29+28)-3))
  
  within_group_varS <- ((SwithinVarFBase1+SwithinVarFBase2+SwithinVarFBase3)/((26+29+28)-3))
  
  within_group_varMn <- ((MnwithinVarFBase1+MnwithinVarFBase2+MnwithinVarFBase3)/((26+29+28)-3))
  
  within_group_varCe <- ((CewithinVarFBase1+CewithinVarFBase2+CewithinVarFBase3)/((26+29+28)-3))
 
  within_group_varZr <- ((ZrwithinVarFBase1+ZrwithinVarFBase2+ZrwithinVarFBase3)/((26+29+28)-3))
  
  within_group_varCa <- ((CawithinVarFBase1+CawithinVarFBase2+CawithinVarFBase)/((26+29+28)-3))
  
  within_group_varY <- ((YwithinVarFBase1+YwithinVarFBase2+YwithinVarFBase3)/((26+29+28)-3))
  
  within_group_varLa <- ((LawithinVarFBase1+LawithinVarFBase2+LawithinVarFBase3)/((26+29+28)-3))
  
  #Calculate F Statistic
  F_Fe <- FeVar/within_group_varFe
  
  F_Ni <- NiVar/within_group_varNi
  
  F_Th <- ThVar/within_group_varTh
  
  F_Nd <- NdVar/within_group_varNd
  
  F_Y <- YVar/within_group_varY
 
  F_Ce <- CeVar/within_group_varCe
  
  F_Ca <- CaVar/within_group_varCa
  
  F_Zr <- ZrVar/within_group_varZr
  
  F_S <- SVar/within_group_varS
  
  F_Mn <- MnVar/within_group_varMn
  
  F_La <- LaVar/within_group_varLa
  
  #Calculate the first degree of freedom
  df1 <- 3-1
  
  #Calculate teh second degree of freedom
  df2 <- (28+29+26)-3
  
  #Calculate the p-value
  pFe <- pf(F_Fe,df1,df2, lower.tail = FALSE)
  pMn <- pf(F_Mn,df1,df2, lower.tail = FALSE)
  pS <- pf(F_S,df1,df2, lower.tail = FALSE)
  pNi <- pf(F_Ni,df1,df2, lower.tail = FALSE)
  pCa <- pf(F_Ca,df1,df2, lower.tail = FALSE)
  pCe <- pf(F_Ce,df1,df2, lower.tail = FALSE)
  pTh <- pf(F_Th,df1,df2, lower.tail = FALSE)
  pZr <- pf(F_Zr,df1,df2, lower.tail = FALSE)
  pY <- pf(F_Y,df1,df2, lower.tail = FALSE)
  pNd <- pf(F_Nd,df1,df2, lower.tail = FALSE)
  pLa <- pf(F_La,df1,df2, lower.tail = FALSE)
  
  #PCA Code
  
  #Call data from Excel file
  library(readxl)
  Bentonite2 <- read_excel("C:/XRF Bentonite Scan Moly Tube Test 2.xlsx")
  BentonitePCA2 <- Bentonite2[,1:57]
  Bentonite <- read_excel("C:/Bentonites- 3.xlsx")
  
  #Limit the rows and columns to those containing the desired data
  BentonitePCA <- Bentonite[1:243,1:58]
  
  BentonitePCA2_Short <- Bentonite2[,1:10]
  install.packages("tidyverse")
  
  Bentonite.pr <- prcomp(BentonitePCA, center = TRUE, scale = TRUE)
  summary(Bentonite.pr)
  
  Bentonite2.pr <- prcomp(BentonitePCA2, center = TRUE, scale = TRUE)
  summary(Bentonite2.pr)
  
  Bentonite2_Short.pr <- prcomp(BentonitePCA2_Short, center = TRUE, scale = TRUE)
  s <- summary(Bentonite2_Short.pr)
  
  biplot(Bentonite2.pr)
  #Plot individuals
  plot(Bentonite2.pr$x[,1], Bentonite2.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  biplot(Bentonite2_Short.pr)
  #Plot individuals
  plot(Bentonite2_Short.pr$x[,1], Bentonite2_Short.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  # Create groups
  pch.group <- c(rep(21, times=16), rep(22, times=14), rep(24, times=3))
  col.group <- c(rep("skyblue2", times=16), rep("gold", times=14), rep("green2", times=3))
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  # Add labels
  ##Bentonite2_Short
  text(Bentonite2_Short.pr$x[,1], Bentonite2_Short.pr$x[,2], labels=row.names(Bentonite2_Short.pr$x), pos=c(1,3,4,2), font=2)
  
  ##Bentonite2
  text(Bentonite2.pr$x[,1], Bentonite2.pr$x[,2], labels=row.names(Bentonite2.pr$x), pos=c(1,3,4,2), font=2)
  
  ##Total Mancos Shale Bentonite Samples
  
  MS1<- BentonitePCA[,apply(BentonitePCA, 2, var, na.rm=TRUE) != 0]
  MS1.pr<- prcomp(MS1[c(1:20)], center = TRUE, scale = TRUE)
  summary(MS1.pr)
  
  text(MS1.pr$x[,1], MS1.pr$x[,2], labels=row.names(MS1.pr$x), pos=c(1,3,4,2), font=2)
  
  # Create groups for Juana Lopez Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  #Plot individuals
  plot(MS1.pr$x[,1], MS1.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(MS1.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("topright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##SST
  SST1<- SST[,apply(SST, 2, var, na.rm=TRUE) != 0]
  SST.pr<- prcomp(SST1[c(1:20)], center = TRUE, scale = TRUE)
  s <- summary(SST.pr)
  biplot(SST.pr)
  
  # Create groups for Mancos Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  text(SST.pr$x[,1], SST.pr$x[,2], labels=row.names(SST.pr$x), pos=c(1,3,4,2), font=2)
  
  #Plot individuals
  plot(SST.pr$x[,1], SST.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(SST.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("bottomright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##Shale
  Shale1<- Shale[,apply(Shale, 2, var, na.rm=TRUE) != 0]
  Shale.pr<- prcomp(Shale1[c(1:55)], center = TRUE, scale = TRUE)
  summary(Shale.pr)
  
  # Create groups for Mancos Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  text(Shale.pr$x[,1], Shale.pr$x[,2], labels=row.names(Shale.pr$x), pos=c(1,3,4,2), font=2)
  
  #Plot individuals
  plot(Shale.pr$x[,1], Shale.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(Shale.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("topright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##UpperGallup
  UpperGallup1<- UpperGallup[,apply(UpperGallup, 2, var, na.rm=TRUE) != 0]
  UpperGallup.pr<- prcomp(UpperGallup1[c(1:58)], center = TRUE, scale = TRUE)
  summary(UpperGallup.pr)
  
  # Create groups for Mancos Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  text(UpperGallup.pr$x[,1], UpperGallup.pr$x[,2], labels=row.names(UpperGallup.pr$x), pos=c(1,3,4,2), font=2)
  
  #Plot individuals
  plot(UpperGallup.pr$x[,1], UpperGallup.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  text(l.x, l.y, labels=row.names(UpperGallup.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("bottomright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##B6
  B61<- B6[,apply(B6, 2, var, na.rm=TRUE) != 0]
  B6.pr<- prcomp(B61[c(1:57)], center = TRUE, scale = TRUE)
  summary(B6.pr)
  
  text(B6.pr$x[,1], B6.pr$x[,2], labels=row.names(B6.pr$x), pos=c(1,3,4,2), font=2)
  
  # Create groups
  pch.group <- c(rep(21, times=8), rep(22, times=10), rep(24, times=30))
  col.group <- c(rep("skyblue2", times=8), rep("gold", times=10), rep("green2", times=30))
  
  #Plot individuals
  plot(B6.pr$x[,1], B6.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(B6.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("topright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##F BASE 1
  FBase1_1<- FBase1[,apply(FBase1, 2, var, na.rm=TRUE) != 0]
  FBase1.pr<- prcomp(FBase1_1[c(1:57)], center = TRUE, scale = TRUE)
  summary(FBase1.pr)
  
  # Create groups for Mancos Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  text(FBase1.pr$x[,1], FBase1.pr$x[,2], labels=row.names(FBase1.pr$x), pos=c(1,3,4,2), font=2)
  
  #Plot individuals
  plot(FBase1.pr$x[,1], FBase1.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(FBase1.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("topright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##F BASE 2
  FBase2_1<- FBase2[,apply(FBase2, 2, var, na.rm=TRUE) != 0]
  FBase2.pr<- prcomp(FBase2_1[c(1:57)], center = TRUE, scale = TRUE)
  summary(FBase2.pr)
  
  # Create groups for Mancos Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  text(FBase2.pr$x[,1], FBase2.pr$x[,2], labels=row.names(FBase2.pr$x), pos=c(1,3,4,2), font=2)
  
  #Plot individuals
  plot(FBase2.pr$x[,1], FBase2.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(FBase2.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("topright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##F BASE 3
  FBase3_1<- FBase3[,apply(FBase3, 2, var, na.rm=TRUE) != 0]
  FBase3.pr<- prcomp(FBase3_1[c(1:56)], center = TRUE, scale = TRUE)
  summary(FBase3.pr)
  
  # Create groups for Mancos Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  text(FBase3.pr$x[,1], FBase3.pr$x[,2], labels=row.names(FBase3.pr$x), pos=c(1,3,4,2), font=2)
  
  #Plot individuals
  plot(FBase3.pr$x[,1], FBase3.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(FBase3.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("topright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##E BASE
  EBase1<- EBase[,apply(EBase, 2, var, na.rm=TRUE) != 0]
  EBase.pr<- prcomp(EBase1[c(1:57)], center = TRUE, scale = TRUE)
  summary(EBase.pr)
  
  text(EBase.pr$x[,1], EBase.pr$x[,2], labels=row.names(EBase.pr$x), pos=c(1,3,4,2), font=2)
  
  # Create groups for Mancos Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  #Plot individuals
  plot(EBase.pr$x[,1], EBase.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(EBase.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("topright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##Torrivio
  Torrivio1<- Torrivio[,apply(Torrivio, 2, var, na.rm=TRUE) != 0]
  Torrivio.pr<- prcomp(Torrivio1[c(1:55)], center = TRUE, scale = TRUE)
  summary(Torrivio.pr)
  
  # Create groups for Mancos Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  text(Torrivio.pr$x[,1], Torrivio.pr$x[,2], labels=row.names(Torrivio.pr$x), pos=c(1,3,4,2), font=2)
  
  #Plot individuals
  plot(Torrivio.pr$x[,1], Torrivio.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(Torrivio.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("topleft", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##Juana Lopez Samples
  
  ##Total Juana Lopez Samples
  JL1<- BentonitePCA2[,apply(BentonitePCA2, 2, var, na.rm=TRUE) != 0]
  JL1.pr<- prcomp(JL1[c(1:20)], center = TRUE, scale = TRUE)
  summary(JL1.pr)
  
  text(JL1.pr$x[,1], JL1.pr$x[,2], labels=row.names(JL1.pr$x), pos=c(1,3,4,2), font=2)
  
  # Create groups for Juana Lopez Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  #Plot individuals
  plot(JL1.pr$x[,1], JL1.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(JL1.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("topright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  #JL_4_A
  JL_4_A1<- JL_4_A[,apply(JL_4_A, 2, var, na.rm=TRUE) != 0]
  JL_4_A.pr<- prcomp(JL_4_A1[c(1:56)], center = TRUE, scale = TRUE)
  summary(JL_4_A.pr)
  
  text(JL_4_A.pr$x[,1], JL_4_A.pr$x[,2], labels=row.names(JL_4_A.pr$x), pos=c(1,3,4,2), font=2)
  
  # Create groups for Juana Lopez Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  #Plot individuals
  plot(JL_4_A.pr$x[,1], JL_4_A.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(JL_4_A.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("bottomleft", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##JL_4_B
  JL_4_B1<- JL_4_B[,apply(JL_4_B, 2, var, na.rm=TRUE) != 0]
  JL_4_B.pr<- prcomp(JL_4_B1[c(1:56)], center = TRUE, scale = TRUE)
  summary(JL_4_B.pr)
  
  text(JL_4_B.pr$x[,1], JL_4_B.pr$x[,2], labels=row.names(JL_4_B.pr$x), pos=c(1,3,4,2), font=2)
  
  # Create groups for Juana Lopez Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  #Plot individuals
  plot(JL_4_B.pr$x[,1], JL_4_B.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(JL_4_B.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("bottomright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##JL_4_C
  JL_4_C1<- JL_4_C[,apply(JL_4_C, 2, var, na.rm=TRUE) != 0]
  JL_4_C.pr<- prcomp(JL_4_C1[c(1:56)], center = TRUE, scale = TRUE)
  summary(JL_4_C.pr)
  
  text(JL_4_C.pr$x[,1], JL_4_C.pr$x[,2], labels=row.names(JL_4_C.pr$x), pos=c(1,3,4,2), font=2)
  
  # Create groups for Juana Lopez Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  #Plot individuals
  plot(JL_4_C.pr$x[,1], JL_4_C.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(JL_4_C.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("topright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##JL_4_D
  JL_4_D1<- JL_4_D[,apply(JL_4_D, 2, var, na.rm=TRUE) != 0]
  JL_4_D.pr<- prcomp(JL_4_D1[c(1:56)], center = TRUE, scale = TRUE)
  summary(JL_4_D.pr)
  
  text(JL_4_D.pr$x[,1], JL_4_D.pr$x[,2], labels=row.names(JL_4_D.pr$x), pos=c(1,3,4,2), font=2)
  
  # Create groups for Juana Lopez Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  #Plot individuals
  plot(JL_4_D.pr$x[,1], JL_4_D.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(JL_4_D.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("bottomleft", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##JL_1
  JL_11<- JL_1[,apply(JL_1, 2, var, na.rm=TRUE) != 0]
  JL_1.pr<- prcomp(JL_11[c(1:56)], center = TRUE, scale = TRUE)
  summary(JL_1.pr)
  
  text(JL_1.pr$x[,1], JL_1.pr$x[,2], labels=row.names(JL_1.pr$x), pos=c(1,3,4,2), font=2)
  
  # Create groups for Juana Lopez Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  #Plot individuals
  plot(JL_1.pr$x[,1], JL_1.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(JL_1.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("bottomright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##JL_2
  JL_21<- JL_2[,apply(JL_2, 2, var, na.rm=TRUE) != 0]
  JL_2.pr<- prcomp(JL_21[c(1:56)], center = TRUE, scale = TRUE)
  summary(JL_2.pr)
  
  text(JL_2.pr$x[,1], JL_2.pr$x[,2], labels=row.names(JL_2.pr$x), pos=c(1,3,4,2), font=2)
  
  # Create groups for Juana Lopez Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  #Plot individuals
  plot(JL_2.pr$x[,1], JL_2.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(JL_2.pr$rotation), col="red", pos=l.pos)
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("bottomright", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  ##JL_3
  JL_31<- JL_3[,apply(JL_3, 2, var, na.rm=TRUE) != 0]
  JL_3.pr<- prcomp(JL_31[c(1:56)], center = TRUE, scale = TRUE)
  summary(JL_3.pr)
  biplot(JL_3.pr)
  
  text(JL_3.pr$x[,1], JL_3.pr$x[,2], labels=row.names(JL_3.pr$x), pos=c(1,3,4,2), font=2)
  
  l.x <- JL_3.pr$rotation[,1]*5
  l.y <- JL_3.pr$rotation[,2]*5
  
  # Create groups for Juana Lopez Bentonite Samples
  pch.group <- c(rep(21, times=5), rep(22, times=7), rep(24, times=11), rep (23, times=33))
  col.group <- c(rep("skyblue2", times=5), rep("gold", times=7), rep("green2", times=11), rep("purple2", times=33))
  
  #Plot individuals
  plot(JL_3.pr$x[,1], JL_3.pr$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)
  
  text(l.x, l.y, labels=row.names(JL_3.pr$rotation), col="red", pos=l.pos, main= )
  
  # Add grid lines
  abline(v=0, lty=2, col="grey50")
  abline(h=0, lty=2, col="grey50")
  
  # Add legend
  legend("topleft", legend=c("Terrigenous", "Redox", "REE", "Miscellaneous"), col="black", pt.bg=c("skyblue2", "gold", "green2","purple2"), pch=c(21, 22, 24, 23), pt.cex=1.5)
  
  # Draw arrows
  arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red", length=2, lwd=1.5)
  
  library(psych)
  # Get co-ordinates of variables (loadings), and multiply by 5
  l.x <- Bentonite2_Short.pr$rotation[,1]*5
  l.y <- Bentonite2_Short.pr$rotation[,2]*5
  
  l.x <- Bentonite2.pr$rotation[,1]*5
  l.y <- Bentonite2.pr$rotation[,2]*5
  
  l.x <- B6.pr$rotation[,1]*5
  l.y <- B6.pr$rotation[,2]*5
  
  l.x <- SST.pr$rotation[,1]*5
  l.y <- SST.pr$rotation[,2]*5
  
  l.x <- Shale.pr$rotation[,1]*5
  l.y <- Shale.pr$rotation[,2]*5
  
  l.x <- UpperGallup.pr$rotation[,1]*5
  l.y <- UpperGallup.pr$rotation[,2]*5
  
  l.x <- Torrivio.pr$rotation[,1]*5
  l.y <- Torrivio.pr$rotation[,2]*5
  
  l.x <- EBase.pr$rotation[,1]*5
  l.y <- EBase.pr$rotation[,2]*5
  
  l.x <- FBase1.pr$rotation[,1]*5
  l.y <- FBase1.pr$rotation[,2]*5
  
  l.x <- FBase2.pr$rotation[,1]*5
  l.y <- FBase2.pr$rotation[,2]*5
  
  l.x <- FBase3.pr$rotation[,1]*5
  l.y <- FBase3.pr$rotation[,2]*5
  
  l.x <- JL_1.pr$rotation[,1]*5
  l.y <- JL_1.pr$rotation[,2]*5
  
  l.x <- JL_2.pr$rotation[,1]*5
  l.y <- JL_2.pr$rotation[,2]*5
  
  l.x <- JL_4A.pr$rotation[,1]*5
  l.y <- JL_4A.pr$rotation[,2]*5
  
  l.x <- JL_4B.pr$rotation[,1]*5
  l.y <- JL_4B.pr$rotation[,2]*5
  
  l.x <- JL_4C.pr$rotation[,1]*5
  l.y <- JL_4C.pr$rotation[,2]*5
  
  l.x <- JL_4D.pr$rotation[,1]*5
  l.y <- JL_4D.pr$rotation[,2]*5
  
  # Draw arrows code
  arrows(x0=0, x1=l.x, y0=0, y1=l.y, col="red", length=5, lwd=1.5)
  
  # Label position
  l.pos <- l.y # Create a vector of y axis coordinates
  lo <- which(l.y < 0) # Get the variables on the bottom half of the plot
  hi <- which(l.y > 0) # Get variables on the top half
  # Replace values in the vector
  l.pos <- replace(l.pos, lo, "1")
  l.pos <- replace(l.pos, hi, "3")
  
  # Variable labels
  text(l.x, l.y, labels=row.names(Bentonite2_Short.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(Bentonite2.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(B6.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(SST.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(Shale.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(UpperGallup.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(FBase1_1.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(FBase2_1.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(FBase3_1.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(EBase1.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(Torrivio1.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(JL_31.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(JL_21.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(JL_11.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(JL_4_A1.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(JL_4_B1.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(JL_4_C1.pr$rotation), col="red", pos=l.pos)
  
  text(l.x, l.y, labels=row.names(JL_4_D1.pr$rotation), col="red", pos=l.pos)
  
  # Add legend
  legend("topleft", legend=c("Terrigenous", "Redox", "REE"), col="black", pt.bg=c("skyblue2", "gold", "green2"), pch=c(21, 22, 24), pt.cex=1.5)
  
  ##Scree plot code
  
  #Total Juana Lopez Bentonite Data
  screeplot(Bentonite2.pr, type = "l", npcs = 58, main = "Screeplot of the 58 PCs for the Juana Lopez Bentonites")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(Bentonite2.pr$sdev^2 / sum(Bentonite2.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
  abline(v = 27, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC27"),
         col=c("blue"), lty=5, cex=0.6)
  biplot(Bentonite2.pr)
  
  #Total Mancos Shale Bentonite Data
  screeplot(Bentonite.pr, type = "l", npcs = 58, main = "Screeplot of the 58 PCs for Mancos Shale Bentonites")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(Bentonite.pr$sdev^2 / sum(Bentonite.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")
  abline(v = 19, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC19"),
         col=c("blue"), lty=5, cex=0.6)
  
  #Sandstone Data
  screeplot(SST.pr, type = "l", npcs = 58, main = "Screeplot of the 58 PCs for Sandstone Reference")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(SST.pr$sdev^2 / sum(SST.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "Sandstone Reference")
  abline(v = 8, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC8"),
         col=c("blue"), lty=5, cex=0.6)
  
  #Shale Data
  screeplot(Shale.pr, type = "l", npcs = 58, main = "Screeplot of the 58 PCs for Shale Reference")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(Shale.pr$sdev^2 / sum(Shale.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "Shale Reference")
  abline(v = 13, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC13"),
         col=c("blue"), lty=5, cex=0.6)
  
  #Torrivio Data
  screeplot(Torrivio.pr, type = "l", npcs = 58, main = "Screeplot of the 58 PCs for Torrivio Sandstone Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(Torrivio.pr$sdev^2 / sum(Torrivio.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "Torrivio Sandstone Bentonite")
  abline(v = 17, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC17"),
         col=c("blue"), lty=5, cex=0.6)
  
  #Upper Gallup Data
  screeplot(UpperGallup.pr, type = "l", npcs = 58, main = "Upper Gallup Sandstone Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(UpperGallup.pr$sdev^2 / sum(UpperGallup.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "Upper Gallup Sandstone Bentonite")
  abline(v = 17, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC17"),
         col=c("blue"), lty=5, cex=0.6)
  
  #JL-4C Data
  screeplot(B6.pr, type = "l", npcs = 58, main = "JL_3C Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(B6.pr$sdev^2 / sum(B6.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "JL_4C Bentonite")
  abline(v = 14, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC14"),
         col=c("blue"), lty=5, cex=0.6)
  biplot(Bentonite2.pr)
  
  #LS5 1 Data
  screeplot(FBase1.pr, type = "l", npcs = 58, main = "LS5 1 Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(FBase1.pr$sdev^2 / sum(FBase1.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "LS5 1 Bentonite")
  abline(v = 14, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC14"),
         col=c("blue"), lty=5, cex=0.6)
  
  #LS5 2 Data
  screeplot(FBase2.pr, type = "l", npcs = 58, main = "LS5 2 Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(FBase2.pr$sdev^2 / sum(FBase2.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "LS5 2 Bentonite")
  abline(v = 15, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC15"),
         col=c("blue"), lty=5, cex=0.6)
  
  #LS5 3 Data
  screeplot(FBase3.pr, type = "l", npcs = 58, main = "LS5 3 Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(FBase3.pr$sdev^2 / sum(FBase3.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "LS5 3 Bentonite")
  abline(v = 15, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC15"),
         col=c("blue"), lty=5, cex=0.6)
  
  #LS4 Data
  screeplot(EBase.pr, type = "l", npcs = 58, main = "LS4 Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(EBase.pr$sdev^2 / sum(EBase.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "LS4 Bentonite")
  abline(v = 14, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC14"),
         col=c("blue"), lty=5, cex=0.6)
  
  #JL-1 Data
  screeplot(JL_1.pr, type = "l", npcs = 58, main = "JL_1 Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(JL_1.pr$sdev^2 / sum(JL_1.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "JL_1 Bentonite")
  abline(v = 15, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC15"),
         col=c("blue"), lty=5, cex=0.6)
  
  #JL-2 Data
  screeplot(JL_2.pr, type = "l", npcs = 58, main = "JL_2 Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(JL_2.pr$sdev^2 / sum(JL_2.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "JL_2 Bentonite")
  abline(v = 14, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC14"),
         col=c("blue"), lty=5, cex=0.6)
  
  #JL-3 Data
  screeplot(JL_3.pr, type = "l", npcs = 58, main = "JL_3 Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(JL_3.pr$sdev^2 / sum(JL_3.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "JL_3 Bentonite")
  abline(v = 15, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC15"),
         col=c("blue"), lty=5, cex=0.6)
  
  #JL-4A Data
  screeplot(JL_4A.pr, type = "l", npcs = 58, main = "JL_4A Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(JL_4A.pr$sdev^2 / sum(JL_4A.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "JL_4A Bentonite")
  abline(v = 15, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC15"),
         col=c("blue"), lty=5, cex=0.6)
  
  #JL-4B Data
  screeplot(JL_4B.pr, type = "l", npcs = 58, main = "JL_4B Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(JL_4B.pr$sdev^2 / sum(JL_4B.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "JL_4B Bentonite")
  abline(v = 15, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC15"),
         col=c("blue"), lty=5, cex=0.6)
  
  #JL-4C Data
  screeplot(JL_4C.pr, type = "l", npcs = 58, main = "JL_4C Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(JL_4C.pr$sdev^2 / sum(JL_4C.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "JL_4C Bentonite")
  abline(v = 15, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC15"),
         col=c("blue"), lty=5, cex=0.6)
  
  #JL-4D Data
  screeplot(JL_4D.pr, type = "l", npcs = 58, main = "JL_4D Bentonite")
  abline(h = 1, col="red", lty=5)
  legend("topright", legend=c("Eigenvalue = 1"),
         col=c("red"), lty=5, cex=0.6)
  cumpro <- cumsum(JL_4D.pr$sdev^2 / sum(JL_4D.pr$sdev^2))
  plot(cumpro[0:58], xlab = "PC #", ylab = "Amount of explained variance", main = "JL_4D Bentonite")
  abline(v = 15, col="blue", lty=5)
  abline(h = 0.88759, col="blue", lty=5)
  legend("topleft", legend=c("Cut-off @ PC15"),
         col=c("blue"), lty=5, cex=0.6)
  
  
  