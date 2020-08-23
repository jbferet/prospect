#############################################################################
# INVERT FOLLOWING LOPEX DATASETS TO ESTIMATE EWT, PROTEINS AND CBC USING PROSPECT-PRO
# - VALIDATION DRY
# - VALIDATION FRESH
#############################################################################
library(prospect)
library(data.table)
# download LOPEX datasets from gitlab repository
gitlab_Rep <- 'https://gitlab.com/jbferet/myshareddata/-/raw/master/LOP/'
# Datasets
dbName <- list('LOPEX_DRY_VAL','LOPEX_FRESH_VAL')
dbName2 <- list('DRY_VAL','FRESH_VAL')
# files available
fileName <- list('DataBioch.txt','ReflectanceData.txt','TransmittanceData.txt')
# download LOPEX data
DataBioch <- Refl <- Tran <- list()
i = 0
for (db in dbName){
  DataBioch[[db]] <- fread(paste(gitlab_Rep,db,'/',fileName[[1]],sep=''))
  Refl[[db]] <- fread(paste(gitlab_Rep,db,'/',fileName[[2]],sep=''))
  Tran[[db]] <- fread(paste(gitlab_Rep,db,'/',fileName[[3]],sep=''))
}

for (db in dbName){
  DataBioch[[db]]$PROT <- DataBioch[[db]]$Proteins_Belgium
  DataBioch[[db]]$CBC <- DataBioch[[db]]$LMA-DataBioch[[db]]$Proteins_Belgium
}

# output directory
Results_Dir <- 'RESULTS'
dir.create(Results_Dir,showWarnings = FALSE)

# invert PROSPECT-PRO for the estimation of EWT, Proteins and CBC based on
# optimal spectral domain identified by Feret et al. (RSE 2020)
Parms2Estimate  <- c('EWT','PROT','CBC','N')
Estimated_LeafChem <- list()

# OptDomain <- c(seq(1700,1799),seq(1860,1879),seq(2020,2399))
load('Optimal_CBC_VAL.RData')
MinRMSE <- which(ListRes_CBC$RMSE == min(ListRes_CBC$RMSE))
OptDomain <- c()
for (i in 1:MinRMSE){
  OptDomain <- c(OptDomain,ListRes_CBC$SpectralFeatures_List[[i]])
}
OptDomain = sort(OptDomain)

for (db in dbName){
  Estimated_LeafChem[[db]] <- list()
  lambda <- unlist(Refl[[db]][,1], use.names=FALSE)
  Refl_tmp <- matrix(unlist(Refl[[db]][,-1], use.names=FALSE),nrow = length(lambda))
  Tran_tmp <- matrix(unlist(Tran[[db]][,-1], use.names=FALSE),nrow = length(lambda))
  # Fit spectral data to match PROSPECT with user optical properties
  SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl_tmp,Tran=Tran_tmp,UserDomain = OptDomain,UL_Bounds = FALSE)
  SubSpecPROSPECT <- SubData$SpecPROSPECT
  Sublambda <- SubData$lambda
  SubRefl <- SubData$Refl
  SubTran <- SubData$Tran
  nbsamples <- ncol(SubRefl)
  # Invert PROSPECT with optimal spectral information
  Estimated_LeafChem[[db]]$CBC = matrix(NA,ncol = 1,nrow = nbsamples)
  for (i in 1:nbsamples){
    print(i)
    res <- Invert_PROSPECT(SubSpecPROSPECT,Refl = SubRefl[,i],Tran = SubTran[,i],PROSPECT_version = 'PRO',Parms2Estimate = Parms2Estimate)
    # only save results for variable of interest
    Estimated_LeafChem[[db]]$CBC[i,1] <- res$CBC
  }
  # save results in a text file
  fileName <- file.path(Results_Dir,paste('CBC_',db,'.txt',sep = ''))
  write.table(x = Estimated_LeafChem[[db]]$CBC,file = fileName,quote = FALSE,row.names = FALSE,col.names = FALSE)
}

load('Optimal_PROT_VAL.RData')
MinRMSE <- which(ListRes_PROT$RMSE == min(ListRes_PROT$RMSE))
OptDomain <- c()
for (i in 1:MinRMSE){
  OptDomain <- c(OptDomain,ListRes_PROT$SpectralFeatures_List[[i]])
}
OptDomain = sort(OptDomain)
for (db in dbName){
  lambda <- unlist(Refl[[db]][,1], use.names=FALSE)
  Refl_tmp <- matrix(unlist(Refl[[db]][,-1], use.names=FALSE),nrow = length(lambda))
  Tran_tmp <- matrix(unlist(Tran[[db]][,-1], use.names=FALSE),nrow = length(lambda))
  # Fit spectral data to match PROSPECT with user optical properties
  SubData <- FitSpectralData(SpecPROSPECT=SpecPROSPECT,lambda=lambda,Refl=Refl_tmp,Tran=Tran_tmp,UserDomain = OptDomain,UL_Bounds = FALSE)
  SubSpecPROSPECT <- SubData$SpecPROSPECT
  Sublambda <- SubData$lambda
  SubRefl <- SubData$Refl
  SubTran <- SubData$Tran
  nbsamples <- ncol(SubRefl)
  # Invert PROSPECT with optimal spectral information
  Estimated_LeafChem[[db]]$PROT = matrix(NA,ncol = 1,nrow = nbsamples)
  for (i in 1:nbsamples){
    print(i)
    res <- Invert_PROSPECT(SubSpecPROSPECT,Refl = SubRefl[,i],Tran = SubTran[,i],PROSPECT_version = 'PRO',Parms2Estimate = Parms2Estimate)
    # only save results for variable of interest
    Estimated_LeafChem[[db]]$PROT[i,1] <- res$PROT
  }
  # save results in a text file
  fileName <- file.path(Results_Dir,paste('PROT_',db,'.txt',sep = ''))
  write.table(x = Estimated_LeafChem[[db]]$PROT,file = fileName,quote = FALSE,row.names = FALSE,col.names = FALSE)
}

