# Libraries
library(prospect) # provide the dataset: a dataframe called babynames
library(future) # provide the dataset: a dataframe called babynames
library(future.apply) # provide the dataset: a dataframe called babynames
#############################################################################
# PROCESS LOPEX DATA #
#############################################################################

library(data.table)
library(pracma)
# source('Library_OptimalFeatures_parallel.R')
gitlab_Rep <- 'https://gitlab.com/jbferet/myshareddata/-/raw/master/LOP/'
# Datasets
dbName <- list('LOPEX_DRY_VAL','LOPEX_FRESH_VAL')
dbName2 <- list('DRY_VAL','FRESH_VAL')
# files available
fileName <- list('DataBioch.txt','ReflectanceData.txt','TransmittanceData.txt')
# download LOPEX data
DataBioch <- Refl <- Tran <- lambda <- list()
i = 0
for (db in dbName){
  DataBioch[[db]] <- fread(paste(gitlab_Rep,db,'/',fileName[[1]],sep=''))
  Refl[[db]] <- fread(paste(gitlab_Rep,db,'/',fileName[[2]],sep=''))
  Tran[[db]] <- fread(paste(gitlab_Rep,db,'/',fileName[[3]],sep=''))
  # extract lambda from R and T
  lambda[[db]] <- unlist(Refl[[db]][,1], use.names=FALSE)
  # delete lambda from R and T
  Refl[[db]] <- matrix(unlist(Refl[[db]][,-1], use.names=FALSE),nrow = length(lambda[[db]]))
  Tran[[db]] <- matrix(unlist(Tran[[db]][,-1], use.names=FALSE),nrow = length(lambda[[db]]))
  if (db =='LOPEX_DRY_VAL'){
    Refl[[db]] <- Refl[[db]][,-25]
    Tran[[db]] <- Tran[[db]][,-25]
    DataBioch[[db]] <- DataBioch[[db]][-25,]
  }
}

# Merge LOPEX_DRY_CAL and LOPEX_FRESH_CAL
Refl_VAL <- cbind(Refl[[1]],Refl[[2]])
Tran_VAL <- cbind(Tran[[1]],Tran[[2]])
DataBioch_VAL <- rbind(DataBioch[[1]],DataBioch[[2]])

# define width of individual spectral domains for the identification of optimal domains (in nm)
spectral_width <- 20
# define spectral domain to be investigated (min and max spectral bands)
spectral_domain <- c(1400,2399)
# define number of spectral domains to identify
nbCPU <- 6
# variables to estimate during inversion
Parms2Estimate <- c('EWT','PROT','CBC','N')
# variable of interest for this particular inversion
Target <- 'PROT'
InitValues <- data.frame(EWT=0.01, PROT=0.001,CBC=0.009, N=1.5)
BiochTruth <- c(DataBioch_VAL[,12])[[1]]

# maximum number of features is defined by spectral_width and spectral domain.
# Here, max = 50  --> width(spectral_domain)/spectral_width
number_features <- 50
ListRes_PROT <- optimal_features_SFS(Refl = Refl_VAL, Tran = Tran_VAL, lambda = lambda[[1]], BiochTruth = BiochTruth, Target = Target,
                                    Parms2Estimate = Parms2Estimate, InitValues = InitValues, SpecPROSPECT = SpecPROSPECT,
                                    spectral_domain = spectral_domain, spectral_width = spectral_width, number_features = number_features,
                                    PROSPECT_version = 'PRO',nbCPU = nbCPU)

#  save list containing all information about SFS
# - optimal spectral features for each step
# - corresponding RMSE
# - corresponding estimated parameter of interest ("Target") for each additional feature
save(list = 'ListRes_PROT',file = 'Optimal_PROT_VAL.RData')
