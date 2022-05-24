## Below is a main code for simulation scenario 3.
library(SeqNet)
library(dnapath)
library(parallel)
library(robustbase)

RNGkind(sample.kind = "Rounding")

data = "/Users/seungjunahn/Desktop/University of Florida/RA/Pseudovalue/Simulation/Code_Recent_03312021/data"
dir_COPDdata = "/Users/seungjunahn/Desktop/University of Florida/RA/Pseudovalue/COPD_Analysis/Data"

dir_sim = "/Users/seungjunahn/Desktop/PRANA_GitHub/Simulation"
AgeEffectNetworks = "/Users/seungjunahn/Desktop/University of Florida/RA/Pseudovalue/Simulation/Code_05052021/Multivariable/AgeEffectNetworks/PertubAE_Network_WIthinLoop/AssumeAE_03212022_Datta/2hubs/"

genesize = 20 # This will be defined when networks are created.
desired.totalN = 40  # Must be an even number; desired number of samples we want to simulate (maximum 1000)
N = 1000 # Total maximum number of samples we want to simulate
simnum = 1 # Simulation size
cores = 4 # This to use the parallel computing (mclapply)

##################################################################################
#  STEP 1. Generate random network and assign a partial correlation to edges 
#          to obtain weighted networks
##################################################################################

# Load in the topological measures
source(file.path(dir_sim, "TotalConnectivity.R"))

# Load the data from the COPDGene study:
combinedCOPDdat = readRDS(file.path(dir_COPDdata, "combinedCOPDdat_validation.rds"))
finalCOPDdat_expr = combinedCOPDdat[ , 9:ncol(combinedCOPDdat)] # gene expression data
finalCOPDdat_expr = as.data.frame(apply(finalCOPDdat_expr, 2, as.numeric)) 
# 

est_method <- run_aracne

sim_results <- vector("list", simnum)