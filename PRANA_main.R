library(dnapath) # To obtain mutual information (MI) estimate via ARACNE
library(dplyr) 
library(parallel) # To use mclapply() when reestimating the association matrix.
library(robustbase) # To fit a robust regression

dir_main = "your file directory where you saved three R codes."
dir_COPD_RGO_data = "your file directory where you saved RDS file from the repository."

source(file.path(dir_main, "EmpiricalBayes_Datta_2005.R"))  # This will load the code calculating adjusted p-values for each genes.
source(file.path(dir_main, "TotalConnectivity.R")) # This is to calculate the total connectivity (thetahats)

combinedCOPDdat_RGO = readRDS(file.path(dir_COPD_RGO_data, "combinedCOPD_RelatedGenesOnly.rds")) # Combined phenotype and expression data


est_method = run_aracne # ARACNE
rnaseqdat = combinedCOPDdat_RGO[ , 8:ncol(combinedCOPDdat_RGO)] # gene expression data part of the downloaded data.
rnaseqdat = as.data.frame(apply(rnaseqdat, 2, as.numeric)) 

## Additional covariates (phenotypic data) sorted by current smoking groups:
phenodat = combinedCOPDdat_RGO[order(combinedCOPDdat_RGO$currentsmoking), 2:7] # first column is ID, so not using it.


#######################################
#                                     #
#    Apply PRANA to the COPD data     #
#                                     #
#######################################

#############################################################################################
#  STEP 1. Estimate an association matrix via ARACNE from the RNA-seq expression data.
#############################################################################################
newindex_A = which(combinedCOPDdat_RGO$currentsmoking == 0) # not current smoker (namely Group A)
newindex_B = which(combinedCOPDdat_RGO$currentsmoking == 1) # current smoker (namely Group B)
rnaseqdatA = rnaseqdat[newindex_A, ] 
rnaseqdatB = rnaseqdat[newindex_B, ]
nw_est_grpA = est_method(rnaseqdatA, verbose = F) # Association matrix for Group A
nw_est_grpB = est_method(rnaseqdatB, verbose = F) # Association matrix for Group B

n_A <- length(newindex_A) # sample size for Group A
n_B <- length(newindex_B) # sample size for Group B


#############################################################################################
#  STEP 2. Calculate \hat{\theta}_{k} (total connectivity) by taking the column sum of the
#          association matrix to obtain the total connectivity of for each gene.
#############################################################################################
thetahat_grpA = thetahats(nw_est_grpA) 
thetahat_grpB = thetahats(nw_est_grpB)


#############################################################################################
#  STEP 3. Re-estimate association matrix using the expression data without i-th subject.
#          Then, calculate \hat{\theta}_{k(i)} from the reestimated association matrix.
#############################################################################################
# Re-estimation part
nw_est_drop_grpA <- mclapply(newindex_A, function(j) est_method(rnaseqdatA[-j, ], verbose = F))
nw_est_drop_grpB <- mclapply(newindex_B, function(j) est_method(rnaseqdatB[-j, ], verbose = F))
# thetahat_{-i} for each gene
thetahat_drop_grpA <- sapply(nw_est_drop_grpA, thetahats)
thetahat_drop_grpB <- sapply(nw_est_drop_grpB, thetahats)



#############################################################################################
#  STEP 4. Calculate jackknife pseudovalues (\tilde{\theta}_{ik} using \hat{\theta}_{k} 
#          and \hat{\theta}_{k(i)}.
#############################################################################################
thetatildefun <- function(thetahatinput, thetahatdropinput, sizegroup) {
        thetatildeout = matrix(NA, ncol=length(thetahatinput), nrow=sizegroup)
        thetatildeout = sapply(1:nrow(thetahatdropinput), function(k) {
                sizegroup * thetahatinput[k] - (sizegroup - 1) * thetahatdropinput[k, ]
        })
        return(thetatildeout)
}

thetatilde_grpA = thetatildefun(thetahat_grpA, thetahat_drop_grpA, n_A)
thetatilde_grpB = thetatildefun(thetahat_grpB, thetahat_drop_grpB, n_B)
thetatilde = rbind(thetatilde_grpA, thetatilde_grpB)
colnames(thetatilde) = colnames(rnaseqdat) # Map the column names (gene names) 


#############################################################################################
#  STEP 5. Fit a robust regression model
##############################################################################################
pseudo.beta_list <- lapply(1:ncol(thetatilde), function(i) {
        m <- thetatilde[, i]
        df <- data.frame(phenodat,
                         m = m)
        fit <- ltsReg(m ~ currentsmoking + packyrs + age + gender + race + FEV1perc, data = df, mcd=FALSE) # include a set of covariates in this model
        return(fit)
})

### Obtain p-values for each genes:
beta_hat = vector(mode = "list", ncol(thetatilde)) 
p_values = vector(mode = "list", ncol(thetatilde)) 
k = NULL
for(k in 1:ncol(thetatilde)) {
        # Estimates for the beta coefficients:
        beta_hat[[k]] <- summary(pseudo.beta_list[[k]])$coef[-1, "Estimate"]
        p_values[[k]] <- summary(pseudo.beta_list[[k]])$coef[-1, "Pr(>|t|)"]
        
}

beta_hat = as.data.frame(bind_rows(beta_hat))  # Convert list into data.frame
rownames(beta_hat) <- colnames(rnaseqdat) # Map the gene names to the data.frame for betahats

p_values = as.data.frame(bind_rows(p_values))  # Convert list into data.frame
rownames(p_values) <- colnames(rnaseqdat) # Map the gene names to the data.frame for p-values
current_smoke_pval = p_values[, 1] # p-values for current smoking status (binary group variable)

# Compute the adjusted p-values via empirical Bayes approach proposed by the reference below:
# NOTE: EBS() is code from Datta S and Datta S (2005). Empirical Bayes screening of many p-values with applications to microarray studies. Bioinformatics, 21(9), 1987-1994
adjp_values = EBS(pvo = current_smoke_pval, alpha = 0.05, B = 500, h = 1)
names(adjp_values) <- colnames(rnaseqdat) # Map the gene names to the data.frame for adj p-values


# Lastly, return the gene names of the significantly differentially connected (DC) genes from PRANA at the 0.05 significance level.
sigDCpseudo = adjp_values[which(adjp_values < 0.05)]
names(sigDCpseudo) 

