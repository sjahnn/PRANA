# NOTE: EBS() is from Datta S and Datta S (2005). Empirical Bayes screening of many p-values with applications to microarray studies. Bioinformatics, 21(9), 1987-1994
EBS <- function(pvo,alpha=0.05,B=500,h=1) {
        
        Kh <- function(x)  dnorm(x/h)/h
        ngene <- length(pvo)
        
        #======================================================================================#
        # This part calculates B batches of monotonited Empirical Bayes estimates of the       #
        # location of the z-transformed p-values under the null hypothesis asuming tests       #
        # are independent; h is a user selectable bandwidth; ngene denotes the total number    # 
        # of p-values to be screened.                                                          #  
        # Depending on the power of your computer and ngene = the length of pvo, and B, this   #
        # part may # need several hours. However if you are going to screen multiple samples   #
        # each of same length, this part needs to be computed only once and can be stored for  #
        # future screening.                                                                    #                                         
        #                                                                                      #
        u.n <- u.no <- matrix(0,B,ngene)                                                    #
        for (iii in 1:B)                                                                    #
        {                                                                                 #
                zo.n <- rnorm(ngene)                                                             #
                z.n <- rep(0,ngene)                                                              #
                for (i in 1:ngene) {                                                             #
                        sum1 <- sum(Kh(zo.n[i]-zo.n))                                           #
                        sum2 <- sum((zo.n[i]-zo.n)*Kh(zo.n[i]-zo.n))                            #
                        z.n[i] <- zo.n[i] - (h^-2)*(sum2/sum1)                                  #
                }                                                             #
                u.n[iii,] <- z.n                                                                 #
                for (ii in (ngene-1):1) u.n[iii,ii] <- min(u.n[iii,ii+1],u.n[iii,ii])            #
        }                                                                                 #
        #======================================================================================#
        
        
        
        
        
        # ====================================================================================#
        #We Z-transform the p-values                                            
        
        zo <- qnorm(pvo)
        zo[zo==-Inf] <- 2*min(zo[(zo!=-Inf) | (zo != Inf)])
        zo[zo==Inf] <- 2*max(zo[(zo!=-Inf) | (zo != Inf)])
        
        #======================================================================================#
        #This part calculates the EB estimates of location of the transfomed p-values
        
        z <- rep(0,ngene)
        for (i in 1:ngene) {
                sum1 <- sum(Kh(zo[i]-zo))
                sum2 <- sum((zo[i]-zo)*Kh(zo[i]-zo))
                z[i] <- zo[i] - (h^-2)*(sum2/sum1)                 
        }
        #=====================================================================================#
        #This part computes the step-down p-values for the EB estimates
        
        rk <- order(pvo)
        u <- z[rk]
        ct <-  rep(0,ngene)
        for (ii in 1:ngene) ct[ii] <- sum(u.n[,ii] <= u[ii])  
        p.n <- ct/B
        for (ii in 2:ngene) p.n[ii] <- max(p.n[ii],p.n[ii-1])
        
        #=====================================================================================#
        # Flag the genes based on the step-down p-values
        
        pos <- rk[p.n<alpha]
        
        ##### 10/30/2021: Below is to obtain the EB p-value
        res = p.n[order(rk)]
        
        return(res)
}