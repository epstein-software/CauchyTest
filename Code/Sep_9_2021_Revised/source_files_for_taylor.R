#Function generate gene expression data
gen_express_cov <-function(preds,
                           int,
                           main_beta_g,
                           var_beta_g,
                           Num_off_diag,
                           cov_beta_g,
                           sigma_e,
                           main_beta_z_bin,
                           main_beta_z_con,
                           var_beta_z_bin,
                           var_beta_z_con,
                           cov_beta_z_bin,
                           cov_beta_z_con) {
    g <- preds[1]
    bz <- preds[2] # binary exposure
    cz <- preds[3] # continuous exposure
    
    # create mean vector of expression levels
    mu_vec <- int + main_beta_g*g + main_beta_z_bin*bz + main_beta_z_con*cz  
    
    # create diagonal elements of variance matrix
    Sigma_1 <- diag(sigma_e + var_beta_g*g + var_beta_z_bin*bz+var_beta_z_con*cz) 

    # create covariance elements of variance matrix
    off_diag<-cov_beta_g*g+cov_beta_z_bin*bz+cov_beta_z_con*cz 

upperTriangle(Sigma_1)<-off_diag
Sigma_1[lower.tri(Sigma_1)]<-t(Sigma_1)[lower.tri(Sigma_1)] # make the matrix symmetric

y<-mvrnorm(1,mu_vec,Sigma_1) # generate outcome based on MVN distribution

return(y)

}


#Function to rescale expression data using double GLM
scale_dglm<-function(x,sim_g,sim_z){
    out2 <- dglm::dglm(x~sim_g+sim_z,
                 ~(sim_g+sim_z),
                 family = gaussian)
    scaled <- residuals(out2)/sqrt(predict.glm(out2$dispersion.fit,type='response'))
    return(scaled)
}


#Function to construct scaled marginal association test

SMAT_cor_cov_dglm<-function(sim_y,sim_g,sim_z_bin){

obs_sub_y<-nrow(sim_y)
obs_gene_y<-ncol(sim_y)

scaled_sim_y<-apply(sim_y,2,function(x) scale_dglm(x,sim_g,sim_z_bin))

zsum_temp<-apply(scaled_sim_y,1,function(z)combn(z,2,function(z)z[1]*z[2])) # Create all pairwise combinations of products

Z<-t(zsum_temp)

Xcov<-cbind(rep(1,obs_sub_y),sim_z_bin)

SMAT_pvalue<-SMAT(Z,Xcov,sim_g)$pvalue # Code for this fucntion found in SMAT_code.R

return(SMAT_pvalue)

}



GAMuT_cor_cov_dglm<-function(sim_y,sim_g,sim_z_bin){
    
    obs_sub_y <- nrow(sim_y)
    obs_gene_y <- ncol(sim_y)
    
    scaled_sim_y<-apply(sim_y,2,function(x) scale_dglm(x = x,
                                                       sim_g = sim_g,
                                                       sim_z = sim_z_bin))
    
    # Create all pairwise combinations of products
    zsum_temp<-apply(scaled_sim_y,1,function(z)combn(z,2,function(z)z[1]*z[2])) 
    
    Z<-t(zsum_temp)
    
    Xcov<-cbind(rep(1,obs_sub_y),sim_z_bin)
    
    adjM<-diag(obs_sub_y)-Xcov%*%solve(t(Xcov)%*%Xcov)%*%t(Xcov)
    
    P0<-apply(Z,2,scale,scale=T,center=T)
    P<-as.matrix(P0)
    linear_pheno <- linear_GAMuT_pheno_cov(P,adjM) 
    Yc = linear_pheno$Kc           # Linear kernel similarity matrix
    lambda_Y = linear_pheno$ev_Kc 
    
    project_pheno<-proj_GAMuT_pheno_cov(P,adjM)
    YcP=project_pheno$Kc
    lambda_Yp=project_pheno$ev_Kc
    
    G  = as.matrix(scale(sim_g,center=T,scale=T))
    
    linear_geno <- linear_GAMuT_geno_cov(G,adjM) 
    Xc <- linear_geno$Lc            # Linear kernel similarity matrix
    lambda_X <- linear_geno$ev_Lc   # Eigenvalues of Xc
    
    GAMuT_pvalue_lin= TestGAMuT(Yc,lambda_Y,Xc,lambda_X)
    GAMuT_pvalue_proj= TestGAMuT(YcP,lambda_Yp,Xc,lambda_X)
    
    return(c(GAMuT_pvalue_lin,GAMuT_pvalue_proj))
    
}



## GAMuT-functions.R
##
## 2016-March-13
##
## this script contains functions used in the main program
## for GAMuT analysis
##
##----------------------------------------------------------------------
## descriptions of individual functions:
##----------------------------------------------------------------------
##
## * proj_GAMuT_pheno
##   constructs projection matrix and corresponding eigenvalues for phenotypes
##
## * linear_GAMuT_pheno
##   constructs linear kernel matrix and corresponding eigenvalues for phenotypes
##
## * linear_GAMuT_geno
##   constructs linear kernel matrix and corresponding eigenvalues for genotypes 
##
## * TestGAMuT
##   constructs GAMuT statistic and returns p-value


##----------------------------------------------------------------------
## phenotypic similarity:  projection matrix
##----------------------------------------------------------------------
proj_GAMuT_pheno_cov <- function(X,adjM){

    out = vector("list", 2)
    names(out) = c("Kc", "ev_Kc")

    Kc_temp<-X %*% solve(t(X) %*% X) %*% t(X)

    Kc_adj = adjM%*%Kc_temp%*%adjM
    ev_Kc = eigen(Kc_adj, symmetric=T,only.values=T)$values

    out$Kc = Kc_adj   # projection matrix
    out$ev_Kc = ev_Kc[ev_Kc > 1e-08]
    
#    out$ev_Kc = rep(1, ncol(X))                 # find eigenvalues
    return(out)	
}

##----------------------------------------------------------------------
## phenotypic similarity:  linear kernel
##----------------------------------------------------------------------
linear_GAMuT_pheno_cov <- function(X,adjM){
    out = vector("list", 2)
    names(out) = c("Kc", "ev_Kc")
    
#    Kc_temp = t(X) %*% X   # transposed kernel to find eigenvalues
     Kc_temp = X%*%t(X)
     Kc_adj = adjM%*%Kc_temp%*%adjM
    ev_Kc = eigen(Kc_adj, symmetric=T, only.values=T)$values  
    
    out$Kc = Kc_adj # similiarity kernel for test
    out$ev_Kc = ev_Kc[ev_Kc > 1e-08]
    return(out)
}

##----------------------------------------------------------------------
## genotypic similarity:  linear kernel
##----------------------------------------------------------------------
linear_GAMuT_geno_cov <- function(X,adjM){
    out = vector("list", 2)
    names(out) = c("Lc", "ev_Lc")
    
#    Lc_temp = t(X) %*% X   # transposed kernel to find eigenvalues
     Lc_temp = X %*% t(X)
Lc_adj=adjM%*%Lc_temp%*%adjM
    ev_Lc = eigen(Lc_adj, symmetric=T, only.values=T)$values  
    
    out$Lc = Lc_adj # similiarity kernel for test
    out$ev_Lc = ev_Lc[ev_Lc > 1e-08]
    return(out)	 
}



##----------------------------------------------------------------------
## constructing the GAMuT statistic and deriving p-value:
##----------------------------------------------------------------------
TestGAMuT <- function(Yc, lambda_Y, Xc, lambda_X) {

    ## test statistic:
    m = nrow(Yc) # number of subjects in study
    GAMuT = (1/m) * sum(sum(t(Yc) * Xc))  
    

    ## populate vector of all pairwise combination of eigenvalues
    ## from the phenotype and genotype similarity matrices:
    Z <- (as.matrix(lambda_Y)) %*% t(as.matrix(lambda_X))
    Zsort <- sort(Z, decreasing=T)

    ## derive p-value of GAMuT statistic:
    scoredavies = GAMuT*m^2
    results_score <- davies(scoredavies, Zsort)
    davies_pvalue <- (results_score$Qq)

    return(davies_pvalue)
} 
