# BCV2 Analysis

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(GGally)
library(dplyr) 
library(tidyverse)
library(dendextend)
library(MASS)
library(lubridate)
library(Matrix)
library(ggpubr)
library(factoextra)
library(ggfortify)
library(ISLR)


# reading in data
BCV2 = read.csv("BCV2.csv")
      # Change Label column to factor
BCV2 = BCV2 %>%
  mutate( Label = factor( Label ))

# subset of data
BCV2_2 = BCV2[,-1]

### Q1
## (a)
# answer from r environment: 375 obs of 28 variables

## (b)
# parallel coordinate plot - raw
ggparcoord(BCV2_2, scale = "globalminmax")
# from plot - variable with largest range: 13 or X13

# parallel coordinate plot - scaled
BCV2_2_scaled = scale(BCV2_2, center = TRUE, scale = TRUE)
ggparcoord(BCV2_2_scaled, scale = "globalminmax")

## (c)
# (i)
# covariance matrix
cov_mat = cov(BCV2_2)

# eigenvalues of covariance matrix with plot
eigenvals = eigen(cov_mat)$value
plot(eigenvals, type="b")

# (ii)
# condition number of covariance matrix
kappa(cov_mat, exact = TRUE)

## (d)
# scaled covariance matrix
cov_mat_scaled = cov(BCV2_2_scaled)

# trace
trace = sum(diag(cov_mat_scaled))
trace
# trace is 28

# entry S12
cov_mat_scaled[1,2]
# S12 is 0.343472

## (e)
# (i)
# PCA
scaled.pca = prcomp(BCV2_2_scaled)

# score plots
# PC1 against PC2
plot(scaled.pca$x[,1], scaled.pca$x[,2], xlab='PC1', ylab="PC2")
# PC1 against PC3
plot(scaled.pca$x[,1], scaled.pca$x[,3], xlab='PC1', ylab="PC3")
# PC2 against PC3
plot(scaled.pca$x[,2], scaled.pca$x[,3], xlab='PC2', ylab="PC3")

# BETTER WAY
autoplot(scaled.pca, x=1, y=2, data = scaled.pca, main = "PC1 against PC2 - scaled data")
autoplot(scaled.pca, x=1, y=3, data = scaled.pca, main = "PC1 against PC3 - scaled data")
autoplot(scaled.pca, x=2, y=3, data = scaled.pca, main = "PC2 against PC3 - scaled data")

# (ii)
# density plots
plot(density(scaled.pca$x[,1]), xlab = "PC1")
plot(density(scaled.pca$x[,2]), xlab = "PC2")
plot(density(scaled.pca$x[,3]), xlab = "PC3")

## (f)
# (i)
# biplot of first two eigenvectors (note: col part shows arrows only in blue)
biplot(scaled.pca, choices = 1:2, col = c("white", "blue") )

### Q2
## (a)
# (i)
# dimension from r env is 28 (BCV2_2)

# (ii) 
# biplot for first two factor loadings
eig_mat_scaled = eigen(cov_mat_scaled)
sqrt_lambda = diag(eig_mat_scaled$values^(0.5))
      # use dimension as upper bound for first part
Factor.two = eig_mat_scaled$vectors[1:28,1:2]  %*%  sqrt_lambda[1:2,1:2]
biplot(Factor.two, Factor.two)

## (b)
# (i)
# no rotation - factor analysis, loadings and biplot
fa.none = factanal(BCV2_2_scaled, 2, rotation = "none", scores = "regression" ) # this is the k=2 case, plain option (likelihood?)
fa.none$loadings
biplot( fa.none$scores, fa.none$loadings, col = c("white", "blue") )

# (ii)
# varimax rotation - factor analysis, loadings and biplot
fa.varimax = factanal(BCV2_2_scaled, 2, rotation = "varimax", scores = "regression" )
fa.varimax$loadings
biplot( fa.varimax$scores, fa.varimax$loadings, col = c("white", "dark green") )

## (c)
# (i)
# hypothesis test of basic factor analysis case (k=2)
fa.none
# The chi square statistic is 10973.49 on 323 degrees of freedom.
# The p-value is 0

# (ii)
# hypothesis test for k=3 => only works when increase lower bound to 0.01 (default is 0.005)
fa.k3 = factanal(BCV2_2_scaled, 3, rotation = "none", scores = "regression", lower = 0.01)
fa.k3

### Q3
## (a)
# (i)
# dendrogram - complete linkage and Euclidean distance
dend = as.dendrogram( hclust( dist((BCV2_2) ), 
                                      method = "complete" ) 
)
plot(dend, main = "Dendrogram", label = FALSE, ylab = "Height")
# Table from dendrogram
dend_classes <- cutree(dend, k = 1:12) # cuts the tree at various levels to get maximum of k clusters

dend_classes.fp = as.data.frame( dend_classes ) %>%
  pivot_longer(cols = 1:12,
               names_to = "max_clust",
               values_to = "cluster"
  )
# cluster table
clustable = xtabs( ~ cluster + max_clust, dend_classes.fp )
clustable

# (ii)
# dendrogram - complete linkage and maximum distance
dend2 = as.dendrogram( hclust( dist((BCV2_2), method = "maximum" ), 
                              method = "complete" )
)
plot(dend2, main = "Dendrogram 2", label = FALSE, ylab = "Height")
# Table from dendrogram
dend_classes2 <- cutree(dend2, k = 1:12) # cuts the tree at various levels to get maximum of k clusters

dend_classes.fp2 = as.data.frame( dend_classes2 ) %>%
  pivot_longer(cols = 1:12,
               names_to = "max_clust",
               values_to = "cluster"
  )
# cluster table
clustable2 = xtabs( ~ cluster + max_clust, dend_classes.fp2 )
clustable2

# (iii)
# dendrogram - single linkage and Euclidean distance
dend3 = as.dendrogram( hclust( dist((BCV2_2) ), 
                              method = "single" ) 
)
plot(dend3, main = "Dendrogram 3", label = FALSE, ylab = "Height")
# Table from dendrogram
dend_classes3 <- cutree(dend3, k = 1:12) # cuts the tree at various levels to get maximum of k clusters

dend_classes.fp3 = as.data.frame( dend_classes3 ) %>%
  pivot_longer(cols = 1:12,
               names_to = "max_clust",
               values_to = "cluster"
  )
# cluster table
clustable3 = xtabs( ~ cluster + max_clust, dend_classes.fp3 )
clustable3

# (v)


# single linkage and Euclidean distance table
single_linkage = xtabs( ~ dend_classes3[,7] + BCV2$Label)
single_linkage

# complete linkage and maximum distance table
max_linkage = xtabs( ~ dend_classes2[,7] + BCV2$Label)
max_linkage

## (b)
# (i)
# likelihood-based LDA with estimated prior probabilities
lda.0 = lda( Label ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25+X26+X27+X28, data = BCV2 )
      # confusion matrix 
preds.0 = predict( lda.0 )$class
table(preds.0, BCV2$Label)

# (ii)
# with uniform priors - depends on number of factors (2 factors here)
# (scaling doesn't matter)

# likelihood-based LDA with uniform priors
lda.1 = lda( Label ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25+X26+X27+X28, data = BCV2, prior = c(1,1)/2)
      # confusion matrix
preds.1 = predict( lda.1 )$class
table(preds.1, BCV2$Label)

## (c)
# (i)
# qda with estimated prior probabilities
qda.0 = qda( Label ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25+X26+X27+X28, data = BCV2 )
      # confusion matrix 
preds.2 = predict( qda.0 )$class
table(preds.2, BCV2$Label)

# (ii)
# qda with uniform priors
qda.1 = qda( Label ~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+X25+X26+X27+X28, data = BCV2, prior = c(1,1)/2)
      # confusion matrix 
preds.3 = predict( qda.1 )$class
table(preds.3, BCV2$Label)
mean(preds.3 != BCV2$Label)
