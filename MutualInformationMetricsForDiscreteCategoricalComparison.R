# Functions for calculating entropy and mutual information from vectors of discrete categorical labels

# Calculate mutual information between two sets of discrete categorical labels
calcMutualInformationTwoLabelSets <- function(label1, label2, base=exp(1)) {
  class1Count <- table(label1)
  class2Count <- table(label2)
  jointCount <- table(label1, label2)
  class1Names <- names(class1Count)
  class2Names <- names(class2Count)
  
  px   <- class1Count / sum(class1Count)
  py   <- class2Count / sum(class2Count)
  pxy  <-  jointCount / sum( jointCount)
  
  sumTerms <- matrix(data = NA, nrow = nrow(pxy), ncol = ncol(pxy))
  for(i in 1:length(class1Names)) {
    for(j in 1:length(class2Names)) {
      sumTerms[i,j] <- ifelse(
        pxy[class1Names[i],class2Names[j]] == 0, 
        0, 
        pxy[class1Names[i],class2Names[j]] * log(pxy[class1Names[i],class2Names[j]] / px[class1Names[i]] / py[class2Names[j]]) / log(base)
      )
    }
  }
  out <- sum(sumTerms)
  
  return(out)
}

# Calculate self entropy for one set of discrete categorical labels
calcSelfEntropyOneLabelSet <- function(label1, base=exp(1)) {
  class1Count <- table(label1)
  class1Names <- names(class1Count)
  
  px   <- class1Count / sum(class1Count)
  
  sumTerms <- rep(NA, length(class1Count))
  for(i in 1:length(class1Names)) {
    sumTerms[i] <- ifelse(
      px[class1Names[i]] == 0,
      0,
      -px[class1Names[i]] * log(px[class1Names[i]]) / log(base)
    )
  }
  out <- sum(sumTerms)
  
  return(out)
}

# Calculate conditional entropy between two sets of discrete categorical labels
calcConditionalEntropyTwoLabelSets <- function(label1, label2, base=exp(1)) {
  class1Count <- table(label1)
  class2Count <- table(label2)
  jointCount <- table(label1, label2)
  class1Names <- names(class1Count)
  class2Names <- names(class2Count)
  
  px   <- class1Count / sum(class1Count)
  py   <- class2Count / sum(class2Count)
  pxy  <-  jointCount / sum( jointCount)
  
  sumTerms <- matrix(data = NA, nrow = nrow(pxy), ncol = ncol(pxy))
  for(i in 1:length(class1Names)) {
    for(j in 1:length(class2Names)) {
      sumTerms[i,j] <- ifelse(
        pxy[class1Names[i],class2Names[j]] == 0, 
        0, 
        -pxy[class1Names[i],class2Names[j]] * log(pxy[class1Names[i],class2Names[j]] / px[class1Names[i]]) / log(base)
      )
    }
  }
  out <- sum(sumTerms)
  
  return(out)
}

# Calculate the information quality ratio between two discrete categorical sets of labels
# Ref: https://en.wikipedia.org/wiki/Mutual_information
# Ref: https://www.sciencedirect.com/science/article/abs/pii/S0169743916304907
calcIQR <- function(label1, label2, base=exp(1)) {
  IXY  <- calcMutualInformationTwoLabelSets(label1, label2, base=base)
  HXY  <- calcSelfEntropyOneLabelSet(paste(label1, label2), base=base)
  IQR <- IXY/HXY
  
  return(IQR)
}

# Simulate permutation distribution of IQR between vectors x and y
iqrPermutationTest <- function(x, y, nSim = 1e3) {
  simIQR <- rep(NA, nSim)
  for(i in 1:nSim) {
    idxPermutation <- sample(length(y), replace = FALSE)
    y.perm <- y[idxPermutation]
    simIQR[i] <- calcIQR(x, y.perm)
  }
  return(simIQR)
}

# Simulate permutation distribution of IQR difference between a comparitor and reference
iqrPermutationDiffTest <- function(ref.clust, ref.class, comp.clust, comp.class, nSim = 1e3) {
  simIQRDiff <- rep(NA, nSim)
  for(i in 1:nSim) {
    idxPermutation  <- sample(length(comp.class) + length(ref.class), replace = FALSE)
    class.perm      <- c(ref.class, comp.class)[idxPermutation]
    ref.class.perm  <- class.perm[   1:length(ref.class) ]
    comp.class.perm <- class.perm[- (1:length(ref.class))]
    simIQRDiff[i]   <- calcIQR(comp.clust, comp.class.perm) - calcIQR(ref.clust, ref.class.perm)
  }
  return(simIQRDiff)
}

# Simulate permutation distribution of IQR difference between a comparitor and reference; condition on class (cell type) fractions remaining constant between comparitor and reference in case this biases analysis
iqrPermutationDiffCondTest <- function(ref.clust, ref.class, comp.clust, comp.class, nSim = 1e3) {
  simIQRDiff <- rep(NA, nSim)
  for(i in 1:nSim) {
    idxPermutation  <- sample(length(ref.class), replace = FALSE)
    ref.class.perm  <- ref.class[idxPermutation]
    idxPermutation  <- sample(length(comp.class), replace = FALSE)
    comp.class.perm <- comp.class[idxPermutation]
    simIQRDiff[i]   <- calcIQR(comp.clust, comp.class.perm) - calcIQR(ref.clust, ref.class.perm)
  }
  return(simIQRDiff)
}

# Calculate two-sided p-values from a permutation test simulation according to Kulinskaya 2008 (arXiv:0810.2124)
calcPvalueFromPermutationTest <- function(actual.stat, perm.sim, two.sided.center=NULL, two.sided=TRUE, lower.tail=FALSE) {
  if(two.sided == FALSE) {
    if(lower.tail == TRUE) {
      p <- mean(perm.sim <= actual.stat)
    }
    if(lower.tail == FALSE) {
      p <- mean(perm.sim >= actual.stat)
    }
    if(p == 0) {
      p <- sprintf("<%.2g", 1/length(perm.sim))
    }
    return(p)
  }
  
  if(two.sided == TRUE) {
    if(is.null(two.sided.center)) {
      two.sided.center <- quantile(perm.sim, 0.5)
    }
    p <- min(
      mean(perm.sim <= actual.stat) / mean(perm.sim <= two.sided.center), 
      mean(perm.sim >= actual.stat) / mean(perm.sim >= two.sided.center), 
      1
    )
    if(p == 0) {
      p <- sprintf("<%.2g", 2/length(perm.sim))
    }
    return(p)
  }
}

###### Demo code starts here


label1 <- c(rep('a', 10), rep('b', 8), rep('c', 6))
label2 <- label1
label2[label2 == 'a'] <- 'd'
label2[label2 == 'b'] <- 'e'
label2[label2 == 'c'] <- 'f'
#label2[1:2] <- "g"
label1[1:2] <- "g"

base <- 2 # Base of information metric, "2" corresponds to "bits"

# Two label set example, label 1 exactly predicts label 2
table(label1)
table(label2)
table(label1,label2)
HX   <- calcSelfEntropyOneLabelSet(label1, base=base)
HY   <- calcSelfEntropyOneLabelSet(label2, base=base)
HYgX <- calcConditionalEntropyTwoLabelSets(label1, label2, base=base)
HXgY <- calcConditionalEntropyTwoLabelSets(label2, label1, base=base)
IXY  <- calcMutualInformationTwoLabelSets(label1, label2, base=base)
HXY  <- calcSelfEntropyOneLabelSet(paste(label1, label2), base=base)
IQR <- IXY/HXY
CXY <- IXY / HY
CYX <- IXY / HX
setNames(c(HX, HY, HYgX, HXgY, IXY, HXY, IQR, CXY, CYX), c("HX", "HY", "HYgX", "HXgY", "IXY", "HXY", "IQR", "CXY", "CYX"))

# Two label set example, label 1 and label 2 correlation broken by permutation
label2Permuted <- sample(label2)
table(label2Permuted)
table(label1)
table(label1,label2Permuted)
HX   <- calcSelfEntropyOneLabelSet(label1, base=base)
HY   <- calcSelfEntropyOneLabelSet(label2Permuted, base=base)
HYgX <- calcConditionalEntropyTwoLabelSets(label1, label2Permuted, base=base)
HXgY <- calcConditionalEntropyTwoLabelSets(label2Permuted, label1, base=base)
IXY  <- calcMutualInformationTwoLabelSets(label1, label2Permuted, base=base)
HXY  <- calcSelfEntropyOneLabelSet(paste(label1, label2Permuted), base=base)
IQR <- IXY/HXY
setNames(c(HX, HY, HYgX, HXgY, IXY, HXY, IQR), c("HX", "HY", "HYgX", "HXgY", "IXY", "HXY", "IQR"))

# Two label set example, label 1 and label 2 correlation broken by permutation and label 2 with one additional category "g"
label2Permuted_4class <- label2Permuted
label2Permuted_4class[5] <- "g"
table(label2Permuted_4class)
table(label1)
table(label1,label2Permuted_4class)
HX   <- calcSelfEntropyOneLabelSet(label1, base=base)
HY   <- calcSelfEntropyOneLabelSet(label2Permuted_4class, base=base)
HYgX <- calcConditionalEntropyTwoLabelSets(label1, label2Permuted_4class, base=base)
HXgY <- calcConditionalEntropyTwoLabelSets(label2Permuted_4class, label1, base=base)
IXY  <- calcMutualInformationTwoLabelSets(label1, label2Permuted_4class, base=base)
HXY  <- calcSelfEntropyOneLabelSet(paste(label1, label2Permuted_4class), base=base)
IQR <- IXY/HXY
setNames(c(HX, HY, HYgX, HXgY, IXY, HXY, IQR), c("HX", "HY", "HYgX", "HXgY", "IXY", "HXY", "IQR"))

# Two label set example, label 1 and label 2 correlation preserved, but label 2 with one additional category "g"
label2_4class <- label2
label2_4class[5] <- "g"
table(label2_4class)
table(label2)
table(label1,label2_4class)
HX   <- calcSelfEntropyOneLabelSet(label1, base=base)
HY   <- calcSelfEntropyOneLabelSet(label2_4class, base=base)
HYgX <- calcConditionalEntropyTwoLabelSets(label1, label2_4class, base=base)
HXgY <- calcConditionalEntropyTwoLabelSets(label2_4class, label1, base=base)
IXY  <- calcMutualInformationTwoLabelSets(label1, label2_4class, base=base)
HXY  <- calcSelfEntropyOneLabelSet(paste(label1, label2_4class), base=base)
IQR <- IXY/HXY
setNames(c(HX, HY, HYgX, HXgY, IXY, HXY, IQR), c("HX", "HY", "HYgX", "HXgY", "IXY", "HXY", "IQR"))

# Example permutation test 
set.seed(3948752)
permutationDistribution <- iqrPermutationTest(label1, label2Permuted_4class)
calcPvalueFromPermutationTest(calcIQR(label1, label2Permuted_4class), permutationDistribution)

