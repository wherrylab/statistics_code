# Turn a data frame into a model matrix, expanding factors and/or characters with more than two levels into dummy variables
makeFullModelMatrix <- function(df) {
  
  # Class of each column
  colClass <- unlist(lapply(df, class))
  
  # Convert character to factor
  idxCharacter <- which(colClass == "character")
  for(i in idxCharacter) df[,i] <- as.factor(df[,i])
  
  # Recompute colClass
  colClass <- unlist(lapply(df, class))
  
  # Make the model matrix
  idxFactorCol <- colClass == "factor"
  nLevelsPerCol <- unlist(lapply(df, function(el) length(levels(el))))
  options(na.action="na.pass")
  mm <- model.matrix(
    object = ~ .,
    data = df,
    contrasts.arg = lapply(df[idxFactorCol & nLevelsPerCol >= 3], contrasts, contrasts=F)
  )
  mm <- mm[,colnames(mm) != "(Intercept)", drop=F]
  
  # Binary variables get a "1" appended to the column name unnecessarily. The following code identifies these cases and removes the "1" 
  idxBinaryFactors <- rep(NA, ncol(df))
  idxBinaryFactors[!idxFactorCol] <- FALSE
  idxBinaryFactors[nLevelsPerCol != 2] <- FALSE
  idxBinaryFactors[nLevelsPerCol == 2] <- unlist(lapply(df[,nLevelsPerCol == 2, drop=F], function(el) all(levels(el) == c("0", "1"))))
  for(colname in colnames(df[,idxBinaryFactors, drop=F])) {
    colnames(mm)[colnames(mm) == paste0(colname, "1")] <- colname
  }

  # Annotate dummy columns
  idxDummy <- !(colnames(mm) %in% colnames(df))
  
  # Add attributes to model matrix
  mm.attr <- setNames(object = rep(NA, ncol(mm)), nm = colnames(mm))
  mm.attr[colnames(mm) %in% colnames(df)[colClass=="factor"]]  <- "binary"
  mm.attr[colnames(mm) %in% colnames(df)[colClass %in% c("numeric", "integer")]] <- "continuous"
  mm.attr[idxDummy] <- "binary"
  attr(mm, 'varclass') <- mm.attr
  
  # Return final model matrix
  return(mm)
}

# Reproduce some of the behavior of cor.mtest, running non-parametric statistics and making accommodations for nonordinal categorical versus ordered data
cor.mtest.nonparam <- function(df1, df2, var1.ordfactor = NULL, var2.ordfactor = NULL) {
  df1.varclass <- unlist(lapply(df1, class))
  idxCharacter <- which(df1.varclass == "character")
  for(i in idxCharacter) df1[,i] <- as.factor(df1[,i])
  
  df2.varclass <- unlist(lapply(df2, class))
  idxCharacter <- which(df2.varclass == "character")
  for(i in idxCharacter) df2[,i] <- as.factor(df2[,i])
  
  mm1 <- makeFullModelMatrix(df1)
  mm2 <- makeFullModelMatrix(df2)
  nr <- ncol(mm1)
  nc <- ncol(mm2)
  
  cor.mat <- cor(mm1, mm2, use = "pairwise.complete.obs", method = "spearman")
  p.mat <- matrix(NA, nr, nc); rownames(p.mat) <- rownames(cor.mat); colnames(p.mat) <- colnames(cor.mat); 
  n.mat <- matrix(NA, nr, nc); rownames(n.mat) <- rownames(cor.mat); colnames(n.mat) <- colnames(cor.mat); 

  for(i in 1:nr) {
    for(j in 1:nc) {
      print(c(i,j))
      class1 <- attr(mm1, which = "varclass")[i]
      class2 <- attr(mm2, which = "varclass")[j]
      
      if(class1 == "continuous" & class2 == "continuous") {
        x <- mm1[,i]; y <- mm2[,j]
        n.mat[i,j] <- sum(!is.na(x) & !is.na(y))
        isErr <- try(expr = test <- cor.test(x, y, method = "spearman", use = "pairwise.complete.obs"))
        if(class(isErr) == "try-error") {
          p.mat[i,j] <- NA
        } else {
          p.mat[i,j] <- test$p.value
        }
      }
      
      if(class1 == "continuous" & class2 == "binary") {
        y <- mm1[,i]; x <- mm2[,j]
        n.mat[i,j] <- sum(!is.na(x) & !is.na(y))
        isErr <- try(expr = test <- wilcox.test(y ~ x))
        if(class(isErr) == "try-error") {
          p.mat[i,j] <- NA
        } else {
          p.mat[i,j] <- test$p.value
        }
      }
      
      if(class1 == "binary" & class2 == "continuous") {
        y <- mm2[,j]; x <- mm1[,i]
        n.mat[i,j] <- sum(!is.na(x) & !is.na(y))
        isErr <- try(expr = test <- wilcox.test(y ~ x))
        if(class(isErr) == "try-error") {
          p.mat[i,j] <- NA
        } else {
          p.mat[i,j] <- test$p.value
        }
      }
      
      if(class1 == "binary" & class2 == "binary") {
        y <- mm2[,j]; x <- mm1[,i]
        n.mat[i,j] <- sum(!is.na(x) & !is.na(y))
        tab <- table(x,y)
        isErr <- try(expr = test <- fisher.test(tab))
        if(class(isErr) == "try-error") {
          p.mat[i,j] <- NA
        } else {
          p.mat[i,j] <- test$p.value
        }
      }
    }
  }
  
  # add appropriate p-value and FDR behavior for self-correlation and cross-correlation
  fdr.mat <- matrix(NA, nr, nc); colnames(fdr.mat) <- colnames(cor.mat); rownames(fdr.mat) <- rownames(cor.mat); 
  if(nr == nc) {
    if(all(df1 == df2, na.rm = T)) {
      idxUpper <- upper.tri(p.mat)
      fdr.mat[idxUpper] <- p.adjust(p.mat[idxUpper], method = "BH")
      idxLower <- lower.tri(p.mat)
      fdr.mat[idxLower] <- p.adjust(p.mat[idxLower], method = "BH")
      diag(p.mat) <- NA
    }
  } else {
    fdr.mat <- matrix(p.adjust(p.mat, method = "BH"), nr, nc)
    colnames(fdr.mat) <- colnames(cor.mat); rownames(fdr.mat) <- rownames(cor.mat)
  }
  
  ret <- list(p = p.mat, fdr = fdr.mat, n = n.mat, scor = cor.mat)
  return(ret)
}

# Capitalize the first character in a string, can handle vectors
CapStr <- function(str) {
  splitList <- strsplit(str, " ")
  splitList <- lapply(splitList, function(el) {
    paste(toupper(substring(el, 1,1)), substring(el, 2),
          sep="", collapse=" ")
    
  })
  return(unlist(splitList))
}

# Calculate all unique permutations of a binary vector (note that this can take a very long time for n >= 25)
allBinaryVecPermutation <- function(n,m) {
  t(combn(n,m,tabulate,nbins=n))
}

# Calculate significance of difference in Spearman correlation coefficient by exact permutation test (note that this can take a very long time for n >= 25)
corrDiffPermTestExact <- function(var1, var2, label, matPerm=NULL) {
  if(is.null(matPerm)) {
    matPerm <- allBinaryVecPermutation(length(label), min(table(label)))
  }
  
  corrOrig <- tapply(1:length(label), label, function(idx) cor.test(var1[idx], var2[idx], method = "spearman")) 
  corrOrig <- corrOrig[[2]]$estimate - corrOrig[[1]]$estimate
  
  refCorr <- rep(NA, nrow(matPerm))
  for(i in 1:nrow(matPerm)) {
    if(i %% 1000 == 0) print(sprintf("Completed %d of %d permutations", i, nrow(matPerm)))
    corrPerm <- tapply(1:length(label), matPerm[i,], function(idx) cor.test(var1[idx], var2[idx], method = "spearman"))
    refCorr[i] <- corrPerm[[2]]$estimate - corrPerm[[1]]$estimate
  }
  
  
  out <- list(
    ref.dist = refCorr,
    corr.diff = corrOrig,
    p.val = min(mean(refCorr <= corrOrig) * 2, mean(refCorr >= corrOrig) * 2, 1)
  )
  
  return(out)
}

# Calculate significance of difference in Spearman correlation coefficient by simmulated permutation test (preferable for large "n")
corrDiffPermTestSim <- function(var1, var2, label, nSim=1e4, seed=2938457) {
  set.seed(seed)
  corrOrig <- tapply(1:length(label), label, function(idx) cor.test(var1[idx], var2[idx], method = "spearman")) 
  corrOrig <- corrOrig[[2]]$estimate - corrOrig[[1]]$estimate
  
  refCorr <- rep(NA, nSim)
  for(i in 1:nSim) {
    if(i %% 1000 == 0) print(sprintf("Completed %d of %d permutations", i, nSim))
    corrPerm <- tapply(1:length(label), sample(label), function(idx) cor.test(var1[idx], var2[idx], method = "spearman"))
    refCorr[i] <- corrPerm[[2]]$estimate - corrPerm[[1]]$estimate
  }
  
  
  out <- list(
    ref.dist = refCorr,
    corr.diff = corrOrig,
    p.val = min(ifelse(sign(corrOrig) == 1, mean(refCorr >= corrOrig) * 2, mean(refCorr <= corrOrig) * 2), 1)
  )
  
  return(out)
}



