# Shared Statistics Functions in R

## Purpose

This repository allows for sharing common R statistics functions for use within the Wherry lab and with external users. This includes code that was used to compute non-parametric correlation heatmaps as in Mathew *et al.* Science, 2020.

The cor.mtest.nonparam() function included here is meant to serve as a non-parametric substitute for the cor.mtest() function of the *corrplot* library that substitutes appropriate non-parametric statistical association tests depending on the underlying data types (e.g. discrete unordered, discrete ordered, integer, or continuous scale). The end user is responsible to make sure that the data are appropriately encoded prior to using this function.

## Technical Notes

The non-parametric correlation code will compute a correlation matrix of Spearman rank correlation coefficients along with associated P-values and FDR-values based on all pairwise column combinations between two input data frames. Assuming that the data are appropriately coded (i.e. ordered data including continuous, integer, or ordinal factor data should be encoded as "numeric" and unordered discrete data should be encoded as "factor"), the following statistical tests will be performed by default:

1. For pairwise ordered vs ordered associations: Spearman rank correlation test

2. For pairwise unordered discrete vs ordered associations: Unpaired Wilcoxon test

3. For pairwise unordered discrete vs unordered discrete associations: Fisher's exact test

Note that "factor" or "character" variables will automatically be expanded to binary "dummy" variables for the purpose of statistical testing. If this is not desired behavior, this code may not be appropriate or may need to be modified. FDR-correction is performed on the resulting matrix of P-values (on-diagonal and redundant tests are removed in the case of self-correlation) using the method of Benjamini-Hochberg.

## Usage

The relevant functions can be accessed by downloading and running source() on the R code file:

source("./NonParametricCorrelationCode.R")

To calculate self-correlation between all columns of a data frame "df" with itself:

out <- cor.mtest.nonparam(df1 = df, df2 = df)

To calculate non-self-correlation between all columns of a data frame "df1" with another data frame "df2":

out <- cor.mtest.nonparam(df1 = df1, df2 = df2)

The output is a list with matrix elements "scor", "p", "fdr", and "n" referring to matrices of pairwise Spearman rank correlation coefficients, P-values, FDR-values, and number of pairwise complete observations, respectively.

## Disclaimer

The code is provided as is. Use and modify at your own risk.

## Contact

Questions about the code can be directed to the principle code author: Derek A. Oldridge

Email: derek.oldridge@pennmedicine.upenn.edu
