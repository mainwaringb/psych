# Issue with Psych Polychoric function
Ben Mainwaring

28/02/2020

Submitted to William Revelle (psych package maintainer) on 28 Feb 2020

# Introduction

`Psych` is one of my most-used R package. However, there are two related issues with the package's `polychoric` function. While these issues should not affect *most* calculations of polychoric correlations, they may affect a fair number of correlations in batteries with uneven numbers of response items (for example, some 4-pt scales and some 6-pt scales). In addition, they might quasi-randomly affect outcomes or computational times for a small number of other analyses.

The problem was identified using psych 1.19.12.31, with R v.3.6.0 and Mac OS 10.14.6. However, the bug should likely replicate across R versions and operating systems (the bug and suggested fix have been reported, so the issue should be addressed in future releases of `psych`). 

# Illustration of Issues

## Issue A
The polychoric function is designed to coerce `global = FALSE` when it detects an uneven number of response alternatives. At present, the code to do so is broken. For example:

```{r}
library(psych)
set.seed(6171729)
mydf <- data.frame(x1 = sample(1:6, 1000, replace = TRUE), x2 = sample(1:6, 1000, replace = TRUE), x3 = sample(1:4, 1000, replace = TRUE), x4 = sample(1:4, 1000, replace = TRUE))
polychoric(mydf, global = TRUE)
```

No error message occurs (although, in the code above, an unrelated warning about matpLower does occur, related to the continuity correction).
 
In the compiled source (accessed via `print(polychoric)`) the relevant code is lines 123-129:
 
    xmax <- apply(x, 2, function(x) max(x, na.rm = TRUE))
    xmax <- max(xmax)
    gmaxx <- gmaxy <- xmax
    if (min(xmax) != max(xmax)) {
      global <- FALSE
      warning("The items do not have an equal number of response alternatives, global set to FALSE.")
    }
 
However, the `if` condition here is broken. The second line of this code (line 124 in the source) converts `xmax` from a vector (of column maximums) to a scalar (the global maximum). `min(any_scalar) != max(any_scalar)` will always evaluate to `FALSE`, so the if conditions will never run. Simply reordering the lines of code, so the `if` statement goes above `xmax <- max(xmax)` should resolve the issue 
 
## Issue B
The polychoric function with `global = FALSE` and `correct` not set to `FALSE` sometimes return an error like the following (or slightly different if `multicore` is turned off). Most commonly the error occurs if variables have different numbers of levels.

```{r, error = TRUE}
polychoric(mydf, global = FALSE, correct = TRUE)
```

The issue seems to lie in the `polyc` function (an internal function used in the computation of `polychoric`), specifically the following section (spacing and comments edited for clarity):

     	 	zerorows <- apply(tab, 1, function(x) all(x == 0))
    		zerocols <- apply(tab, 2, function(x) all(x == 0))
    		zr <- sum(zerorows)
    		zc <- sum(zerocols)
    		tab <- tab[!zerorows, ,drop=FALSE]  
    		tab <- tab[, !zerocols, drop=FALSE] 
        csum <- colSums(tab)
        rsum <- rowSums(tab)
        	#if(correct > 0) tab[tab==0] <- correct/tot
     		if(min(dim(tab)) < 2) {rho <- list(objective = NA) } else {
      	 	 cc <-  qnorm(cumsum(csum)[-length(csum)])
       		 rc <-  qnorm(cumsum(rsum)[-length(rsum)])
      	 	 rho <- optimize(polyF,interval=c(-1,1),rc=rc, cc=cc,tab)
      	 	
          	}

        cc <- qnorm(cumsum(csum)[-length(csum)])
        rc <- qnorm(cumsum(rsum)[-length(rsum)])
        rho <- optimize(polyF, interval = c(-1, 1), 
          rc = rc, cc = cc, tab)

Essentially, the code above generates a two-way crosstabulation of data, gets the cumulative row and column percentiles, and finds the z-score associated with the percentile. 
The `qnorm` returns z-scores from percentiles. The function needs a value between 0 and 1, and will return `NA` if it receives larger or smaller values. In most cases, `cumsum(rsum)[-length(rsum)]` should be bounded by 0 and 1, as it's constructed from a two-way crosstabulation (`tab`) which sums to 1. 
 
However, the continuity correction changes the sum of `tab` to > 1.This can introduce values greater than 1.0 into `cumsum(rsum)[-length(rsum)]`. Values greater than 1.0 generate an `NA` in `qnorm`, and causes an error when `polyF` is called with an `NA` value. The continuity correction can also introduce an value equal to 1.0, which will evaluate to `qnorm(1.0) = Inf` getting passed to `polyF`.

# Proposed Solutions

## Issue A

Simply reordering the lines of code, so the `if` statement goes above `xmax <- max(xmax)` should resolve the issue:

    xmin <- apply(x,2,function(x) min(x,na.rm=TRUE))
    #if(global)  { xmin <- min(xmin)} 
    xmin <- min(xmin)
    x <- t(t(x) - xmin +1)  #all numbers now go from 1 to nvalues
    
    gminx <- gminy <- 1  #allow for different minima if minmax is null
    xmax <- apply(x,2,function(x)  max(x,na.rm=TRUE))
    #if(global) xmax <- max(xmax)     
    if (min(xmax) != max(xmax)) {
      global <- FALSE
      warning("The items do not have an equal number of response alternatives, global set to FALSE.")
    }
    xmax <- max(xmax)  #don't test for globality xmax

    gmaxx <- gmaxy <- xmax #check for different maxima

## Issue B

There are probably two components to the solution. The first is to move the following lines of `polyc` earlier in the script, before introducing the continuity correction:
 
      zerorows <- apply(tab, 1, function(x) all(x == 0))
      zerocols <- apply(tab, 2, function(x) all(x == 0))
 
In other words, we want to identify zero rows so that they can be dropped, rather than add a continuity correction to zero rows so they are no longer flagged.

The second component is rescaling the table (`tab`) so it sums to 1.0 after adding a continuity correction: `tab <- tab/sum(tab)`. I think this is in line with how continuity corrections are supposed to work, although my theoretical background on polychoric/latent correlations is limited. Together, these changes should produce the following code:
 
     
    if(!is.na(sum(tab)) ){
      zerorows <- apply(tab, 1, function(x) all(x == 0))
      zerocols <- apply(tab, 2, function(x) all(x == 0))
      zr <- sum(zerorows)
      zc <- sum(zerocols)
      tab <- tab[!zerorows, ,drop=FALSE]  
      tab <- tab[, !zerocols, drop=FALSE] 
    }
    
     if(correct > 0) {if(any(tab[]==0)) {
       fixed <- 1
       tab[tab==0] <- correct/tot 
       tab <- tab / sum(tab, na.rm = TRUE) #not sure if this will work with NAs
       }
    } #moved from below 16.6.22
  
There are other possible ways to handle rescaling the data table. It would be possible (although probably not advised) to leave the table summing to 1.0
