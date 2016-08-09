mediation.scan <- function(target, mediator, annotation, qtl.geno, covar = NULL,
                           int,
                           method = c("double-lod-diff", "ignore", "lod-diff", "lod-ratio"), 
                           verbose = TRUE) 
{
  LL <- function(y, X) {
    -length(y)/2 * log10(sum(qr.resid(qr(cbind(X, 1)), y)^2))
  }
  stopifnot(NROW(target) == NROW(mediator))
  stopifnot(NROW(annotation) == NCOL(mediator))
  stopifnot(NROW(qtl.geno) == NROW(target))
  stopifnot(is.null(covar) | NROW(target) == NROW(covar))
  stopifnot(!any(is.na(covar)))
  stopifnot(!any(is.na(qtl.geno)))
  stopifnot(all(is.numeric(target)))
  stopifnot(all(is.numeric(mediator)))
  stopifnot(all(is.numeric(qtl.geno)))
  stopifnot(all(is.numeric(covar)))
  stopifnot(c("CHR", "POS") %in% toupper(names(annotation)))
  method = match.arg(method)
  mediator <- cbind(mediator)
  N <- ncol(mediator)
  if (is.null(covar)) 
    covar <- cbind(rep(1, N))
  LOD <- rep(NA, N)
  qtl.geno.int <- cbind(qtl.geno, qtl.geno*covar[,int])
  if (method == "double-lod-diff") {
    no.na <- !is.na(target)
    LOD0 <- LL(target[no.na], cbind(covar, qtl.geno.int)[no.na, 
                                                         ]) - LL(target[no.na], covar[no.na, ])
  }
  for (i in 1:N) {
    if (verbose & i%%1000 == 0) 
      print(i)
    no.na <- !is.na(target) & !is.na(mediator[, i])
    loglik0 <- LL(target[no.na], cbind(covar[no.na, ], mediator[no.na, 
                                                                i]))
    loglik1 <- LL(target[no.na], cbind(covar[no.na, ], mediator[no.na, 
                                                                i], qtl.geno.int[no.na, ]))
    if (method == "ignore" | (method == "double-lod-diff" & 
                              all(no.na))) {
      LOD[i] <- loglik1 - loglik0
    }
    else {
      loglik2 <- LL(target[no.na], covar[no.na, ])
      loglik3 <- LL(target[no.na], cbind(covar[no.na, ], 
                                         qtl.geno.int[no.na, ]))
      if (method == "lod-diff") {
        LOD[i] <- loglik3 - loglik2 - (loglik1 - loglik0)
      }
      else if (method == "double-lod-diff") {
        LOD[i] <- LOD0 - (loglik3 - loglik2 - (loglik1 - 
                                                 loglik0))
      }
      else if (method == "lod-ratio") {
        LOD[i] <- (10^loglik1 - 10^loglik0)/(10^loglik3 - 
                                               10^loglik2)
      }
    }
  }
  output <- annotation
  output$LOD <- LOD
  class(output) <- c("mediation", "data.frame")
  return(output)
}