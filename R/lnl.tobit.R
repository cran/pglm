lnl.tobit <- function(param, y, X, id, model, link, rn, other = NULL, start.sigma = FALSE, trunc = NULL){

    if (is.null(trunc)){
        lower <- 0
        upper <- + Inf
    }
    else{
        if (length(trunc) == 1){
            lower <- trunc
            upper <- + Inf
        }
        if (length(trunc) == 2){
            lower <- trunc[1]
            upper <- trunc[2]
        }
    }   
        
    if (is.null(other)) other <- "sd"
    mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
    N <- length(y)
    n <- length(unique(id))
    K <- ncol(X)
    beta <- param[1L:K]
    sigma <- param[K+1L]
    if (other == "var") sigma <- sqrt(sigma)
    if (other == "lsd") sigma <- exp(sigma)
    Xb <- as.numeric(crossprod(t(X), beta))
    if (start.sigma){
        ez <- - Xb[y == 0] / sigma
        ep <- (y - Xb)[y > 0] / sigma
        mz <- mills(ez)
        fp <- fs <- numeric(length(y))
        fs[y == 0] <- - (ez + mz) * mz / sigma ^ 2
        fs[y >  0] <- - 1 / sigma ^ 2
        fp[y == 0] <- - mz / sigma
        fp[y >  0] <- ep / sigma
        fp <- tapply(fp, id, sum)
        fs <- tapply(fs, id, sum)
        return(sqrt(2) * sd(- fp / fs))
    }
    YLO <- (y == lower)
    YUT <- (y > lower) & (y < upper)
    YUP <- y == upper
    if (model == "pooling"){
        lnL <- numeric(length = N)
        e <- (y - Xb) / sigma
        mz <- mills(e)
        lnL[YLO] <- pnorm(  e[YLO], log.p = TRUE)
        lnL[YUT] <- dnorm(  e[YUT], log   = TRUE) - log(sigma)
        lnL[YUP] <- pnorm(- e[YUP], log.p = TRUE)
        lnL <-  sum(lnL)
        mz <- mills(e)
        mmz <- mills(- e)
        gradi <- matrix(0, nrow = nrow(X), ncol = ncol(X) + 1L)
        gradi[YLO, 1L:K] <- - mz[YLO] *  X[YLO, , drop = FALSE] / sigma
        gradi[YLO, K+1L] <- -  e[YLO] * mz[YLO]  / (2 * sigma ^ 2)
        gradi[YUT, 1L:K] <-    e[YUT]  * X[YUT, , drop = FALSE] / sigma
        gradi[YUT, K+1L] <- - (1 - e[YUT] ^ 2) / (2 * sigma ^ 2)
        
        gradi[YUP, 1L:K] <-   mmz[YLO] * X[YLO, , drop = FALSE] / sigma
        gradi[YUP, K+1L] <- - mmz[YLO] * e[YLO] / (2 * sigma ^ 2)
        
        if (other == "sd")  gradi[, K+1L] <- gradi[, K+1L] * (2 * sigma)
        if (other == "lsd") gradi[, K+1L] <- gradi[, K+1L] * (2 * sigma ^ 2)
        hbb <- hbs <- hss <- numeric(length = N)
        hbb[YLO] <- - (e[YLO] + mz[YLO]) * mz[YLO] / sigma ^ 2
        hbs[YLO] <- mz[YLO] * (1 - (e[YLO] + mz[YLO]) * e[YLO])/(2 * sigma ^ 3)
        hss[YLO] <- e[YLO] * mz[YLO] * (3 - (e[YLO] + mz[YLO]) * e[YLO]) / (4 * sigma ^ 4)
        hbb[YUT] <- - 1 / sigma ^ 2
        hbs[YUT] <- - e[YUT] / sigma ^ 3
        hss[YUT] <- (1 - 2 * e[YUT] ^ 2) / (2 * sigma ^ 4)
        hbb <- crossprod(hbb * X, X)
        if (other == "sd"){
            hbs <- hbs * (2 * sigma)
            hss <- 4 * sigma ^ 2 * hss + gradi[, K+1L] / sigma
        }
        if (other == "lsd"){
            hbs <- hbs * (2 * sigma ^ 2)
            hss <- 2 * gradi[, K+1L] + 4 * sigma ^ 4 * hss
        }
        hbs <- apply(hbs * X, 2, sum)
        hss <- sum(hss)
        H <-  rbind(cbind(hbb, hbs), c(hbs, hss))
    }
    if (model == "random"){
        smu <- param[K + 2L]
        Pitr <- lapply(rn$nodes,
                       function(z){
                           result <- numeric(length = N)
                           e <- (y - Xb - sqrt(2) * smu * z) / sigma
                           result[YLO] <- pnorm(e[YLO])
                           result[YUT] <- dnorm(e[YUT]) / sigma
                           result
                       }
                       )
        Pir <- lapply(Pitr, function(x) tapply(x, id, prod))
        Li <- Reduce("+", mapply("*", Pir, rn$weights, SIMPLIFY = FALSE)) / sqrt(pi)
        lnL <- sum(log(Li))
        gradi <- Reduce("+",
                        mapply(
                            function(w, x, p){
                                e <- (y - Xb - sqrt(2) * smu * x) / sigma
                                mz <- mills(e)
                                gradi <- matrix(0, nrow = N, ncol = 2)
                                gradi[YLO, 1] <- - mz[YLO] / sigma
                                gradi[YLO, 2] <- - e[YLO] * mz[YLO]  / (2 * sigma ^ 2)
                                gradi[YUT, 1] <- e[YUT]  / sigma
                                gradi[YUT, 2] <- - (1 - e[YUT] ^ 2) / (2 * sigma ^ 2)
                                gradi <- cbind(gradi, gradi[, 1] * sqrt(2) * x)
                                w * as.numeric(p[as.character(id)]) * gradi
                            },
                            rn$weights, rn$nodes, Pir, SIMPLIFY = FALSE
                        )
                        )
        ogradi <- gradi
        if (other == "sd") gradi[, 2L] <- gradi[, 2L] * (2 * sigma)
        if (other == "lsd") gradi[, 2L] <- gradi[, 2L] * (2 * sigma ^ 2)
        gradi <- cbind(gradi[, 1] * X, gradi[, 2:3]) /
            as.numeric(Li[as.character(id)])/ sqrt(pi)
        ogradi <- cbind(ogradi[, 1] * X, ogradi[, 2:3]) /
            as.numeric(Li[as.character(id)])/ sqrt(pi)
        H <- mapply(
            function(w, x, p){
                P <- as.numeric(( p / Li)[as.character(id)])
                sp <- as.numeric(p / Li)
                e <- (y - Xb - sqrt(2) * smu * x) / sigma
                mz <- mills(e)
                gradi <- matrix(0, nrow = N, ncol = 2)
                gradi[YLO, 1] <- - mz[YLO] / sigma
                gradi[YLO, 2] <- - e[YLO] * mz[YLO]  / (2 * sigma ^ 2)
                gradi[YUT, 1] <- e[YUT]  / sigma
                gradi[YUT, 2] <- - (1 - e[YUT] ^ 2) / (2 * sigma ^ 2)
                gradi <- cbind(gradi[, 1] * X, gradi[, 2], gradi[, 1] * sqrt(2) * x)
                gradi <- apply(gradi, 2, tapply, id, sum)
                H1 <- crossprod(sqrt(sp) * gradi)
                hbb <- hbs <- hss <- numeric(length = N)
                hbb[YLO] <- - (e[YLO] + mz[YLO]) * mz[YLO] / sigma ^ 2
                hbs[YLO] <- mz[YLO] * (1 - (e[YLO] + mz[YLO]) * e[YLO])/(2 * sigma ^ 3)
                hss[YLO] <- e[YLO] * mz[YLO] * (3 - (e[YLO] + mz[YLO]) * e[YLO]) / (4 * sigma ^ 4)
                hbb[YUT] <- - 1 / sigma ^ 2
                hbs[YUT] <- - e[YUT] / sigma ^ 3
                hss[YUT] <- (1 - 2 * e[YUT] ^ 2) / (2 * sigma ^ 4)
                hbb <- crossprod(hbb * cbind(X, sqrt(2) * x) * P,
                                 cbind(X, sqrt(2) * x))
                hbs <- apply(hbs * cbind(X, sqrt(2)* x) * P, 2, sum)
                hss <- sum(hss * P)
                H2 <- rbind(cbind(hbb, hbs), c(hbs, hss))
                mX <- H2[1L:K, K+1L]
                sX <- H2[1L:K, K+2L]
                mm <- H2[K+1L, K+1L]
                ss <- H2[K+2L, K+2L]
                H2[K+1L, K+1L] <- ss
                H2[K+2L, K+2L] <- mm
                H2[1L:K, K+1L] <- H2[K+1L, 1L:K] <- sX
                H2[1L:K, K+2L] <- H2[K+2L, 1L:K] <- mX
                (H1 + H2) * w / sqrt(pi)
            },
            rn$weights, rn$nodes, Pir, SIMPLIFY = FALSE
        )
        H <- Reduce("+", H) - crossprod(apply(ogradi, 2, tapply, id, sum))
        if (other == "sd"){
            H[K+1, c(1:K, K+2)] <- H[c(1:K, K+2), K+1] <- H[K+1, c(1:K, K+2)] * (2 * sigma)
            H[K+1, K+1] <- 4 * sigma ^ 2 * H[K+1, K+1] + sum(gradi[, K+1L]) / sigma
        }
        if (other == "lsd"){
            H[K+1, c(1:K)] <- H[c(1:K), K+1] <- H[K+1, c(1:K)] * (2 * sigma ^ 2)
            H[K+2, K+1] <- H[K+1, K+2] <- H[K+2, K+1] * (2 * sigma ^ 2)
            H[K+1, K+1] <- 4 * sigma ^ 4 * H[K+1, K+1] + 2 * sum(gradi[, K+1L])
        }
    }
    attr(lnL, "gradient") <- apply(gradi, 2, sum)
    attr(lnL, "hessian") <- H
    lnL
}

