piToTime <- function(timing_param, type=c("Mono-allelic Gain","CN-LOH", "Bi-allelic Gain (WGD)")){
    type <- match.arg(type)
    evo <- NA
    w <- timing_param[,"s"]==1 &! timing_param[,"mixFlag"] ## Clonal states, not mixture (CN subclonal)
    M <- sum(w) ## Max multiplicity = Major copy number
    m <- sum(w & timing_param[,"n.m.s"]==2) ## Minor CN
    #evo <- paste0("1:1", paste0(2,":",m), if(M>2) paste0(3,":",) else "", collapse="->")
    t <- timing_param[(M:(M-1)),c("T.m.sX","T.m.sX.lo","T.m.sX.up"), drop=FALSE] ## Timing M and M-1
    if(nrow(t)==1) t <- rbind(t, NA)
    if(!any(is.na(t))){
        if(M==4) {
            if(timing_param[M-1,"P.m.sX"] < 0.02){ ## No MCN 3
                if(type=="CN-LOH"){ ## Hotfix for 4+0, treated at 1:1 -> 2:0 -> 4:0
                    t <- timing_param[c(M,M-2),c("T.m.sX","T.m.sX.lo","T.m.sX.up"), drop=FALSE]*c(1,0.5) ## Timing M and M-2
                    evo <- "1:1->2:0->4:0"
                }
                else if(type=="Bi-allelic Gain (WGD)"){         
                    if(m==2) {## Hotfix for 4+2 regions, treated at 1:1 -> 2:1 -> 4:2
                        t <- timing_param[c(M,M-2),c("T.m.sX","T.m.sX.lo","T.m.sX.up"), drop=FALSE] ## Timing M and M-2
                        t[2,] <- pmax(0,2*t[2,] - t[1,])/3
                        evo <- "1:1->2:1->4:2"
                    } else if(m==4) {## Hotfix for 4+4 regions, treated at 1:1 -> 2:2 -> 4:4
                        t <- timing_param[c(M,M-2),c("T.m.sX","T.m.sX.lo","T.m.sX.up"), drop=FALSE]*c(1,0.5) ## Timing M and M-2
                        evo <- "1:1->2:2->4:4"
                    }
                }           
            } else if(type=="Bi-allelic Gain (WGD)"){ ## Can't uniquely time second event
                t[2,] <- NA
            } 
            if(m==3) {
                t[2,] <- NA ## Don'time secondary 4+3 for now, needs more work
            } 
        } else {
            if(M==3 & type=="Bi-allelic Gain (WGD)") {## Hotfix for 3+2 regions, treated at 1:1 -> 2:2 -> 3:2
                t[2,] <- pmax(0,2*t[2,] - t[1,])
                evo <- "1:1->2:2->3:2"
            }
        }
    }
    colnames(t) <- c("","lo","up")
    t[,2] <- pmin(t[,1],t[,2])
    t[,3] <- pmax(t[,1],t[,3])
    t <- pmin(apply(t,2,cumsum),1) ## Times are actually deltas
    if(M < 3) t[2,] <- NA
    t[is.infinite(t)] <- NA
    rownames(t) <- c("",".2nd")
    return(c(t[1,],`.2nd`=t[2,])) ## Linearise
}

bbToTime <- function(bb, timing_param = bb$timing_param, pseudo.count=5){
    sub <- duplicated(bb) 
    covrg <- countQueryHits(findOverlaps(bb, bb)) 
    maj <- sapply(timing_param, function(x) if(length(x) > 0) x[1, "majCNanc"] else NA) #bb$major_cn
    min <- sapply(timing_param, function(x) if(length(x) > 0) x[1, "minCNanc"] else NA) #bb$minor_cn
    type <- sapply(seq_along(bb), function(i){
                if(maj[i] < 2 | is.na(maj[i]) | sub[i] | (maj[i] > 4 & min[i] >= 2)) return(NA)
                type <- if(min[i]==1){ "Mono-allelic Gain" 
                        }else if(min[i]==0){"CN-LOH"}
                        else "Bi-allelic Gain (WGD)"
                return(type)
            })
    time <- t(sapply(seq_along(bb), function(i){
                        if(sub[i] | is.na(type[i])) return(rep(NA,6)) 
                        else piToTime(timing_param[[i]],type[i])
                    }))
    
    res <- data.frame(type=factor(type, levels=c("Mono-allelic Gain","CN-LOH","Bi-allelic Gain (WGD)")), time=time)
    colnames(res) <- c("type","time","time.lo","time.up","time.2nd","time.2nd.lo","time.2nd.up")
    
    # posthoc adjustment of CI's
    res$time.up <- (pseudo.count + bb$n.snv_mnv * res$time.up)/(pseudo.count + bb$n.snv_mnv)
    res$time.lo <- (0 + bb$n.snv_mnv * res$time.lo)/(pseudo.count + bb$n.snv_mnv)
    res$time.2nd.up <- (pseudo.count + bb$n.snv_mnv * res$time.2nd.up)/(pseudo.count + bb$n.snv_mnv)
    res$time.2nd.lo <- (0 + bb$n.snv_mnv * res$time.2nd.lo)/(pseudo.count + bb$n.snv_mnv)
    
    res$time.star <- factor((covrg == 1) + (min < 2 & maj <= 2 | min==2 & maj==2) * (covrg == 1), levels=0:2, labels = c("*","**","***")) ## ***: 2+0, 2+1, 2+2; **: 2<n+1, {3,4}+2; *: subclonal gains
    res$time.star[is.na(res$time)] <- NA
    return(res)
}

dtrbinom <- function(x, size, prob, xmin=0) dbinom(x,size,prob) / pbinom(xmin-1, size, prob, lower.tail=FALSE)
pbetabinom <- function(x, size, prob, rho){
    if(rho!=0)
        VGAM::pbetabinom(x, size, prob, rho)
    else
        pbinom(x, size, prob)
}

dbetabinom <- function(x, size, prob, rho){
    if(rho!=0)
        VGAM::dbetabinom(x, size, prob, rho)
    else
        dbinom(x, size, prob)
}

dtrbetabinom <- function(x, size, prob, rho, xmin=0) {y <- dbetabinom(x,size,prob,rho) / (1-pbetabinom(xmin-1, size, prob, rho))
    y[x<xmin] <- 0
    return(y)}

ptrbetabinom <- function(x, size, prob, rho, xmin=0) {
    pmin <- pbetabinom(xmin-1, size, prob, rho)
    pmax(0,(pbetabinom(x,size,prob,rho) - pmin) / (1-pmin))}

defineMcnStates <- function(bb, clusters, purity, gender='female', isWgd= FALSE){
    P <- vector(mode='list', length(bb))
    uniqueBB <- unique(bb)
    overlaps <- findOverlaps(uniqueBB, bb)
    
    majorCN <- split(bb$major_cn[subjectHits(overlaps)], queryHits(overlaps))
    m <- bb$minor_cn #hack: minor_cn > 0 in male samples - Battenberg bug?
    if(gender=='male')
        m[as.character(seqnames(bb)) %in% c('X','Y')] <- 0
    minorCN <- split(m[subjectHits(overlaps)], queryHits(overlaps)) 
    h <- selectHits(overlaps, "first")
    H <- selectHits(overlaps, "last")
    
    cnNormal <- 2 - (gender=='male' & seqnames(uniqueBB)=="X" | seqnames(uniqueBB)=="Y")
    
    # Fix cluster and purity discrepancies
    clusters$proportion[which.max(clusters$proportion)] <- purity
    
    cloneFreq <- split(bb$clonal_frequency[subjectHits(overlaps)], queryHits(overlaps))
    cnStates <- matrix(0, nrow=10000, ncol=6)
    colnames(cnStates) <- c("state","m","f","n.m.s","pi.m.s","s")
    
    power.c <- rep(0, nrow(clusters))
    
    deltaFreq <- 0.05 # merge clusters withing deltaFreq
    
    
    for( i in seq_along(uniqueBB)){
        if(!i %in% names(majorCN)) next
        majcni <- majorCN[[as.character(i)]]
        mincni <- minorCN[[as.character(i)]]
        cfi <- cloneFreq[[as.character(i)]]
        effCnTumor <- sum((majcni + mincni)*cfi)
        effCnNormal <- as.numeric(cnNormal[i]) * (1-purity)
        
        if(any(is.na(majcni))) next
        
        mixFlag <- FALSE
        subclonalGainFlag <- FALSE
        clonalFlag <- TRUE
        majdelta <- 0
        mindelta <- 0
        
        majanc <- majder <- majcni
        minanc <- minder <- mincni
        
        if(length(cfi)>1){ # multiple (subclonal) CN states, if so add clonal option (ie. mixture of both states), subclonal states only change by 1..delta(CN)
            d <- colSums(abs(rbind(majcni, mincni) - c(1,1) * (1+ isWgd)))
            derived <- d == max(d) ## derived state further from diploid/tetraploid
            if(all(derived)) next 
            majanc <- majcni[!derived]
            minanc <- mincni[!derived]
            majder <- majcni[derived]
            minder <- mincni[derived]
            majdelta <- majcni[derived] - majcni[!derived]
            mindelta <- mincni[derived] - mincni[!derived]
            majcni <- c(min(majcni), # clonal, sub on allele that doesn't change
                    (majcni[!derived] + cfi[derived]/purity*majdelta), # clonal, sub on allele that does change
                    (majcni[derived] >0) + (majdelta > 0)) # subclonal, subs after or before CNA, m=1,delta
            mincni <- c(min(mincni), (mincni[!derived] + cfi[derived]/purity*mindelta), (mincni[derived] >0) + (mindelta > 0))
            majdelta <- c(0, cfi[derived]/purity*majdelta, majdelta)
            mindelta <- c(0, cfi[derived]/purity*mindelta, mindelta)
            cfi <- c(purity, purity,  cfi[derived])
            mixFlag <- c(FALSE, TRUE, FALSE)
            clonalFlag <- c(TRUE,TRUE,FALSE)
            subclonalGainFlag <- c(FALSE, FALSE, TRUE)
        }
        
        
        a <- sapply(clusters$proportion, function(p) all(abs(p-cfi) > deltaFreq)) # subclone(s) not coinciding with CN change
        if(any(a)){ # assume subclones have derived from most abundant CN state
            majcni <- c(majcni, rep(majcni[which.max(cfi)]>0, sum(a))+0)
            mincni <- c(mincni, rep(mincni[which.max(cfi)]>0, sum(a))+0)
            cfi <- c(cfi, clusters$proportion[a])
            mixFlag <- c(mixFlag, rep(FALSE, sum(a)))
            clonalFlag <- c(clonalFlag, rep(FALSE, sum(a)))
            subclonalGainFlag <- c(subclonalGainFlag, rep(FALSE, sum(a)))
            majdelta <- c(majdelta, rep(0,sum(a)))
            mindelta <- c(mindelta, rep(0,sum(a)))
        }
        totcni <- majcni+mincni
        if(all(totcni==0)) next
        
        k <- 0
        for( j in seq_along(majcni)){
            if(majcni[j]==0 & mincni[j]==0) {
                f <- m <- 0 # allele frequency
                pi.m.s <- n.m.s <- 1
                l <- 1
            }else{
                l <- 1:max(majcni[j], mincni[j]) # mincni>majcni can occur if minor allele changes subclonally
                m <- l
                n.m.s <- rep(1, length(l)) #multiplicity of cn states
                if(!mixFlag[j]){ # single subclone, or no subclonal cn
                    f <- m * cfi[j] / (effCnTumor + effCnNormal)
                    if(mincni[j] > 0)
                        n.m.s[1:min(majcni[j], mincni[j])] <- 2
                    pi.m.s <- n.m.s/sum(n.m.s)
                }else{ # coexisting cn subclones, use mixture
                    d <- if(majdelta[j] != 0) majdelta[j] else mindelta[j] 
                    M <- max(mincni[j]*(mindelta[j]!=0), majcni[j]*(majdelta[j]!=0))
                    m <- (d > 0):M + (M-floor(M))
                    l <- seq_along(m)
                    f <- m *cfi[j] / (effCnTumor + effCnNormal) # Variant on major allele
                    n.m.s <- rep(1, length(l))
                    pi.m.s <- n.m.s/sum(n.m.s)
                }
            }
            cnStates[k + l,"state"]=j
            cnStates[k + l,"m"]=m
            cnStates[k + l,"f"]=f
            cnStates[k + l,"pi.m.s"]=pi.m.s
            cnStates[k + l,"n.m.s"]=n.m.s
            k <- k + length(l)
        }
        whichStates <- (1:k)[cnStates[1:k,"f"]>0]
        
        # State probabilities - based on cell fractions
        cfi.s <- unique(cfi)
        pi.s <- sapply(cfi.s, function(p) ifelse(min(abs(clusters$proportion - p)) < deltaFreq, clusters$n_ssms[which.min(abs(clusters$proportion - p))], 1))
        pi.s <- pi.s/sum(pi.s)
        
        c.to.s <- sapply(cfi.s, function(p) ifelse(min(abs(clusters$proportion - p)) < deltaFreq, which.min(abs(clusters$proportion - p)), NA)) # map to cluster
        s.to.c <- sapply(clusters$proportion, function(c) which.min(abs(cfi.s - c)))
        
        cnStates[1:k,"s"] = as.numeric(factor(cfi, levels=cfi.s))[cnStates[1:k,"state"]]
        
        timing_param <- cbind(cnStates[whichStates,,drop=FALSE], cfi=cfi[cnStates[whichStates,"state"]], pi.s=pi.s[cnStates[whichStates,"s"]], P.m.sX=NA,P.m.sX.lo=NA, P.m.sX.up=NA, T.m.sX=NA, T.m.sX.lo=NA, T.m.sX.up=NA, power.s=NA, power.m.s = NA,
                majCN=majcni[cnStates[whichStates,"state"]], minCN=mincni[cnStates[whichStates,"state"]], 
                majDelta = majdelta[cnStates[whichStates,"state"]], minDelta = mindelta[cnStates[whichStates,"state"]], 
                clonalFlag=clonalFlag[cnStates[whichStates,"state"]], subclonalGainFlag=subclonalGainFlag[cnStates[whichStates,"state"]], mixFlag=mixFlag[cnStates[whichStates,"state"]], majCNanc=majanc, minCNanc=minanc, majCNder=majder, minCNder=minder)
        P[[h[i]]] <- timing_param
        if(H[i] != h[i]) P[[H[[i]]]] <- P[[h[i]]]
        
    }
    return(P)
} 

computeMutCn <- function(vcf, bb, clusters, purity, gender='female', isWgd= FALSE, xmin=3, rho=0, n.boot=200){
    n <- nrow(vcf)
    D <- DataFrame(MutCN=rep(NA,n), MutDeltaCN=rep(NA,n), MajCN=rep(NA,n), MinCN=rep(NA,n), MajDerCN=rep(NA,n), MinDerCN=rep(NA,n), CNF=rep(NA,n), CNID =as(vector("list", n),"List"), pMutCN=rep(NA,n), pGain=rep(NA,n),pSingle=rep(NA,n),pSub=rep(NA,n), pMutCNTail=rep(NA,n))    
    P <- defineMcnStates(bb,clusters, purity, gender, isWgd)
    if(n==0)
        return(list(D=D, P=P))
    
    altCount <- getAltCount(vcf)
    tumDepth <- getTumorDepth(vcf)
    names(altCount) <- names(tumDepth) <- NULL
    
    # Match VCF and CN
    overlaps <- findOverlaps(vcf, bb)
    D[["CNID"]] <- as(overlaps, "List")
    majorCN <- split(bb$major_cn[subjectHits(overlaps)], queryHits(overlaps))
    m <- bb$minor_cn #hack: minor_cn > 0 in male samples - Battenberg bug?
    if(gender=='male')
        m[as.character(seqnames(bb)) %in% c('X','Y')] <- 0
    minorCN <- split(m[subjectHits(overlaps)], queryHits(overlaps)) 
    h <- selectHits(overlaps, "first")
    H <- selectHits(overlaps, "last")
    
    cnNormal <- 2 - (gender=='male' & seqnames(vcf)=="X" | seqnames(vcf)=="Y")
    
    # Fix cluster and purity discrepancies
    clusters$proportion[which.max(clusters$proportion)] <- purity
    
    cloneFreq <- split(bb$clonal_frequency[subjectHits(overlaps)], queryHits(overlaps))
    
    power.c <- rep(0, nrow(clusters))
    
    deltaFreq <- 0.05 # merge clusters withing deltaFreq
    
    boundaryHits <- countOverlaps(vcf, unique(bb)) > 1 # indels overlapping with CN boundaries
    
    for(globalIt in 1:2){ # 2 iterations, fist local (ie per segment) fit, then global
        for( i in which( (diff(c(-1, h)) !=0 | is.na(diff(c(-1, h)) !=0) ) & ! boundaryHits )){
            if(!i %in% names(majorCN)) next
            if(is.na(h[i])) next
            if(if(i>1) h[i] != h[i-1] | is.na(h[i-1]) else TRUE){ #ie. new segment
                
                hh <- which(h==h[i] & !is.na(altCount) &! is.na(tumDepth))
                if(length(hh)==0) next
                
                if(is.null(bb$timing_param[[h[i]]])){
                    cnStates <- P[[h[i]]]
                    if(is.null(cnStates)) next
                    whichStates <- 1:nrow(cnStates)
                    
                    # State probabilities - based on cell fractions
                    cfi.s <- unique(cnStates[,"cfi"])
                    pi.s <- sapply(cfi.s, function(p) ifelse(min(abs(clusters$proportion - p)) < deltaFreq, clusters$n_ssms[which.min(abs(clusters$proportion - p))], 1))
                    pi.s <- pi.s/sum(pi.s)
                    
                    c.to.s <- sapply(cfi.s, function(p) ifelse(min(abs(clusters$proportion - p)) < deltaFreq, which.min(abs(clusters$proportion - p)), NA)) # map to cluster
                    s.to.c <- sapply(clusters$proportion, function(c) which.min(abs(cfi.s - c)))
                    
                    
                    # Likelihood
                    L <- matrix(sapply(pmin(cnStates[whichStates,"f"],1), function(pp) dtrbetabinom(altCount[hh],tumDepth[hh],pp, rho=rho, xmin=pmin(altCount[hh],0)) + .Machine$double.eps), ncol=length(whichStates))
                    
                    # Power
                    power.sm <- colMeans(matrix(sapply(pmin(cnStates[whichStates,"f"],1), function(pp) 1-ptrbetabinom(pmin(altCount[hh],xmin-1),tumDepth[hh],pp, rho=rho, xmin=0)), ncol=length(whichStates)), na.rm=TRUE)
                    if(globalIt==2){
                        P.m.sX <- P[[h[i]]][,"P.m.sX"]
                        power.s <- sapply(split(power.sm * P.m.sX, cnStates[whichStates,"s"]), sum) # Power for state
                        power.s[!is.na(c.to.s)] <- power.c[na.omit(c.to.s)]
                        power.m.s <- power.sm # Relative power of m states (local) conditioned on s (global).
                        for(s in unique(cnStates[whichStates,"s"])) power.m.s[cnStates[whichStates,"s"]==s] <- power.m.s[cnStates[whichStates,"s"]==s] / max(power.m.s[cnStates[whichStates,"s"]==s]) #power.s[s]
                    }else{ # First iteration
                        power.m.s <- power.sm
                        power.s <- rep(1,length(whichStates))
                    }
                    
                    mm <- function(x) {
                        x <- factor(x)
                        if(nlevels(x) > 1) t(model.matrix( ~ x + 0)) else matrix(1, ncol=length(x))
                    }
                    
                    # EM algorithm (mixture fitting) for pi
                    P.m.sX <- cnStates[whichStates,"pi.m.s"]
                    s.from.m <- mm(cnStates[whichStates,"s"]) # indicator matrix to map
                    for(em.it in 1:100){
                        P.xsm <- L * rep(pi.s[cnStates[whichStates,"s"]] * P.m.sX / power.m.s / power.s[cnStates[whichStates,"s"]], each=nrow(L)) # P(X,s,m)
                        P.sm.x <- P.xsm/rowSums(P.xsm) # P(s,m|Xi)
                        P.sm.X <- pmax(.Machine$double.xmin,colMeans(P.sm.x)) # P(s,m|X) / piState[cnStates[1:k,"state"]] / cnStates[1:k,"pi.m.s"]
                        if(em.it==100) break
                        P.s.X <- s.from.m %*% P.sm.X 
                        P.m.sX <- P.sm.X / P.s.X[cnStates[whichStates,"s"]]
                    }
                    
                    toTime <- function(cnStates, P.m.sX, s.m) {
                        mAnc <- cnStates[,"m"] - cnStates[,"minDelta"] - cnStates[,"majDelta"]
                        mAnc.s <- factor(paste(mAnc, cnStates[,"s"], sep="."))
                        n <- (mAnc <= cnStates[,"majCNanc"]) + (mAnc <= cnStates[,"minCNanc"] )
                        mAnc.s.from.m <- mm(x = mAnc.s)# indicator matrix to map
                        return((mAnc.s.from.m[mAnc.s,] %*% P.m.sX) / (s.m[cnStates[,"s"],] %*% (P.m.sX * mAnc)) *  (cnStates[,"majCNanc"] + cnStates[,"minCNanc"]) / n)
                    }
                    
                    T.m.sX <- toTime(cnStates, P.m.sX, s.from.m) 
                    
                    if(globalIt==1){
                        p <- (sapply(split(power.sm * P.m.sX, cnStates[whichStates,"s"]), sum) * nrow(L))[s.to.c]
                        if(!any(is.na(p) | is.nan(p)))
                            power.c <- power.c + p 
                    }
                    
                    # Bootstrapping for CIs
                    if(globalIt==2){
                        b.m.sX <- if(n.boot>0) sapply(1:n.boot, function(foo){
                                                L <- rbind(L, rep(1e-3, each=ncol(L))) #add an uniformative row
                                                L <- L[sample(1:nrow(L), replace=TRUE),,drop=FALSE]
                                                P.m.sX <- cnStates[whichStates,"pi.m.s"]
                                                for(em.it in 1:50){
                                                    P.xsm <- L * rep(pi.s[cnStates[whichStates,"s"]] * P.m.sX / power.m.s / power.s[cnStates[whichStates,"s"]], each=nrow(L)) # P(X,s,m)
                                                    P.sm.x <- P.xsm/rowSums(P.xsm) # P(s,m|Xi)
                                                    P.sm.X <- colMeans(P.sm.x) # P(s,m|X) / piState[cnStates[1:k,"state"]] / cnStates[1:k,"pi.m.s"]
                                                    P.s.X <- s.from.m %*% P.sm.X 
                                                    P.m.sX <- P.sm.X / P.s.X[cnStates[whichStates,"s"]]
                                                }
                                                return(P.m.sX)
                                            }) else NA
                        if(n.boot>0) try({
                                        CI.m.sX <- apply(b.m.sX, 1, quantile, c(0.025, 0.975))
                                        cnStates[,"P.m.sX.lo"] <- CI.m.sX[1,] 
                                        cnStates[,"P.m.sX.up"] <- CI.m.sX[2,]
                                        B.m.sX <- toTime(cnStates = cnStates, P.m.sX = b.m.sX, s.m = s.from.m)
                                        C.m.sX <- apply(B.m.sX, 1, quantile, c(0.025, 0.975))
                                        cnStates[,"T.m.sX.lo"] <- C.m.sX[1,] 
                                        cnStates[,"T.m.sX.up"] <- C.m.sX[2,]
                                    }, silent=TRUE)
                    }
                    
                    P.sm.x[apply(is.na(P.sm.x)|is.nan(P.sm.x),1,any),] <- NA
                    cnStates[,"P.m.sX"] <- P.m.sX
                    cnStates[,"T.m.sX"] <- T.m.sX
                    cnStates[,"power.s"] <- power.s[cnStates[whichStates,"s"]]
                    cnStates[,"power.m.s"] <- power.m.s
                    
                }else{
                    cnStates <- bb$timing_param[[h[i]]]
                    if(is.null(cnStates)) next
                    P.sm.x <- posteriorMutCN(x=altCount[hh],n=tumDepth[hh], cnStates=cnStates, xmin=0, rho=rho)
                }
                
                # Tail probs
                pMutCNTail <- matrix(sapply(pmin(cnStates[,"f"],1), function(pp) ptrbetabinom(altCount[hh],tumDepth[hh],pp, rho=rho, xmin=pmin(altCount[hh],xmin))), ncol=nrow(cnStates)) #%*% c(pi.s[cnStates[whichStates,"state"]] * P.m.sX)              
                
                P[[h[i]]] <- cnStates
                if(H[i] != h[i]) P[[H[[i]]]] <- P[[h[i]]]
                
                w <- apply(P.sm.x, 1, function(x) if(any(is.na(x))) NA else which.max(x) )
                if(all(is.na(w))) next
                
                D[hh, "pSub"] <- rowSums(P.sm.x[, !cnStates[,"clonalFlag"], drop=FALSE])
                D[hh, "pGain"] <- rowSums(P.sm.x[, cnStates[,"clonalFlag"] & cnStates[,"m"] > 1.00001 + cnStates[,"majDelta"] + cnStates[,"minDelta"], drop=FALSE])
                #D[hh, "pSingle"] <- rowSums(P.sm.x[, cnStates[1:k,"state"] %in% which(clonalFlag) & cnStates[1:k,"m"]<=1, drop=FALSE])
                D[hh, "pSingle"] <-  1 - D[hh, "pSub"] - D[hh, "pGain"]         
                
                D[hh,"MutCN"]  <- cnStates[w,"m"]
                D[hh,"MutDeltaCN"]  <- cnStates[w,"majDelta"] + cnStates[w,"minDelta"]
                D[hh,"MinCN"] <- cnStates[w,"minCNanc"]
                D[hh,"MajCN"] <- cnStates[w,"majCNanc"]
                D[hh,"MinDerCN"] <- cnStates[w,"minCNder"]
                D[hh,"MajDerCN"] <- cnStates[w,"majCNder"]
                
                D[hh,"CNF"]  <- cnStates[w,"cfi"]
                D[hh,"pMutCN"] <- sapply(seq_along(w), function(i) P.sm.x[i,w[i]])
                D[hh,"pMutCNTail"] <- sapply(seq_along(w), function(i) pMutCNTail[i,w[i]])
            }       
        }
        if(globalIt==1){
            power.c <- power.c / sum(!is.na(D[,"pSub"]))
            if(any(is.na(power.c) | power.c==0)) {
                break # Cancel 2nd iteration
                warning("Power calculation yielded NA, aborting.")
            }
            if(any(power.c > 1)) {
                warning("Calculated power > 1, rounding down.")
                power.c <- pmin(1, power.c)
            }
        }
    }
    return(list(D=D,P=P, power.c=power.c))
}

averagePloidy <- function(bb) {
    c <- if(!is.null(bb$copy_number)) bb$copy_number else bb$total_cn
    sum(width(bb) * c * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
}

averageHom <- function(bb){
    sum(width(bb) * (bb$minor_cn == 0) * bb$clonal_frequency, na.rm=TRUE) / sum(width(bb) * bb$clonal_frequency, na.rm=TRUE)
}

.classWgd <- function(ploidy, hom) 2.9 -2*hom <= ploidy

fractionGenomeWgdCompatible <- function(bb, min.dist=0.05){
    m <- findMainCluster(bb)
    l <- pmin(bb$time.lo, bb$time - min.dist)
    u <- pmax(bb$time.up, bb$time + min.dist)
    w <- which(l <= m & u >= m)
    avgCi <- weighted.mean(bb$time.up- bb$time.lo, width(bb), na.rm=TRUE)
    sd.wgd <- sqrt(weighted.mean((bb$time[w] - m)^2, width(bb)[w], na.rm=TRUE))
    sd.all <- sqrt(weighted.mean((bb$time - m)^2, width(bb), na.rm=TRUE))
    c(nt.wgd=sum(as.numeric(width(bb))[w]), nt.total=sum(as.numeric(width(bb))[!is.na(bb$time)]), time.wgd=m, n.wgd=length(w), n.all = sum(!is.na(bb$time)), chr.wgd = length(unique(seqnames(bb)[w])), chr.all = length(unique(seqnames(bb)[!is.na(bb$time)])), sd.wgd=sd.wgd, avg.ci=avgCi, sd.all=sd.all) 
}
