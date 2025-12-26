maha.test.robust <- function(xy,d,g,k,mu.x,mu.y,Sigma.x,Sigma.y) {
    x <- xy[,d==0]
    y <- xy[,d==1]
    p <- nrow(x)
    n.x <- ncol(x); n.y <- ncol(y)
    G <- nlevels(g)
    group.size <- (n.x+n.y)/G
    stopifnot(all.equal(group.size-floor(group.size), 0))
    x.1 <- x[1:k,]; y.1 <- y[1:k,]
    x.2 <- x[(k+1):p,,drop=FALSE]; y.2 <- y[(k+1):p,,drop=FALSE];
    mu.x.1 <- mu.x[1:k]; mu.y.1 <- mu.y[1:k]
    mu.x.2 <- mu.x[(k+1):p]; mu.y.2 <- mu.y[(k+1):p]
    Sigma.x.11 <- Sigma.x[1:k,1:k]; Sigma.x.21 <- Sigma.x[(k+1):p,1:k]; Sigma.x.22 <- Sigma.x[(k+1):p,(k+1):p]
    Sigma.y.11 <- Sigma.y[1:k,1:k]; Sigma.y.21 <- Sigma.y[(k+1):p,1:k]; Sigma.y.22 <- Sigma.y[(k+1):p,(k+1):p]
    Sigma <- Sigma.x+Sigma.y
    Sigma.11 <- Sigma[1:k,1:k]; Sigma.21 <- Sigma[(k+1):p,1:k]; Sigma.22 <- Sigma[(k+1):p,(k+1):p]
    A.y <- (Sigma.x.21+Sigma.y.21+mu.y.2%*%t(-mu.y.1+mu.x.1))%*%solve(Sigma.11)
    q.y <- c(t(y.1-mu.y.1)%*%solve(Sigma.11)%*%(mu.y.1-mu.x.1))
    A.x <- (Sigma.x.21+Sigma.y.21+mu.x.2%*%t(mu.y.1-mu.x.1))%*%solve(Sigma.11)
    q.x <- c(t(x.1-mu.x.1)%*%solve(Sigma.11)%*%(mu.y.1-mu.x.1))
    infl.y <- (n.x+n.y)/n.y * (y.2 - t(t(y.2)*q.y) - A.y%*%(y.1-mu.y.1) + Sigma.21%*%solve(Sigma.11)%*%( t(t(y.1-mu.y.1)*q.y) - c(Sigma.y.11%*%solve(Sigma.11)%*%(mu.y.1-mu.x.1)) ) )
    infl.x <- (n.x+n.y)/n.x * (-x.2 - t(t(x.2)*q.x) + A.x%*%(x.1-mu.x.1) + Sigma.21%*%solve(Sigma.11)%*%( t(t(x.1-mu.x.1)*q.x) - c(Sigma.x.11%*%solve(Sigma.11)%*%(mu.y.1-mu.x.1)) ) )
    infl <- matrix(nrow=p-k,ncol=n.x+n.y)
    infl[,d==0] <- infl.x; infl[,d==1] <- infl.y
    group.sums <- lapply(split(as.data.frame(t(infl)),g), colSums)
    group.sums <- do.call(cbind,group.sums)
    mu.infl <- rowMeans(infl)
    var.hat <- (group.sums%*%t(group.sums) / (n.x+n.y) - group.size*mu.infl%*%t(mu.infl)) / (n.x+n.y)
    theta <- mu.y.2 - mu.x.2 - (Sigma.21)%*%solve(Sigma.11)%*%(mu.y.1-mu.x.1)
    z.stat <- expm::sqrtm(solve(var.hat))%*%theta
    test.stat <- sum(z.stat^2)
    p.val <- pchisq(test.stat,df=p-k,lower.tail=FALSE)
    c(test.stat=test.stat,p.val=p.val,df=p-k)
}

maha.test.exact <- function(x,y,k) {
    p <- nrow(x)
    n.x <- ncol(x); n.y <- ncol(y)
    m <- n.x+n.y-2
    df1 <- p-k; df2 <- m-p+1
    M <- ((n.x-1)*var(t(x)) + (n.y-1)*var(t(y))) * (1/n.x+1/n.y)
    u <- rowMeans(y) - rowMeans(x)
    D2.p <- m*t(u)%*%solve(M)%*%u
    D2.k <- m*t(u[1:k])%*%solve(M[1:k,1:k])%*%(u[1:k])
    test.stat <- (m-p+1)/(p-k)*(D2.p-D2.k)/(m+D2.k)
    p.val <- pf(test.stat,df1,df2,lower.tail=FALSE)
    return(c(test.stat=test.stat,p.val=p.val,df1=df1,df2=df2))
} 

## diff between full AUC and AUC after removing last covariate--gaussian
auc.diff.normal <- function(mu,Sigma) {
    p <- length(mu)
    k <- p-1
    mu.1 <- mu[1:k]; mu.2 <- mu[(k+1):p]
    Sigma.11 <- Sigma[1:k,1:k]
    delta.1 <- t(mu.1)%*%solve(Sigma.11)%*%mu.1
    auc.reduced <- pnorm(0,mean=delta.1,sd=sqrt(delta.1))
    delta <- t(mu)%*%solve(Sigma)%*%mu
    auc.full <- pnorm(0,mean=delta,sd=sqrt(delta))
    return(abs(auc.full-1/2)-abs(auc.reduced-1/2))
    auc.full - auc.reduced
}


## diff between full AUC and AUC after removing last covariate--Students t
auc.diff.t <- function(mu,Sigma,nu) {
    auc.t <- function(mu,Sigma,nu) {
        Delta <- c(t(mu)%*%solve(Sigma)%*%mu)
        f <- function(t)t^(nu-1)*besselK(t*sqrt((nu-2)*Delta),nu=nu/2)^2 * sin(-t*Delta)
        1/2 - 1/pi*((nu-2)*Delta)^(nu/2)*2^(2-nu)/gamma(nu/2)^2 * integrate(f, 0, Inf)$val
    }
    p <- length(mu)
    k <- p-1
    mu.1 <- mu[1:k]; mu.2 <- mu[(k+1):p]
    Sigma.11 <- Sigma[1:k,1:k]
    delta.1 <- t(mu.1)%*%solve(Sigma.11)%*%mu.1
    auc.reduced <- auc.t(mu.1,Sigma.11,nu)
    auc.full <- auc.t(mu,Sigma,nu)
    return(abs(auc.full-1/2)-abs(auc.reduced-1/2))
}
