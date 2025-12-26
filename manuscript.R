require(mvtnorm)
source("utils.R")

power.approx <- function(n, p, group.size, alpha,
                         reps, sampler,
                         mu.x, mu.1, mu.2s, Sigma, var.effect) {
  k <- p - 1
  sapply(mu.2s, function(mu.2) {
    mu.y <- c(mu.1, mu.2)
    n.x <- n.y <- n

    p.vals <- replicate(reps, {
      x <- sampler(n.x, mu.x, Sigma, var.effect)
      y <- sampler(n.y, mu.y, Sigma, var.effect)

      xy <- cbind(x, y)
      d <- rep(c(0, 1), c(n.x, n.y))
      G <- (n.x + n.y) / group.size
      g <- sample(gl(G, group.size))

      latent <- t(rmvnorm(G, sigma = diag(p) * var.effect))
      xy <- xy + latent[, g]

      x.final <- xy[, d == 0]
      y.final <- xy[, d == 1]
      mu.x.hat <- rowMeans(x.final)
      mu.y.hat <- rowMeans(y.final)
      Sigma.x.hat <- (n.x - 1) / n.x * var(t(x.final))
      Sigma.y.hat <- (n.y - 1) / n.y * var(t(y.final))

      robust <- maha.test.robust(xy, d, g, k, mu.x.hat, mu.y.hat, Sigma.x.hat, Sigma.y.hat)["p.val"]
      exact <- maha.test.exact(x.final, y.final, k)["p.val"]

      mu.hat <- mu.y.hat - mu.x.hat
      Sigma.hat.pooled <- (n.x * Sigma.x.hat + n.y * Sigma.y.hat) / (n.x + n.y)
      gamma.hat <- solve(Sigma.hat.pooled) %*% mu.hat
      beta.hat <- solve(Sigma.hat.pooled[1:k, 1:k]) %*% mu.hat[1:k]

      delong <- suppressMessages(
        pROC::roc.test(d, c(t(beta.hat) %*% xy[1:k, ]), c(t(gamma.hat) %*% xy),
          method = "delong", smooth = FALSE
        )$p.value
      )

      c(robust = robust, exact = exact, delong = delong)
    })
    rowMeans(p.vals <= alpha)
  })
}

## launcher
power.sim <- function(n, p, group.size, alpha,
                      reps, resolution, sampler, param.to.delta,
                      rho, effect.scale, init.mu.sd) {
  cat(n, group.size, rho, effect.scale, reps, "\n", sep = "/")
  set.seed(0)
  stopifnot(n %% group.size == 0)
  k <- p - 1

  Sigma <- with(list(S = matrix(runif(p^2), p)), S %*% t(S))
  Sigma <- diag(1 / sqrt(diag(Sigma))) %*% Sigma %*% diag(1 / sqrt(diag(Sigma)))
  var.effect <- rho / (1 - rho)
  Sigma <- Sigma + diag(p) * var.effect

  Sigma.11 <- Sigma[1:k, 1:k]
  Sigma.21 <- Sigma[(k + 1):p, 1:k]

  mu.x <- rep(0, p)
  mu.1 <- rnorm(k) * init.mu.sd
  mu.2.null <- c((Sigma.21) %*% solve(Sigma.11) %*% mu.1)

  auc.diffs <- seq(0, 1, len = resolution) / effect.scale
  mu.2s <- sapply(auc.diffs, function(auc.diff) {
    obj <- function(mu.2) (param.to.delta(c(mu.1, mu.2), Sigma) - auc.diff)^2
    interval <- if (auc.diff < 0) c(mu.2.null - 1, mu.2.null) else c(mu.2.null, mu.2.null + 1)
    opt0 <- optimize(obj, interval)
    if (opt0$obj > 1e-8) {
      return(NA)
    }
    opt0$min
  })
  stopifnot(!anyNA(mu.2s))

  power <- power.approx(
    n = n, p = p, group.size = group.size, alpha = alpha,
    reps = reps, sampler = sampler,
    mu.x = mu.x, mu.1 = mu.1, mu.2s = mu.2s,
    Sigma = Sigma, var.effect = var.effect
  )

  data.frame(t(power), auc.diffs = auc.diffs)
}



plot.sim <- function(by.settings, alpha, filename = NULL) {
  # if (!is.null(filename)) png(filename, width = 1024, height = 768, pointsize = 15, family = "CM Roman")
  if (!is.null(filename)) tikzDevice::tikz(filename, width = 6.5, height = 7)
  op <- par(mfrow = c(3, 3))
  lapply(1:length(by.settings), function(i) {
    power <- by.settings[[i]]
    n <- settings$n[i]
    rho <- settings$rho[i]
    plot(power$auc.diffs, power$robust.p.val, type = "l", ylim = c(0, 1), xlab = "", ylab = "", main = paste0("n=", n, ", $\\rho=$", rho)) # ,xlab='difference in AUCs',ylab='power')
    lines(power$auc.diffs, power$exact.p.val, lty = 2) # ,type='o')
    lines(power$auc.diffs, power$delong, lty = 3)
    abline(h = alpha, lty = 4)
  })
  mtext("power", side = 2, line = -2, outer = TRUE)
  mtext("difference in AUCs", side = 1, line = -2, outer = TRUE)
  par(op)
  if (!is.null(filename)) dev.off()
}



alpha <- .1
p <- 5
init.mu.sd <- .2
reps <- 1e3
resolution <- 30
rho <- c(0, .3, .6)

## gaussian data
sampler <- function(n, mu, var.indiv, var.effect) {
  p <- nrow(var.indiv)
  t(rmvnorm(n, mu, var.indiv - diag(p) * var.effect))
}
param.to.delta <- auc.diff.normal
settings <- data.frame(n = rep(c(5e1, 2e2, 3e2), each = 3), group.size = 25, effect.scale = rep(c(10, 20, 20), each = 3), latent.signal = c(0, .2, .4), rho = rho)
by.settings <- apply(settings, 1, function(c) {
  power.sim(
    n = c["n"],
    p = p,
    group.size = c["group.size"],
    alpha = alpha,
    reps = reps,
    resolution = resolution,
    sampler = sampler,
    param.to.delta = param.to.delta,
    rho = c["rho"],
    effect.scale = c["effect.scale"],
    init.mu.sd = init.mu.sd
  )
}, simplify = FALSE)
plot.sim(by.settings, alpha = alpha, filename = "figs/power_normal.tex")
## plot.sim(by.settings,alpha=alpha)


## Student's t
nu <- 5
sampler <- function(n, mu, var.indiv, var.effect) {
  p.dim <- nrow(var.indiv)
  t(rmvt(
    n = n, delta = mu,
    sigma = (nu - 2) / nu * (var.indiv - diag(p.dim) * var.effect),
    df = nu
  ))
}
param.to.delta <- function(mu, Sigma) auc.diff.t(mu, Sigma, nu)
settings <- data.frame(n = rep(c(5e1, 2e2, 3e2), each = 3), group.size = 25, effect.scale = rep(c(10, 20, 100), each = 3), latent.signal = c(0, .2, .4), rho = rho)
by.settings <- apply(settings, 1, function(c) {
  power.sim(
    n = c["n"],
    p = p,
    group.size = c["group.size"],
    alpha = alpha,
    reps = reps,
    resolution = resolution,
    sampler = sampler,
    param.to.delta = param.to.delta,
    rho = c["rho"],
    effect.scale = c["effect.scale"],
    init.mu.sd = init.mu.sd
  )
}, simplify = FALSE)
## plot.sim(by.settings, alpha = .1)
plot.sim(by.settings, alpha = alpha, filename = "figs/power_t.tex")
