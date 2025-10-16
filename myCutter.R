library(lpSolve)
library(RSpectra)

n <- 100
m <- 15
set.seed(12345)
w <- crossprod(matrix(rnorm(n * m), n, m)) / n

myCutter <- function(w,
                     f.obj = rep(1, nrow(w)),
                     itmax = 1000,
                     eps = 1e-6,
                     verbose = TRUE) {
  m <- nrow(w)
  f.con <- diag(m)
  f.rhs <- diag(w)
  d.old <- f.rhs
  f.dir <- rep(">=", m)
  e.old <- sum(f.obj * d.old)
  itel <- 1
  repeat {
    dd <- diag(d.old)
    df <- dd - w
    ev <- eigs_sym(df, 1, which = "SA")
    lv <- ev$values
    xx <- drop(ev$vectors)
    f.con <- rbind(f.con, xx^2)
    f.rhs <- c(f.rhs, sum(xx * (w %*% xx)))
    f.dir <- c(f.dir, ">=")
    ll <- lp("min", f.obj, f.con, f.dir, f.rhs)
    d.new <- ll$solution
    e.new <- ll$objval
    if (verbose) {
      cat(
        "itel",
        formatC(itel, width = 4, format = "d"),
        "oldval",
        formatC(
          e.old,
          width = 10,
          digits = 6,
          format = "f"
        ),
        "newval",
        formatC(
          e.new,
          width = 10,
          digits = 6,
          format = "f"
        ),
        "eigval",
        formatC(
          lv,
          width = 10,
          digits = 6,
          format = "f"
        ),
        "\n"
      )
    }
    if ((itel == itmax) || (abs(e.new - e.old) < eps)){
      break
    }
    
    itel <- itel + 1
    e.old <- e.new
    d.old <- d.new
  }
  return(list(f = e.new, d = d.new, itel = itel))
}
