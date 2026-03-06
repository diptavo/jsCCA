soft_threshold <- function(x, delta) {
  sign(x) * pmax(abs(x) - delta, 0)
}

weighted_soft_threshold <- function(x, delta, weights) {
  sign(x) * pmax(abs(x) - delta * weights, 0)
}

find_delta_for_l1 <- function(x, l1_bound, tol = 1e-8, max_iter = 100) {
  if (l1_bound <= 0) {
    return(max(abs(x)) + 1)
  }

  if (sum(abs(x)) <= l1_bound) {
    return(0)
  }

  lo <- 0
  hi <- max(abs(x))

  for (i in seq_len(max_iter)) {
    mid <- 0.5 * (lo + hi)
    s <- sum(abs(soft_threshold(x, mid)))

    if (abs(s - l1_bound) < tol) {
      return(mid)
    }

    if (s > l1_bound) {
      lo <- mid
    } else {
      hi <- mid
    }
  }

  0.5 * (lo + hi)
}

find_delta_for_l1_weighted <- function(x, l1_bound, weights, tol = 1e-8, max_iter = 100) {
  if (l1_bound <= 0) {
    return(max(abs(x)) + 1)
  }

  if (sum(abs(x)) <= l1_bound) {
    return(0)
  }

  w <- pmax(as.numeric(weights), .Machine$double.eps)
  lo <- 0
  hi <- max(abs(x) / w)

  for (i in seq_len(max_iter)) {
    mid <- 0.5 * (lo + hi)
    s <- sum(abs(weighted_soft_threshold(x, mid, w)))

    if (abs(s - l1_bound) < tol) {
      return(mid)
    }

    if (s > l1_bound) {
      lo <- mid
    } else {
      hi <- mid
    }
  }

  0.5 * (lo + hi)
}

sparse_unit_vector <- function(x, l1_bound, tol = 1e-8, max_iter = 100) {
  sparse_unit_vector_penalized(
    x = x,
    l1_bound = l1_bound,
    penalty = "l1",
    penalty_params = list(),
    tol = tol,
    max_iter = max_iter
  )
}

sparse_unit_vector_penalized <- function(
  x,
  l1_bound,
  penalty = c("l1", "mcp", "adaptive_lasso"),
  penalty_params = list(),
  tol = 1e-8,
  max_iter = 100
) {
  penalty <- match.arg(penalty)
  x <- as.numeric(x)

  if (!any(is.finite(x)) || all(abs(x) < .Machine$double.eps)) {
    return(rep(0, length(x)))
  }

  if (penalty == "l1") {
    delta <- find_delta_for_l1(x, l1_bound = l1_bound, tol = tol, max_iter = max_iter)
    sx <- soft_threshold(x, delta)
  } else if (penalty == "mcp") {
    gamma <- if (is.null(penalty_params$gamma)) 3 else as.numeric(penalty_params$gamma)
    if (gamma <= 1) {
      stop("For MCP penalty, gamma must be > 1.")
    }
    lambda <- max(0, as.numeric(l1_bound))
    ax <- abs(x)
    sx <- rep(0, length(x))
    idx_mid <- ax <= gamma * lambda
    denom <- (1 - 1 / gamma)
    sx[idx_mid] <- sign(x[idx_mid]) * pmax(ax[idx_mid] - lambda, 0) / denom
    sx[!idx_mid] <- x[!idx_mid]
  } else {
    # Adaptive lasso: heavier shrinkage on smaller coefficients.
    gamma <- if (is.null(penalty_params$gamma)) 1 else as.numeric(penalty_params$gamma)
    eps <- if (is.null(penalty_params$eps)) 1e-4 else as.numeric(penalty_params$eps)
    w <- 1 / ((abs(x) + eps)^gamma)
    delta <- find_delta_for_l1_weighted(x, l1_bound = l1_bound, weights = w, tol = tol, max_iter = max_iter)
    sx <- weighted_soft_threshold(x, delta, w)
  }

  nrm <- sqrt(sum(sx^2))

  if (nrm < .Machine$double.eps) {
    return(rep(0, length(x)))
  }

  sx / nrm
}

standardize_matrix <- function(X) {
  X <- as.matrix(X)
  X <- scale(X, center = TRUE, scale = TRUE)
  X[, !is.finite(colSums(X, na.rm = TRUE))] <- 0
  X[is.na(X)] <- 0
  X
}

jscca_component <- function(
  Zcm,
  Zme,
  lambda_u,
  lambda_w,
  lambda_v,
  penalty_u = c("l1", "mcp", "adaptive_lasso"),
  penalty_w = c("l1", "mcp", "adaptive_lasso"),
  penalty_v = c("l1", "mcp", "adaptive_lasso"),
  penalty_params_u = list(),
  penalty_params_w = list(),
  penalty_params_v = list(),
  u_init = NULL,
  w_init = NULL,
  v_init = NULL,
  tol = 1e-6,
  max_iter = 500,
  update_mode = c("paper", "alternating")
) {
  update_mode <- match.arg(update_mode)
  penalty_u <- match.arg(penalty_u)
  penalty_w <- match.arg(penalty_w)
  penalty_v <- match.arg(penalty_v)

  p <- nrow(Zcm)
  q <- ncol(Zcm)
  r <- ncol(Zme)

  if (nrow(Zme) != q) {
    stop("Dimension mismatch: ncol(Zcm) must equal nrow(Zme).")
  }

  if (is.null(u_init)) {
    u <- rep(0, p)
    u[which.max(rowSums(abs(Zcm)))] <- 1
  } else {
    u <- as.numeric(u_init)
  }

  if (is.null(w_init)) {
    w <- rep(0, q)
    w[which.max(rowSums(abs(Zme)))] <- 1
  } else {
    w <- as.numeric(w_init)
  }

  if (is.null(v_init)) {
    v <- rep(0, r)
    v[which.max(colSums(abs(Zme)))] <- 1
  } else {
    v <- as.numeric(v_init)
  }

  if (sqrt(sum(u^2)) > 0) u <- u / sqrt(sum(u^2))
  if (sqrt(sum(w^2)) > 0) w <- w / sqrt(sum(w^2))
  if (sqrt(sum(v^2)) > 0) v <- v / sqrt(sum(v^2))

  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    u_old <- u
    w_old <- w
    v_old <- v

    if (update_mode == "paper") {
      # As written in the appendix: each update uses (i-1) iterate.
      w <- sparse_unit_vector_penalized(
        t(Zcm) %*% u_old,
        l1_bound = lambda_w,
        penalty = penalty_w,
        penalty_params = penalty_params_w
      )
      v <- sparse_unit_vector_penalized(
        t(Zme) %*% w_old,
        l1_bound = lambda_v,
        penalty = penalty_v,
        penalty_params = penalty_params_v
      )
      u <- sparse_unit_vector_penalized(
        Zcm %*% w_old,
        l1_bound = lambda_u,
        penalty = penalty_u,
        penalty_params = penalty_params_u
      )
    } else {
      # Standard alternating refinement uses the newest available updates.
      w <- sparse_unit_vector_penalized(
        t(Zcm) %*% u_old,
        l1_bound = lambda_w,
        penalty = penalty_w,
        penalty_params = penalty_params_w
      )
      v <- sparse_unit_vector_penalized(
        t(Zme) %*% w,
        l1_bound = lambda_v,
        penalty = penalty_v,
        penalty_params = penalty_params_v
      )
      u <- sparse_unit_vector_penalized(
        Zcm %*% w,
        l1_bound = lambda_u,
        penalty = penalty_u,
        penalty_params = penalty_params_u
      )
    }

    max_delta <- max(abs(u - u_old), abs(v - v_old), abs(w - w_old))
    if (!is.finite(max_delta)) {
      break
    }

    if (max_delta < tol) {
      converged <- TRUE
      break
    }
  }

  q_cm <- as.numeric(t(u) %*% Zcm %*% w)
  q_me <- as.numeric(t(w) %*% Zme %*% v)

  list(
    u = u,
    w = w,
    v = v,
    q_cm = q_cm,
    q_me = q_me,
    converged = converged,
    iterations = iter
  )
}

jscca_fit <- function(
  C,
  M,
  E,
  ncomp = min(ncol(C), ncol(M), ncol(E)),
  lambda_u,
  lambda_w,
  lambda_v,
  penalty_u = c("l1", "mcp", "adaptive_lasso"),
  penalty_w = c("l1", "mcp", "adaptive_lasso"),
  penalty_v = c("l1", "mcp", "adaptive_lasso"),
  penalty_params_u = list(),
  penalty_params_w = list(),
  penalty_params_v = list(),
  deflation = c("hotelling", "submatrix_uv"),
  update_mode = c("paper", "alternating"),
  standardize = TRUE,
  tol = 1e-6,
  max_iter = 500,
  zero_tol = 1e-10,
  seed = NULL,
  verbose = FALSE
) {
  deflation <- match.arg(deflation)
  update_mode <- match.arg(update_mode)
  penalty_u <- match.arg(penalty_u)
  penalty_w <- match.arg(penalty_w)
  penalty_v <- match.arg(penalty_v)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  C <- as.matrix(C)
  M <- as.matrix(M)
  E <- as.matrix(E)

  if (!(nrow(C) == nrow(M) && nrow(M) == nrow(E))) {
    stop("C, M, E must have the same number of rows (samples).")
  }

  if (standardize) {
    C <- standardize_matrix(C)
    M <- standardize_matrix(M)
    E <- standardize_matrix(E)
  }

  p_full <- ncol(C)
  q_full <- ncol(M)
  r_full <- ncol(E)

  U <- matrix(0, nrow = p_full, ncol = ncomp)
  W <- matrix(0, nrow = q_full, ncol = ncomp)
  V <- matrix(0, nrow = r_full, ncol = ncomp)

  q_cm <- rep(NA_real_, ncomp)
  q_me <- rep(NA_real_, ncomp)
  converged <- rep(FALSE, ncomp)
  iters <- rep(NA_integer_, ncomp)

  selected_u <- vector("list", ncomp)
  selected_w <- vector("list", ncomp)
  selected_v <- vector("list", ncomp)

  idx_c <- seq_len(p_full)
  idx_m <- seq_len(q_full)
  idx_e <- seq_len(r_full)

  Zcm <- crossprod(C[, idx_c, drop = FALSE], M[, idx_m, drop = FALSE])
  Zme <- crossprod(M[, idx_m, drop = FALSE], E[, idx_e, drop = FALSE])

  k_done <- 0

  for (k in seq_len(ncomp)) {
    if (nrow(Zcm) == 0 || ncol(Zcm) == 0 || nrow(Zme) == 0 || ncol(Zme) == 0) {
      if (verbose) message("Stopped: active submatrix became empty at component ", k, ".")
      break
    }

    fit_k <- jscca_component(
      Zcm = Zcm,
      Zme = Zme,
      lambda_u = lambda_u,
      lambda_w = lambda_w,
      lambda_v = lambda_v,
      penalty_u = penalty_u,
      penalty_w = penalty_w,
      penalty_v = penalty_v,
      penalty_params_u = penalty_params_u,
      penalty_params_w = penalty_params_w,
      penalty_params_v = penalty_params_v,
      tol = tol,
      max_iter = max_iter,
      update_mode = update_mode
    )

    u_k <- fit_k$u
    w_k <- fit_k$w
    v_k <- fit_k$v

    u_full <- rep(0, p_full)
    w_full <- rep(0, q_full)
    v_full <- rep(0, r_full)

    u_full[idx_c] <- u_k
    w_full[idx_m] <- w_k
    v_full[idx_e] <- v_k

    U[, k] <- u_full
    W[, k] <- w_full
    V[, k] <- v_full

    q_cm[k] <- fit_k$q_cm
    q_me[k] <- fit_k$q_me
    converged[k] <- fit_k$converged
    iters[k] <- fit_k$iterations

    nz_u_local <- which(abs(u_k) > zero_tol)
    nz_w_local <- which(abs(w_k) > zero_tol)
    nz_v_local <- which(abs(v_k) > zero_tol)

    selected_u[[k]] <- idx_c[nz_u_local]
    selected_w[[k]] <- idx_m[nz_w_local]
    selected_v[[k]] <- idx_e[nz_v_local]

    if (verbose) {
      message(
        sprintf(
          "Component %d: q_cm=%.4f q_me=%.4f nz(u,w,v)=(%d,%d,%d)",
          k, q_cm[k], q_me[k], length(nz_u_local), length(nz_w_local), length(nz_v_local)
        )
      )
    }

    if (deflation == "hotelling") {
      # Dimension-consistent Hotelling-style rank-1 deflation.
      Zcm <- Zcm - abs(fit_k$q_cm) * (u_k %*% t(w_k))
      Zme <- Zme - abs(fit_k$q_me) * (w_k %*% t(v_k))
    } else if (deflation == "submatrix_uv") {
      # Optional deflation: remove variables for non-zero u and non-zero v.
      if (length(nz_u_local) > 0) idx_c <- idx_c[-nz_u_local]
      if (length(nz_v_local) > 0) idx_e <- idx_e[-nz_v_local]

      # M remains unchanged by request; w stays in the original methylation space.
      Zcm <- crossprod(C[, idx_c, drop = FALSE], M[, idx_m, drop = FALSE])
      Zme <- crossprod(M[, idx_m, drop = FALSE], E[, idx_e, drop = FALSE])
    }

    k_done <- k
  }

  if (k_done < ncomp) {
    U <- U[, seq_len(k_done), drop = FALSE]
    W <- W[, seq_len(k_done), drop = FALSE]
    V <- V[, seq_len(k_done), drop = FALSE]
    q_cm <- q_cm[seq_len(k_done)]
    q_me <- q_me[seq_len(k_done)]
    converged <- converged[seq_len(k_done)]
    iters <- iters[seq_len(k_done)]
    selected_u <- selected_u[seq_len(k_done)]
    selected_w <- selected_w[seq_len(k_done)]
    selected_v <- selected_v[seq_len(k_done)]
  }

  list(
    U = U,
    W = W,
    V = V,
    q_cm = q_cm,
    q_me = q_me,
    converged = converged,
    iterations = iters,
    selected_u = selected_u,
    selected_w = selected_w,
    selected_v = selected_v,
    args = list(
      ncomp = ncomp,
      lambda_u = lambda_u,
      lambda_w = lambda_w,
      lambda_v = lambda_v,
      penalty_u = penalty_u,
      penalty_w = penalty_w,
      penalty_v = penalty_v,
      deflation = deflation,
      update_mode = update_mode,
      standardize = standardize,
      tol = tol,
      max_iter = max_iter,
      zero_tol = zero_tol
    )
  )
}

# Example usage:
# set.seed(1)
# n <- 120; p <- 60; q <- 80; r <- 70
# C <- matrix(rnorm(n * p), n, p)
# M <- matrix(rnorm(n * q), n, q)
# E <- matrix(rnorm(n * r), n, r)
# fit <- jscca_fit(C, M, E,
#                  ncomp = 3,
#                  lambda_u = 6,
#                  lambda_w = 8,
#                  lambda_v = 7,
#                  deflation = "submatrix_uv",
#                  update_mode = "paper")
