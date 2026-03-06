# Run Instructions for JSCCA

This project contains a standalone JSCCA implementation in:

- `jscca.R`

and a runnable vignette in:

- `vignettes/jscca_vignette.Rmd`

## 1) Prerequisites

- R (>= 4.1 recommended)
- Optional for vignette rendering: `rmarkdown`
- Pandoc (required by `rmarkdown::render`)

Install `rmarkdown` if needed:

```r
install.packages("rmarkdown")
```

Check pandoc availability:

```r
rmarkdown::pandoc_available()
```

If `FALSE`, install Pandoc or render from RStudio (which bundles Pandoc).

## 2) Quick Run (R Console)

From the project root:

```r
source("jscca.R")

set.seed(1)
n <- 120; p <- 50; q <- 70; r <- 60
C <- matrix(rnorm(n * p), n, p)
M <- matrix(rnorm(n * q), n, q)
E <- matrix(rnorm(n * r), n, r)

fit <- jscca_fit(
  C, M, E,
  ncomp = 3,
  lambda_u = 6,
  lambda_w = 8,
  lambda_v = 7,
  penalty_u = "mcp",
  penalty_w = "l1",
  penalty_v = "mcp",
  penalty_params_u = list(gamma = 3),
  penalty_params_v = list(gamma = 3),
  deflation = "submatrix_uv",  # or "hotelling"
  update_mode = "paper",
  standardize = TRUE,
  verbose = TRUE
)

fit$q_cm
fit$q_me
fit$selected_u
fit$selected_w
fit$selected_v
```

## 3) Run via Rscript (non-interactive)

```bash
Rscript -e "source('jscca.R'); set.seed(1); n<-100; p<-40; q<-50; r<-45; C<-matrix(rnorm(n*p),n,p); M<-matrix(rnorm(n*q),n,q); E<-matrix(rnorm(n*r),n,r); fit<-jscca_fit(C,M,E,ncomp=2,lambda_u=5,lambda_w=6,lambda_v=5,deflation='submatrix_uv'); print(fit$q_cm); print(fit$q_me)"
```

## 4) Render the Vignette

From project root:

```bash
Rscript -e "rmarkdown::render('vignettes/jscca_vignette.Rmd', output_file='jscca_vignette.html')"
```

This creates `vignettes/jscca_vignette.html`.

## 5) Main Arguments

- `ncomp`: maximum number of components.
- `lambda_u`, `lambda_w`, `lambda_v`: sparsity bounds.
- `penalty_u`, `penalty_w`, `penalty_v`: penalty type for each loading vector.
  - Options: `"l1"`, `"mcp"`, `"adaptive_lasso"`.
- `penalty_params_u`, `penalty_params_w`, `penalty_params_v`: per-penalty settings.
  - For `"mcp"`: use `list(gamma = 3)` (must have `gamma > 1`).
  - For `"adaptive_lasso"`: use `list(gamma = 1.25, eps = 1e-3)`.
- `deflation`:
  - `"hotelling"`: rank-1 deflation of cross-product matrices.
  - `"submatrix_uv"`: remove variables corresponding to non-zero `u` and non-zero `v` after each component.
- `update_mode`:
  - `"paper"`: appendix-style previous-iterate updates.
  - `"alternating"`: newest-available updates each cycle.
- `standardize`: if `TRUE`, z-score columns before fitting.

## 6) Output Structure

`jscca_fit()` returns a list with:

- `U`, `W`, `V`: loading matrices (full original feature dimensions).
- `q_cm`, `q_me`: canonical correlations per component.
- `selected_u`, `selected_w`, `selected_v`: selected feature indices per component.
- `converged`, `iterations`: optimization diagnostics.
