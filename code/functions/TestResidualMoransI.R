# ==============================================================================
# TestResidualMoransI.R
# Computes inverse-distance spatial weights across plot centroids and performs
# a Monte Carlo permutation test for Moran's I on model residuals.
# ==============================================================================

TestResidualMoransI <- function(
  plots_vect,
  rf_model,
  data_df,
  model_id,
  response_col = "shannon",
  n_perm = 9999,
  seed = 123
) {
  coords_mat <- terra::crds(terra::centroids(plots_vect))
  dist_mat   <- as.matrix(stats::dist(coords_mat))
  inv_dist_w <- 1 / dist_mat
  diag(inv_dist_w) <- 0
  inv_dist_w <- inv_dist_w / rowSums(inv_dist_w)

  obs_vals  <- data_df[[response_col]]
  pred_vals <- stats::predict(rf_model, newdata = data_df)
  resids    <- obs_vals - pred_vals

  n_pts <- length(resids)
  z_res <- resids - mean(resids)
  moran_I_obs <- (n_pts / sum(inv_dist_w)) *
    sum(inv_dist_w * outer(z_res, z_res)) / sum(z_res^2)

  set.seed(seed)
  moran_perm <- replicate(n_perm, {
    z_p <- sample(z_res)
    (n_pts / sum(inv_dist_w)) * sum(inv_dist_w * outer(z_p, z_p)) / sum(z_p^2)
  })

  expected_I <- -1 / (n_pts - 1)
  perm_sd    <- stats::sd(moran_perm)
  z_score    <- (moran_I_obs - mean(moran_perm)) / perm_sd
  p_val_two_sided <- (sum(abs(moran_perm) >= abs(moran_I_obs)) + 1) / (n_perm + 1)

  moran_summary <- tibble::tibble(
    model_id = model_id,
    n_plots = n_pts,
    morans_I_observed = moran_I_obs,
    morans_I_expected = expected_I,
    perm_mean = mean(moran_perm),
    perm_sd = perm_sd,
    z_score = z_score,
    p_value_two_sided = p_val_two_sided,
    n_permutations = n_perm
  )

  residuals_df <- tibble::tibble(
    plot_id = data_df$plot_id,
    model_id = model_id,
    observed = obs_vals,
    predicted = as.numeric(pred_vals),
    residual = as.numeric(resids)
  )

  return(list(
    moran_summary = moran_summary,
    residuals_df = residuals_df
  ))
}
