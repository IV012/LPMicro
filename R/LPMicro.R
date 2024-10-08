#' Microbial Taxa Screening using Polynomial Group Lasso
#'
#' @param x A n by p * q matrix, where each row denotes a subject,
#' and each column denotes a microbial taxon (p) at a single timepoint (q).
#' @param y An outcome vector of n subjects.
#' @param idx A numerical index vector with its length equaling to
#' the number of columns in x. The idx should range from 1 to p, indicating
#' the microbial taxa of the columns in x.
#' @param ord The order of splines. Default 3.
#' @param spline  The spline type. Default natural splines.
#' @param lambda  The penalty parameter. Default NULL.
#' @param intercept Include the intercept or not. Default FALSE.
#'
#' @return Return the index of screened features.
#'
#' @importFrom grplasso grplasso
#' @importFrom grplasso lambdamax
#' @importFrom grplasso LinReg
#' @importFrom grplasso LogReg
#'
#' @export
screen_taxa <- function(x, y, idx, ord = 3,
                           spline = c("natural", "b-spline")[1],
                           lambda = NULL, intercept = FALSE) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  y <- as.numeric(y)
  if (n != length(y)) {
    stop("The number of observations does not match in x and y.")
  }
  if (p != length(idx)) {
    stop("The number of variables does not match in x and idx.")
  }

  if (spline == "natural") {
    for (m in 1:ord){
      x <- cbind(x, x^m)
    }
  } else {
    stop("Spline not supported. Please wait for future updates.")
  }

  if (intercept) {
    idx <- c(0, rep(idx, ord))
  } else {
    idx <- c(NA, rep(idx, ord))
  }
  if (length(unique(y)) == 2) {
    model <- grplasso::LogReg()
  } else {
    model <- grplasso::LinReg()
  }

  if (is.null(lambda)) {
    lambda <- grplasso::lambdamax(x, y = y, index = idx, penscale = sqrt,
                        model = model) * 0.5
  }
  grp_fit <- grplasso::grplasso(x, y, index = idx, model = model,
                          penscale = sqrt, lambda = lambda,
                          control = grplasso::grpl.control(
                            update.hess = "lambda", trace = 0))
  coef.lasso <- grp_fit$coefficients
  feature <- which(coef.lasso != 0) - 1

  return(feature)
}

#' Cumulative Prediction
#'
#' @param taxa_list A list of taxa. Each element denotes the information
#' up to a chosen visit in a list of x (n by p matrix of predictors), y (
#' the outcome vector of length n), and idx (the variable annotation
#' vector of length p, identifying the microbial taxa).
#' @param mod_args The model and its parameters to be used for prediction.
#' @param screen To screen using group lasso or not. Default TRUE.
#' @param ratio A vector of the proportion of training, validating, and testing
#' data. Default c(0.7, 0.1, 0.2). Note that validating set can be set to 0.
#' @param plot.result Whether to plot the result. Default TRUE.
#' @param type Type of prediction, depending on continuous outcomes or
#' binary outcomes. Default "regression".
#' @param ... Other arguments to be passed to screen_taxa.
#'
#' @return Return a dataframe, whether the first two columns record the
#' prediction accuracy (MSE, PCC for regression, ACC and AUC for
#' classification.)
#'
#' @export
cumulative_predict <- function(
      taxa_list, mod_args, screen = TRUE, ratio = c(0.7, 0.1, 0.2),
      plot.result = TRUE,
      type = c("regression", "binary-classification")[1], ...) {

  n_tp <- length(taxa_list)
  result <- matrix(NA, n_tp, 2)
  result <- as.data.frame(result)
  if (type == "binary-classification") {
    colnames(result) <- c("acc", "auc")
  } else {
    colnames(result) <- c("mse", "pcc")
  }

  i <- 1
  for (df in taxa_list) {
    x <- df$x
    y <- df$y
    idx <- df$idx
    n <- dim(x)[1]

    tr_idx <- sample(seq(n), n * ratio[1], replace = FALSE)
    val_idx <- sample(seq(n)[-tr_idx], n * ratio[2], replace = FALSE)
    trainx <- x[tr_idx, ]
    trainy <- y[tr_idx]
    valx <- x[val_idx, ]
    valy <- y[valy]
    testx <- x[-c(tr_idx, val_idx), ]
    testy <- y[-c(tr_idx, val_idx)]
    validate <- ifelse(length(val_idx) >= 3, TRUE, FALSE)

    if (screen) {
      feature_set <- screen_taxa(x, y, idx, ...)
    }else {
      feature_set <- seq(dim(x)[2])
    }

    result[i, ] <- cv_fit(trainx, trainy, valx, valy, testx, testy,
                          validate, feature_set, mod_args, type)[1:2]
    i <- i + 1
  }

  return(result)
}

#' Visitwise Prediction
#'
#' @param taxa_list A list of taxa. Each element denotes the information
#' from one visit in a list of x (n by p matrix of predictors) and y (
#' the outcome vector of length n).
#' @param mod_args The model and its parameters to be used for prediction.
#' @param plot.result Whether to plot the result. Default TRUE.
#' @param ratio The proportion of training data. Default 0.8.
#' @param type Type of prediction, depending on continuous outcomes or
#' binary outcomes. Default "regression".
#'
#' @return Return a dataframe, whether the first two columns record the
#' prediction accuracy (MSE, PCC for regression, ACC and AUC for
#' classification.)
#'
#' @export
visit_predict <-  function(
      taxa_list, mod_args, plot.result = TRUE, ratio = 0.8,
      type = c("regression", "binary-classification")[1]) {

  n_tp <- length(taxa_list)
  result <- matrix(NA, n_tp, 2)
  result <- as.data.frame(result)
  if (type == "binary-classification") {
    colnames(result) <- c("acc", "auc")
  } else {
    colnames(result) <- c("mse", "pcc")
  }

  i <- 1
  for (df in taxa_list) {
    x <- df$x
    y <- df$y
    n <- dim(x)[1]

    tr_idx <- sample(seq(n), n * ratio, replace = FALSE)
    trainx <- x[tr_idx, ]
    trainy <- y[tr_idx]
    testx <- x[-tr_idx, ]
    testy <- y[-tr_idx]

    result[i, ] <- visit_fit(trainx, trainy, testx, testy, mod_args, type)

    i <- i + 1
  }

  return(result)
}

#' Microbial Taxa Interpretation using Permutation Importance Score
#'
#' @param x A n by p * q matrix, where each row denotes a subject,
#' and each column denotes a microbial taxon (p) at a single timepoint (q).
#' @param y An outcome vector of n subjects.
#' @param idx_tx A numerical index vector with its length equaling to
#' the number of columns in x. The idx should range from 1 to p, indicating
#' the microbial taxa of the columns in x.
#' @param idx_tp A numerical index vector with its length equaling to
#' the number of columns in x. The idx should range from 1 to q, indicating
#' the time point of the columns in x.
#' @param mod_args The model and its parameters to be used for prediction.
#' @param k_fold Number of folds for permutation importance score. Default 5.
#' @param screen Whether to screen the taxa. Default TRUE.
#' @param ... Other arguments to be passed to screen_taxa.
#'
#' @return Return a list of three dataframes: (i) grid: grid importance, (ii)
#' taxa: taxa importance, and (iii) visits: visits importance.
#'
#' @importFrom deepTL permfit
#'
#' @export
taxa_interpret <- function(x, y, idx_tx, idx_tp, mod_args,
     k_fold = 5, screen = TRUE, ...) {

  if (screen) {
    feature <- screen_taxa(x, y, idx_tx, ...)
  } else {
    feature <- seq(dim(x)[2])
  }
  xs <- x[, feature]
  idx_tx <- idx_tx[feature]
  idx_tp <- idx_tp[feature]
  trainset <- deepTL::importDnnet(xs, y)

  pathway <- list()
  # create taxa groups
  n_txs <- 0
  for (tx in unique(idx_tx)) {
    tx_name <- paste0("Taxa:", tx)
    if (!(tx_name %in% names(pathway))) {
      if (sum(idx_tx == tx) > 1) {
        pathway[[tx_name]] <- which(idx_tx == tx)
        n_txs <- n_txs + 1
      }
    }
  }
  # create visit groups
  n_tps <- 0
  for (tp in unique(idx_tp)) {
    tp_name <- paste0("Visit:", tp)
    if (!(tp_name) %in% names(pathway)) {
      if (sum(idx_tp == tp)) {
        pathway[[tp_name]] <- which(idx_tp == tp)
        n_tps <- n_tps + 1
      }
    }
  }
  # exclude groups with only one or less member
  pathway <- pathway[sapply(pathway, length) > 1]
  taxa <- data.frame(p.val = rep(NA, n_txs))
  visits <- data.frame(p.val = rep(NA, n_tps))

  nshuffle <- sample(length(trainset@y))
  permod <- do.call(
    deepTL::permfit,
    c(
      list(
        train = trainset,
        k_fold = k_fold,
        n_perm = 100,
        shuffle = nshuffle,
        pathway_list = pathway
      ),
      mod_args
    )
  )
  value <- permod@importance@importance_pval_x
  taxa[, 1] <- permod@block_importance$importance_pval_x[seq(n_txs)]
  visits[, 1] <- permod@block_importance$importance_pval_x[seq(n_txs+1, n_txs+n_tps)]

  # grid importance
  value[is.na(value)] <- 1
  grid <- data.frame(p.val = value, method = mod_args$method)
  grid$logp <- -log(grid$p.val, 10)
  grid$feature <- colnames(xs)
  grid$timepoint <- idx_tp
  grid$taxa <- idx_tx

  # visit importance
  visits[is.na(visits)] <- 1
  visits$logp <- -log(visits[, 1], 10)
  visits$timepoint <- names(pathway)[seq(n_txs + 1, n_txs + n_tps)]

  # taxa importance
  taxa[is.na(taxa)] <- 1
  taxa$logp <- -log(taxa[, 1], 10)
  taxa$taxa <- names(pathway)[seq(n_txs)]
  return(list(grid = grid, visits = visits, taxa = taxa))
}