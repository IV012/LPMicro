#' Helper function for cumulative prediction.
#'
#' @importFrom deepTL importDnnet
#' @importFrom deepTL mod_permfit
#' @importFrom deepTL predict_mod_permfit
cv_fit <- function(trainx, trainy, valx, valy, testx, testy,
                    feature_set, mod_args, validate = TRUE,
                    type = c("regression", "binary-classification")[1]) {
    result <- rep(NA, 3)
    if (type == "binay-classification") {
        trainy <- as.factor(trainy)
        valy <- as.factor(valy)
        testy <- as.factor(testy)
        names(result) <- c("acc", "auc", "feature")
    } else {
        names(result) <- c("mse", "pcc", "feature")
    }

    if (validate) {
        val_acc <- rep(NA, length(feature_set))
        k <- 1
        for (j in feature_set) {
            trainset <- deepTL::importDnnet(trainx[, j], trainy)
            valset <- deepTL::importDnnet(valx[, j], valy)
            full_mod <- do.call(
                deepTL::mod_permfit,
                c(list(model.type = type, object = trainset), mod_args)
            )
            predy <- deepTL::predict_mod_permfit(
                mod = full_mod,
                object = valset,
                method = mod_args$method,
                model.type = type
            )
            if (type == "binary-classification") {
                predy <- ifelse(predy >= .5, levels(trainy)[1], levels(trainy)[2])
                val_acc[k] <- -mean(valy == predy)
            } else {
                val_acc[k] <- mean((valy - predy)^2)
            }
            k <- k + 1
        }
        p_mod <- tail(which(val_acc == min(val_acc)), 1)
    } else {
        p_mod <- 1
    }

    trainset <- deepTL::importDnnet(
        rbind(trainx, valx)[, feature_set[[p_mod]]],
        c(trainy, valy))
    testset <- deepTL::importDnnet(testx[, feature_set[[p_mod]]], testy)
    full_mod <- do.call(
        deepTL::mod_permfit,
        c(list(model.type = type, object = trainset), mod_args)
    )
    predy <- deepTL::predict_mod_permfit(
        mod = full_mod,
        object = testset,
        method = mod_args$method,
        model.type = type
    )

    if (type == "binary-classification") {
        predy <- ifelse(predy >= .5, levels(trainy)[1], levels(trainy)[2])
        result[1] <- mean(testy == predy)
        result[2] <- suppressMessages(pROC::auc(pROC::roc(testy, predy)))
    } else {
        result[1] <- mean((testy - predy)^2)
        result[2] <- cor(testy, predy, method = "pearson")
    }
    result[3] <- p_mod

    return(result)
}
#' Helper function for visitwise prediction.
#'
#' @importFrom deepTL importDnnet
#' @importFrom deepTL mod_permfit
#' @importFrom deepTL predict_mod_permfit
visit_fit <- function(trainx, trainy, testx, testy, mod_args,
                    type = c("regression", "binary-classification")[1]) {
    trainset <- deepTL::importDnnet(trainx, trainy)
    testset <- deepTL::importDnnet(testx, testy)
    result <- rep(NA, 2)
    if (type == "binay-classification") {
        trainy <- as.factor(trainy)
        valy <- as.factor(valy)
        testy <- as.factor(testy)
        names(result) <- c("acc", "auc", "feature")
    } else {
        names(result) <- c("mse", "pcc", "feature")
    }

    full_mod <- do.call(
        deepTL::mod_permfit,
        c(list(model.type = type, object = trainset), mod_args)
    )
    predy <- deepTL::predict_mod_permfit(
        mod = full_mod,
        object = testset,
        method = mod_args$method,
        model.type = type
    )

    if (type == "binary-classification") {
        predy <- ifelse(predy >= .5, levels(trainy)[1], levels(trainy)[2])
        result[1] <- mean(testy == predy)
        result[2] <- suppressMessages(pROC::auc(pROC::roc(testy, predy)))
    } else {
        result[1] <- mean((testy - predy)^2)
        result[2] <- cor(testy, predy, method = "pearson")
    }
    return(result)
}

#' Sample Data Generator
#'
#' @param seed Random seed. Default 1.
#' @param n The number of subjects. Default 200.
#' @param v The number of visits. Default 5.
#' @param p The number of microbial taxa. Default 100.
#'
#' @return Return a list of z (the design matrix) and y (outcomes).
#'
#' @importFrom expm sqrtm
#' @importFrom Matrix bdiag
#'
#' @export
generate_data <- function(seed = 1, n = 200, v = 5, p = c(100, 200, 500)[1]) {
  set.seed(seed)
  lambda <- matrix(rnorm(n * p * v), n, p * v)
  mat_sigma <- expm::sqrtm(diag(0.9, p, p) + 0.1)
  mat_sigma <- Matrix::bdiag(replicate(v, mat_sigma, simplify = FALSE))
  lambda <- lambda %*% mat_sigma
  x <- matrix(rpois(n * p * v, as.vector(exp(lambda))), n, p * v)
  z <- t(t(x) / rowSums(x))
  x <- apply(z, 2, scale)
  y <- 0
  for (i in c(1)){
    active <- (0:(v - 1)) * p + i
    xsub <- x[, active]
    beta <- runif(v, -2, 2)
    y <- y + cos(xsub) %*% beta
  }
  for (i in c(2:5)){
    active <- c((v - 2) * p + i, (v - 1) * p + i)
    xsub <- x[, active]
    gamma <- runif(2, 0, 4)
    y <- y + log((xsub^2) %*% gamma + 1)
  }
  return(list(
    z = z,
    y = scale(y) + rnorm(n, 0, 0.1)
  ))
}