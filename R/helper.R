.subject_split <- function(n, ratio, subject_id = NULL, validation = FALSE) {
    if (is.null(subject_id)) {
        subject_id <- seq_len(n)
    }
    if (length(subject_id) != n) {
        stop("subject_id must have length equal to the number of rows in x.")
    }
    if (anyNA(subject_id)) {
        stop("subject_id cannot contain missing values.")
    }
    if (!is.numeric(ratio) || anyNA(ratio) || any(ratio < 0)) {
        stop("ratio must be a non-missing numeric vector with non-negative values.")
    }

    subjects <- unique(subject_id)
    n_subjects <- length(subjects)
    subject_order <- sample(subjects, n_subjects, replace = FALSE)

    take_subjects <- function(from, to) {
        if (to < from) {
            return(subject_order[0])
        }
        subject_order[seq.int(from, to)]
    }

    if (validation) {
        if (length(ratio) != 3) {
            stop("ratio must contain training, validation, and testing proportions.")
        }
        if (sum(ratio) <= 0 || sum(ratio) > 1 + sqrt(.Machine$double.eps)) {
            stop("ratio must sum to a value in (0, 1].")
        }
        n_train <- floor(n_subjects * ratio[1])
        n_val <- floor(n_subjects * ratio[2])
        n_test <- n_subjects - n_train - n_val
        if (n_train < 1 || n_test < 1) {
            stop("ratio must allocate at least one subject to training and testing.")
        }

        train_subjects <- take_subjects(1, n_train)
        val_subjects <- take_subjects(n_train + 1, n_train + n_val)
        test_subjects <- take_subjects(n_train + n_val + 1, n_subjects)
    } else {
        if (length(ratio) != 1 || ratio <= 0 || ratio >= 1) {
            stop("ratio must be a single training proportion between 0 and 1.")
        }
        n_train <- floor(n_subjects * ratio)
        n_test <- n_subjects - n_train
        if (n_train < 1 || n_test < 1) {
            stop("ratio must allocate at least one subject to training and testing.")
        }

        train_subjects <- take_subjects(1, n_train)
        val_subjects <- subject_order[0]
        test_subjects <- take_subjects(n_train + 1, n_subjects)
    }

    list(
        train = which(subject_id %in% train_subjects),
        validation = which(subject_id %in% val_subjects),
        test = which(subject_id %in% test_subjects)
    )
}

#' Helper function for cumulative prediction.
#'
#' @importFrom deepTL importDnnet
#' @importFrom deepTL mod_permfit
#' @importFrom deepTL predict_mod_permfit
cv_fit <- function(trainx, trainy, valx, valy, testx, testy,
                    feature_set, mod_args, validate = TRUE,
                    type = c("regression", "binary-classification")[1]) {
    result <- rep(NA, 3)
    if (type == "binary-classification") {
        trainy <- as.factor(trainy)
        if (length(levels(trainy)) != 2) {
            stop("Training data must contain two outcome classes.")
        }
        valy <- factor(valy, levels = levels(trainy))
        testy <- factor(testy, levels = levels(trainy))
        names(result) <- c("acc", "auc", "feature")
    } else {
        names(result) <- c("mse", "pcc", "feature")
    }

    if (validate) {
        val_acc <- rep(NA, length(feature_set))
        k <- 1
        for (j in feature_set) {
            if (length(j) == 0) {
                val_acc[k] <- Inf
                k <- k + 1
                next
            }
            trainset <- deepTL::importDnnet(trainx[, j, drop = FALSE], trainy)
            valset <- deepTL::importDnnet(valx[, j, drop = FALSE], valy)
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
                predy <- ifelse(predy >= .5, levels(trainy)[2], levels(trainy)[1])
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

    selected <- feature_set[[p_mod]]
    if (length(selected) == 0) {
        stop("No features were available for model fitting.")
    }
    trainset <- deepTL::importDnnet(
        rbind(trainx, valx)[, selected, drop = FALSE],
        c(trainy, valy))
    testset <- deepTL::importDnnet(testx[, selected, drop = FALSE], testy)
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
        pred_prob <- as.numeric(predy)
        predy <- ifelse(pred_prob >= .5, levels(trainy)[2], levels(trainy)[1])
        result[1] <- mean(testy == predy)
        result[2] <- suppressMessages(
            as.numeric(pROC::auc(pROC::roc(testy, pred_prob, quiet = TRUE)))
        )
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
    result <- rep(NA, 2)
    if (type == "binary-classification") {
        trainy <- as.factor(trainy)
        if (length(levels(trainy)) != 2) {
            stop("Training data must contain two outcome classes.")
        }
        testy <- factor(testy, levels = levels(trainy))
        names(result) <- c("acc", "auc")
    } else {
        names(result) <- c("mse", "pcc")
    }
    trainset <- deepTL::importDnnet(trainx, trainy)
    testset <- deepTL::importDnnet(testx, testy)

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
        pred_prob <- as.numeric(predy)
        predy <- ifelse(pred_prob >= .5, levels(trainy)[2], levels(trainy)[1])
        result[1] <- mean(testy == predy)
        result[2] <- suppressMessages(
            as.numeric(pROC::auc(pROC::roc(testy, pred_prob, quiet = TRUE)))
        )
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
