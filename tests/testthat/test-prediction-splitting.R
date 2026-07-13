source_lpmicro <- function() {
  env <- new.env(parent = globalenv())
  root <- normalizePath(file.path(testthat::test_path("..", "..")),
                        mustWork = TRUE)
  sys.source(file.path(root, "R", "helper.R"), envir = env)
  sys.source(file.path(root, "R", "LPMicro.R"), envir = env)
  env
}

test_that("screen_taxa uses one polynomial block per requested order", {
  env <- source_lpmicro()
  set.seed(1)
  x <- matrix(rnorm(30 * 6), 30, 6)
  y <- rnorm(30)
  idx <- rep(seq_len(3), 2)

  expect_error(
    env$screen_taxa(x, y, idx, lambda = 0.01),
    NA
  )
  expect_error(
    env$screen_taxa(x, y, idx, ord = 1, lambda = 0.01),
    NA
  )
})

test_that("subject split keeps repeated rows from the same subject together", {
  env <- source_lpmicro()
  subject_id <- rep(seq_len(12), each = 2)

  set.seed(10)
  split <- env$.subject_split(
    n = length(subject_id),
    ratio = c(0.5, 0.25, 0.25),
    subject_id = subject_id,
    validation = TRUE
  )

  train_subjects <- unique(subject_id[split$train])
  validation_subjects <- unique(subject_id[split$validation])
  test_subjects <- unique(subject_id[split$test])

  expect_length(intersect(train_subjects, validation_subjects), 0)
  expect_length(intersect(train_subjects, test_subjects), 0)
  expect_length(intersect(validation_subjects, test_subjects), 0)
  expect_equal(length(split$train), length(train_subjects) * 2)
  expect_equal(length(split$validation), length(validation_subjects) * 2)
  expect_equal(length(split$test), length(test_subjects) * 2)
})

test_that("cumulative prediction reuses one subject split and screens on training data", {
  env <- source_lpmicro()
  sim <- env$generate_data(seed = 2, n = 30, v = 5, p = 6)
  x <- sim$z
  y <- as.numeric(sim$y)
  p <- 6
  v <- 5
  idx <- rep(seq_len(p), v)
  taxa_list <- lapply(seq_len(v), function(tp) {
    cols <- seq_len(p * tp)
    list(x = x[, cols, drop = FALSE], y = y, idx = idx[cols])
  })

  env$.screen_calls <- list()
  env$.cv_calls <- list()
  env$screen_taxa <- function(x, y, idx, ...) {
    env$.screen_calls[[length(env$.screen_calls) + 1]] <<- list(
      n = nrow(x),
      y = y,
      idx = idx
    )
    list(seq_len(min(2, ncol(x))))
  }
  env$cv_fit <- function(trainx, trainy, valx, valy, testx, testy,
                         feature_set, mod_args, validate = TRUE,
                         type = "regression") {
    expect_true(is.list(feature_set))
    expect_true(is.logical(validate))
    env$.cv_calls[[length(env$.cv_calls) + 1]] <<- list(
      trainy = trainy,
      valy = valy,
      testy = testy,
      validate = validate
    )
    c(mse = mean(testy), pcc = length(testy), feature = 1)
  }

  set.seed(20)
  result <- env$cumulative_predict(
    taxa_list,
    mod_args = list(method = "mock"),
    ratio = c(0.6, 0.2, 0.2),
    plot.result = FALSE
  )

  expect_equal(nrow(result), v)
  expect_length(env$.screen_calls, v)
  expect_length(env$.cv_calls, v)
  expect_equal(env$.screen_calls[[1]]$n, length(env$.cv_calls[[1]]$trainy))
  for (i in seq_len(v)) {
    expect_identical(env$.cv_calls[[i]]$trainy, env$.cv_calls[[1]]$trainy)
    expect_identical(env$.cv_calls[[i]]$valy, env$.cv_calls[[1]]$valy)
    expect_identical(env$.cv_calls[[i]]$testy, env$.cv_calls[[1]]$testy)
  }
})
