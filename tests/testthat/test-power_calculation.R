test_that("powerCalculation returns expected output structure", {
  result <- powerCalculation(data.type = "gut", method = "amp", nSim = 10, nSam = 50,nOTU = 200)

  # Check that the result is a list
  expect_type(result, "list")

  # Check that the list contains the expected elements
  expect_named(result, c("overall_power", "replicate_power"))

  # Check that overall_power is a numeric value between 0 and 1
  expect_type(result$overall_power, "double")
  expect_true(result$overall_power >= 0 && result$overall_power <= 1)

  # Check that replicate_power is a numeric vector of length nSim
  expect_type(result$replicate_power, "double")
  expect_length(result$replicate_power, 10)
})

test_that("powerCalculation handles invalid nSam values appropriately", {
  # Expect an error when nSam exceeds the maximum allowed for 'gut' data type
  expect_error(powerCalculation(data.type = "gut", method = "amp", nSam = 300),
               "For 'gut' data, the maximum sample size is 206. Please choose nSam <= 206.")
})

test_that("powerCalculation produces consistent results with fixed random seed", {
  set.seed(123)
  result1 <- powerCalculation(data.type = "gut", method = "amp", nSim = 10, nSam = 50)

  set.seed(123)
  result2 <- powerCalculation(data.type = "gut", method = "amp", nSim = 10, nSam = 50)

  # Check that the results are identical
  expect_equal(result1, result2)
})

