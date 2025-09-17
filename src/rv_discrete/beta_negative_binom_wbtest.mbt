///|
test "BetaNegativeBinom::pmf" {
  let dist = BetaNegativeBinom::new(5, 2.0, 3.0, loc=0.0)
  // Test valid integer values - should return positive probabilities
  let pmf_0 = dist.pmf(0.0)
  let pmf_1 = dist.pmf(1.0)
  let pmf_5 = dist.pmf(5.0)
  assert_true((pmf_0 - 0.047619047619047616).abs() < 1.0e-10)
  assert_true((pmf_1 - 0.07142857142857144).abs() < 1.0e-10)
  assert_true((pmf_5 - 0.06293706293706297).abs() < 1.0e-10)
  let dist2 = BetaNegativeBinom::new(5, 2.0, 3.0)
  let pmf_200 = dist2.pmf(200)
  assert_true((pmf_200 - 3.846350701402042e-5).abs() < 1.0e-10)
}

///|
test "BetaNegativeBinom::pmf/boundary_cases" {
  let dist = BetaNegativeBinom::new(3, 1.5, 2.5, loc=1.0)
  // Test non-integer values - should return 0
  inspect(dist.pmf(1.5), content="0")
  inspect(dist.pmf(2.7), content="0")
  // Test negative values (below loc) - should return 0
  inspect(dist.pmf(0.0), content="0")
  inspect(dist.pmf(-1.0), content="0")
  // Test exactly at loc (0 failures) - should be positive
  let pmf_loc = dist.pmf(1.0)
  assert_true((pmf_loc - 0.10937500000000001).abs() < 1.0e-10)
}

///|
test "BetaNegativeBinom::pmf/location_parameter" {
  let dist_loc5 = BetaNegativeBinom::new(2, 1.0, 1.0, loc=5.0)
  // With loc=5.0, k=5.0 represents 0 failures, k=6.0 represents 1 failure
  let pmf_zero_failures = dist_loc5.pmf(5.0)
  let pmf_one_failure = dist_loc5.pmf(6.0)
  // Both should be positive for valid failure counts
  assert_true((pmf_zero_failures - 0.3333333333333333).abs() < 1.0e-10)
  assert_true((pmf_one_failure - 0.16666666666666663).abs() < 1.0e-10)
  // Values below loc should return 0
  inspect(dist_loc5.pmf(4.0), content="0")
  inspect(dist_loc5.pmf(4.9), content="0")
}

///|
test "BetaNegativeBinom::cdf" {
  let dist = BetaNegativeBinom::new(5, 2.0, 3.0, loc=0.0)

  // Test that CDF is monotonically increasing
  let cdf_0 = dist.cdf(0.0)
  let cdf_5 = dist.cdf(5.0)
  let cdf_10 = dist.cdf(10.0)
  assert_true((cdf_0 - dist.pmf(0.0)).abs() < 0.001)
  assert_true(cdf_5 >= cdf_0)
  assert_true(cdf_10 >= cdf_5)

  // Test that CDF approaches 1 for large values
  let cdf_large = dist.cdf(100.0)
  assert_true(cdf_large > 0.9)
}

///|
test "BetaNegativeBinom::cdf/boundary_cases" {
  let dist = BetaNegativeBinom::new(3, 1.5, 2.5, loc=2.0)

  // Test x below location parameter
  let cdf_below = dist.cdf(1.0)
  assert_true((cdf_below - 0.0).abs() < 0.001)

  // Test x exactly at location parameter
  let cdf_at_loc = dist.cdf(2.0)
  assert_true((cdf_at_loc - dist.pmf(2.0)).abs() < 0.001)

  // Test fractional x values
  let cdf_frac = dist.cdf(2.7)
  assert_true((cdf_frac - dist.cdf(2.0)).abs() < 0.001)
}

///|
test "BetaNegativeBinom::cdf/properties" {
  let dist = BetaNegativeBinom::new(2, 3.0, 4.0, loc=1.0)

  // Test that CDF at x equals sum of PMF up to floor(x - loc)
  let x = 5.5
  let expected = dist.pmf(1.0) +
    dist.pmf(2.0) +
    dist.pmf(3.0) +
    dist.pmf(4.0) +
    dist.pmf(5.0)
  let actual = dist.cdf(x)
  assert_true((actual - expected).abs() < 0.001)

  // Test CDF is non-decreasing
  assert_true(dist.cdf(3.0) <= dist.cdf(4.0))
  assert_true(dist.cdf(4.0) <= dist.cdf(5.0))
}

///|
test "BetaNegativeBinom::ppf" {
  let dist = BetaNegativeBinom::new(5, 2.0, 3.0, loc=10.0)

  // Test q = 0.0 returns loc - 1.0
  inspect(dist.ppf(0.0), content="9")

  // Test normal cases - ppf should return values where cdf >= q
  let q1 = 0.1
  let result1 = dist.ppf(q1)
  assert_true((11.0 - result1).abs() < 1.0e-10)
  let q2 = 0.5
  let result2 = dist.ppf(q2)
  assert_true((17.0 - result2).abs() < 1.0e-10)
}

///|
test "panic BetaNegativeBinom::ppf/invalid_q_negative" {
  let dist = BetaNegativeBinom::new(5, 2.0, 3.0)
  ignore(dist.ppf(-0.1))
}

///|
test "panic BetaNegativeBinom::ppf/invalid_q_one" {
  let dist = BetaNegativeBinom::new(5, 2.0, 3.0)
  ignore(dist.ppf(1.0))
}

///|
test "BetaNegativeBinom::mean" {
  let bnb1 = BetaNegativeBinom::new(5, 2.0, 3.0, loc=0.0)
  assert_true((bnb1.mean() - 15.0).abs() < 1.0e-15)
  let bnb2 = BetaNegativeBinom::new(10, 3.5, 2.5, loc=5.0)
  assert_true((bnb2.mean() - 15.0).abs() < 1.0e-15)
  let bnb3 = BetaNegativeBinom::new(1, 5.0, 1.0, loc=2.0)
  assert_true((bnb3.mean() - 2.25).abs() < 1.0e-15)
}

///|
test "BetaNegativeBinom::varc" {
  let dist1 = BetaNegativeBinom::new(5, 3.0, 2.0, loc=10.0)
  let expected1 = 70.0
  assert_true((dist1.varc() - expected1).abs() < 1.0e-15)
  let dist2 = BetaNegativeBinom::new(10, 5.0, 3.0, loc=0.0)
  let expected2 = 61.25
  assert_true((dist2.varc() - expected2).abs() < 1.0e-15)
}

///|
test "BetaNegativeBinom::entropy" {
  let dist1 = BetaNegativeBinom::new(5, 3.0, 2.0, loc=10.0)
  inspect(dist1.entropy(), content="2.66140988208517")
  let dist2 = BetaNegativeBinom::new(10, 5.0, 3.0, loc=0.0)
  inspect(dist2.entropy(), content="3.034913867925415")
}

///|
test "BetaNegativeBinom::interval" {
  let dist = BetaNegativeBinom::new(3, 4.0, 5.0, loc=2.0)
  let alpha = 0.1
  let (lower, upper) = dist.interval(alpha)
  assert_true((lower - 2.0).abs() < 1.0e-10)
  assert_true((upper - 18.0).abs() < 1.0e-10)
}

///|
test "BetaNegativeBinom::skew" {
  let dist = BetaNegativeBinom::new(3, 12.0, 8.0)
  let expected = 2.0184766831204053
  assert_true((dist.skew() - expected).abs() < 1.0e-10)
}

///|
test "BetaNegativeBinom::kurt" {
  let dist = BetaNegativeBinom::new(3, 12.0, 8.0)
  let expected = 7.494987468671679
  assert_true((dist.kurt() - expected).abs() < 1.0e-10)
}
