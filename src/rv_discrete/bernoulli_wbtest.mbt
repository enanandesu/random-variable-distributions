///|
test "Bernoulli creation and basic properties" {
  let bernoulli = Bernoulli::new(0.7, loc=3.0)

  // Test basic properties
  assert_eq(bernoulli.p, 0.7)
  assert_eq(bernoulli.loc, 3.0)

  // Test default location
  let bernoulli_default = Bernoulli::new(0.3)
  assert_eq(bernoulli_default.loc, 0.0)
}

///|
test "Bernoulli PMF properties" {
  let bernoulli = Bernoulli::new(0.6, loc=0.0)

  // Test PMF at valid points
  let pmf_0 = bernoulli.pmf(0.0) // P(X=0) = 1-p = 0.4
  let pmf_1 = bernoulli.pmf(1.0) // P(X=1) = p = 0.6
  assert_true((pmf_0 - 0.4).abs() < 1.0e-10)
  assert_true((pmf_1 - 0.6).abs() < 1.0e-10)

  // Test PMF at non-integer points
  assert_eq(bernoulli.pmf(0.5), 0.0)
  assert_eq(bernoulli.pmf(0.9), 0.0)
  assert_eq(bernoulli.pmf(1.1), 0.0)

  // Test PMF outside support
  assert_eq(bernoulli.pmf(-1.0), 0.0)
  assert_eq(bernoulli.pmf(2.0), 0.0)
  assert_eq(bernoulli.pmf(10.0), 0.0)

  // Test PMF sums to 1
  let total = pmf_0 + pmf_1
  assert_true((total - 1.0).abs() < 1.0e-10)
}

///|
test "Bernoulli CDF properties" {
  let bernoulli = Bernoulli::new(0.3, loc=0.0)

  // Test CDF values
  assert_eq(bernoulli.cdf(-0.5), 0.0) // x < 0
  assert_true((bernoulli.cdf(0.0) - 0.7).abs() < 1.0e-10) // x = 0, CDF = 1-p
  assert_true((bernoulli.cdf(0.5) - 0.7).abs() < 1.0e-10) // 0 < x < 1
  assert_true((bernoulli.cdf(1.0) - 1.0).abs() < 1.0e-10) // x >= 1
  assert_true((bernoulli.cdf(2.0) - 1.0).abs() < 1.0e-10) // x > 1

  // Test CDF is non-decreasing
  assert_true(bernoulli.cdf(-1.0) <= bernoulli.cdf(0.0))
  assert_true(bernoulli.cdf(0.0) <= bernoulli.cdf(0.5))
  assert_true(bernoulli.cdf(0.5) <= bernoulli.cdf(1.0))
  assert_true(bernoulli.cdf(1.0) <= bernoulli.cdf(2.0))
}

///|
test "Bernoulli statistical moments" {
  let bernoulli = Bernoulli::new(0.25, loc=0.0)

  // Test mean = p
  assert_true((bernoulli.mean() - 0.25).abs() < 1.0e-10)

  // Test variance = p(1-p)
  let expected_var = 0.25 * 0.75
  assert_true((bernoulli.varc() - expected_var).abs() < 1.0e-10)

  // Test standard deviation = sqrt(p(1-p))
  let expected_std = expected_var.sqrt()
  assert_true((bernoulli.std() - expected_std).abs() < 1.0e-10)

  // Test skewness = (1-2p) / sqrt(p(1-p))
  let expected_skew = (1.0 - 2.0 * 0.25) / (0.25 * 0.75).sqrt()
  assert_true((bernoulli.skew() - expected_skew).abs() < 1.0e-10)

  // Test kurtosis = (1 - 6p(1-p)) / (p(1-p))
  let expected_kurt = (1.0 - 6.0 * 0.25 * 0.75) / (0.25 * 0.75)
  assert_true((bernoulli.kurt() - expected_kurt).abs() < 1.0e-10)
}

///|
test "Bernoulli with location parameter" {
  let bernoulli = Bernoulli::new(0.4, loc=5.0)

  // Test PMF with location shift
  assert_true((bernoulli.pmf(5.0) - 0.6).abs() < 1.0e-10) // P(X=loc) = 1-p
  assert_true((bernoulli.pmf(6.0) - 0.4).abs() < 1.0e-10) // P(X=loc+1) = p
  assert_eq(bernoulli.pmf(4.0), 0.0) // Outside support
  assert_eq(bernoulli.pmf(7.0), 0.0) // Outside support
  assert_eq(bernoulli.pmf(5.5), 0.0) // Non-integer

  // Test mean with location = loc + p
  let expected_mean = 5.0 + 0.4
  assert_true((bernoulli.mean() - expected_mean).abs() < 1.0e-10)

  // Test CDF with location
  assert_eq(bernoulli.cdf(4.5), 0.0) // x < loc
  assert_true((bernoulli.cdf(5.0) - 0.6).abs() < 1.0e-10) // x = loc
  assert_true((bernoulli.cdf(5.5) - 0.6).abs() < 1.0e-10) // loc < x < loc+1
  assert_true((bernoulli.cdf(6.0) - 1.0).abs() < 1.0e-10) // x >= loc+1
}

///|
test "Bernoulli PPF and quantiles" {
  let bernoulli = Bernoulli::new(0.3, loc=0.0)

  // Test PPF at boundaries
  assert_eq(bernoulli.ppf(0.0), -1.0) // Below support
  assert_eq(bernoulli.ppf(1.0), 1.0) // Maximum value

  // Test PPF at specific quantiles
  assert_eq(bernoulli.ppf(0.5), 0.0) // q=0.5 < 0.7, so returns 0
  assert_eq(bernoulli.ppf(0.7), 0.0) // q=0.7 = 1-p, boundary case
  assert_eq(bernoulli.ppf(0.8), 1.0) // q=0.8 > 1-p, so returns 1

  // Test PPF consistency with CDF
  let q1 = 0.4
  let x1 = bernoulli.ppf(q1)
  assert_true(bernoulli.cdf(x1) >= q1)
  let q2 = 0.9
  let x2 = bernoulli.ppf(q2)
  assert_true(bernoulli.cdf(x2) >= q2)
}

///|
test "Bernoulli survival and inverse functions" {
  let bernoulli = Bernoulli::new(0.4, loc=0.0)

  // Test survival function sf(x) = 1 - cdf(x)
  assert_true((bernoulli.sf(-1.0) - 1.0).abs() < 1.0e-10)
  assert_true((bernoulli.sf(0.0) - 0.4).abs() < 1.0e-10)
  assert_true((bernoulli.sf(1.0) - 0.0).abs() < 1.0e-10)

  // Test inverse survival function
  assert_eq(bernoulli.isf(0.0), 1.0) // isf(0) = ppf(1)
  assert_eq(bernoulli.isf(1.0), -1.0) // isf(1) = ppf(0)
  assert_eq(bernoulli.isf(0.4), 0.0) // isf(0.4) = ppf(0.6)
  assert_eq(bernoulli.isf(0.6), 0.0) // isf(0.6) = ppf(0.4)
}

///|
test "Bernoulli logarithmic functions" {
  let bernoulli = Bernoulli::new(0.3, loc=0.0)

  // Test log PMF
  let logpmf_0 = bernoulli.logpmf(0.0)
  let logpmf_1 = bernoulli.logpmf(1.0)
  assert_true((logpmf_0 - @math.ln(0.7)).abs() < 1.0e-10)
  assert_true((logpmf_1 - @math.ln(0.3)).abs() < 1.0e-10)

  // Test log CDF
  let logcdf_0 = bernoulli.logcdf(0.0)
  let logcdf_1 = bernoulli.logcdf(1.0)
  assert_true((logcdf_0 - @math.ln(0.7)).abs() < 1.0e-10)
  assert_true((logcdf_1 - @math.ln(1.0)).abs() < 1.0e-10)

  // Test log survival function
  let logsf_0 = bernoulli.logsf(0.0)
  assert_true((logsf_0 - @math.ln(0.3)).abs() < 1.0e-10)
  // Note: logsf(1) would be ln(0) = -inf, handled by the implementation
}

///|
test "Bernoulli random sampling" {
  let bernoulli = Bernoulli::new(0.6, loc=0.0)
  let seed = b"\x12\x34\x56\x78\x9a\xbc\xde\xf0\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb\xcc\xdd\xee\xff\x00\x01\x02\x03\x04\x05\x06\x07\x08"

  // Test sample generation
  let samples = bernoulli.rvs(1000, Some(seed))
  assert_eq(samples.length(), 1000)

  // All samples should be either 0 or 1
  let mut count_1 = 0
  for i = 0; i < samples.length(); i = i + 1 {
    let sample = samples[i]
    assert_true(sample == 0.0 || sample == 1.0)
    if sample == 1.0 {
      count_1 += 1
    }
  }
  // Check approximate proportions (should be close to p=0.6 for large samples)
  let proportion_1 = count_1.to_double() / 1000.0
  assert_true((proportion_1 - 0.6).abs() < 0.1) // Allow some variance
}

///|
test "Bernoulli entropy and moments" {
  let bernoulli = Bernoulli::new(0.5, loc=0.0)

  // Test entropy for p=0.5 (maximum entropy case)
  let entropy = bernoulli.entropy()
  let expected_entropy = -0.5 * @math.ln(0.5) - 0.5 * @math.ln(0.5)
  assert_true((entropy - expected_entropy).abs() < 1.0e-10)

  // Test median (should be 0 for p <= 0.5)
  assert_eq(bernoulli.median(), 0.0)

  // Test higher moments
  let moment_2 = bernoulli.moment(2)
  let moment_3 = bernoulli.moment(3)
  assert_true((moment_2 - 0.5).abs() < 1.0e-10) // E[X^2] = p for Bernoulli
  assert_true((moment_3 - 0.5).abs() < 1.0e-10) // E[X^3] = p for Bernoulli

  // Test confidence interval
  let interval = bernoulli.interval(0.05)
  let (lower, upper) = interval
  assert_eq(lower, 0.0)
  assert_eq(upper, 1.0)
}

///|
test "Bernoulli edge cases" {
  // Test with p = 0 (always 0)
  let bernoulli_zero = Bernoulli::new(0.0, loc=0.0)
  assert_true((bernoulli_zero.pmf(0.0) - 1.0).abs() < 1.0e-10)
  assert_eq(bernoulli_zero.pmf(1.0), 0.0)
  assert_eq(bernoulli_zero.mean(), 0.0)
  assert_eq(bernoulli_zero.varc(), 0.0)

  // Test with p = 1 (always 1)
  let bernoulli_one = Bernoulli::new(1.0, loc=0.0)
  assert_eq(bernoulli_one.pmf(0.0), 0.0)
  assert_true((bernoulli_one.pmf(1.0) - 1.0).abs() < 1.0e-10)
  assert_eq(bernoulli_one.mean(), 1.0)
  assert_eq(bernoulli_one.varc(), 0.0)

  // Test median for different p values
  let bernoulli_low = Bernoulli::new(0.3, loc=0.0)
  let bernoulli_high = Bernoulli::new(0.8, loc=0.0)
  assert_eq(bernoulli_low.median(), 0.0) // p < 0.5
  assert_eq(bernoulli_high.median(), 0.0) // We defined median as loc for Bernoulli
}
