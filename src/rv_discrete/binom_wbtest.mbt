///|
test "Binom creation and basic properties" {
  let binom = Binom::new(10, 0.3, loc=2.0)

  // Test basic properties
  assert_eq(binom.n, 10)
  assert_eq(binom.p, 0.3)
  assert_eq(binom.loc, 2.0)

  // Test default location
  let binom_default = Binom::new(5, 0.5)
  assert_eq(binom_default.loc, 0.0)
}

///|
test "Binom PMF properties" {
  let binom = Binom::new(3, 0.5, loc=0.0)

  // Test PMF at valid points
  let pmf_0 = binom.pmf(0.0) // C(3,0) * 0.5^0 * 0.5^3 = 1 * 1 * 0.125 = 0.125
  let pmf_1 = binom.pmf(1.0) // C(3,1) * 0.5^1 * 0.5^2 = 3 * 0.5 * 0.25 = 0.375
  let pmf_2 = binom.pmf(2.0) // C(3,2) * 0.5^2 * 0.5^1 = 3 * 0.25 * 0.5 = 0.375
  let pmf_3 = binom.pmf(3.0) // C(3,3) * 0.5^3 * 0.5^0 = 1 * 0.125 * 1 = 0.125
  assert_true((pmf_0 - 0.125).abs() < 1.0e-10)
  assert_true((pmf_1 - 0.375).abs() < 1.0e-10)
  assert_true((pmf_2 - 0.375).abs() < 1.0e-10)
  assert_true((pmf_3 - 0.125).abs() < 1.0e-10)

  // Test PMF at non-integer points
  assert_eq(binom.pmf(0.5), 0.0)
  assert_eq(binom.pmf(1.7), 0.0)

  // Test PMF outside support
  assert_eq(binom.pmf(-1.0), 0.0)
  assert_eq(binom.pmf(4.0), 0.0)

  // Test PMF sums to 1
  let total = pmf_0 + pmf_1 + pmf_2 + pmf_3
  assert_true((total - 1.0).abs() < 1.0e-10)
}

///|
test "Binom CDF properties" {
  let binom = Binom::new(3, 0.5, loc=0.0)

  // Test CDF is non-decreasing
  let cdf_0 = binom.cdf(0.0)
  let cdf_1 = binom.cdf(1.0)
  let cdf_2 = binom.cdf(2.0)
  let cdf_3 = binom.cdf(3.0)
  assert_true(cdf_0 <= cdf_1)
  assert_true(cdf_1 <= cdf_2)
  assert_true(cdf_2 <= cdf_3)

  // Test CDF at boundaries
  assert_eq(binom.cdf(-1.0), 0.0)
  assert_true((binom.cdf(3.0) - 1.0).abs() < 1.0e-10)
  assert_true((binom.cdf(10.0) - 1.0).abs() < 1.0e-10)

  // Test CDF at specific points
  assert_true((cdf_0 - 0.125).abs() < 1.0e-10)
  assert_true((cdf_1 - 0.5).abs() < 1.0e-10)
  assert_true((cdf_2 - 0.875).abs() < 1.0e-10)
}

///|
test "Binom statistical moments" {
  let binom = Binom::new(10, 0.3, loc=0.0)

  // Test mean = n * p
  let expected_mean = 10.0 * 0.3
  assert_true((binom.mean() - expected_mean).abs() < 1.0e-10)

  // Test variance = n * p * (1-p)
  let expected_var = 10.0 * 0.3 * 0.7
  assert_true((binom.varc() - expected_var).abs() < 1.0e-10)

  // Test standard deviation
  let expected_std = expected_var.sqrt()
  assert_true((binom.std() - expected_std).abs() < 1.0e-10)

  // Test skewness = (1-2p) / sqrt(np(1-p))
  let expected_skew = (1.0 - 2.0 * 0.3) / (10.0 * 0.3 * 0.7).sqrt()
  assert_true((binom.skew() - expected_skew).abs() < 1.0e-10)

  // Test kurtosis = (1 - 6p(1-p)) / (np(1-p))
  let expected_kurt = (1.0 - 6.0 * 0.3 * 0.7) / (10.0 * 0.3 * 0.7)
  assert_true((binom.kurt() - expected_kurt).abs() < 1.0e-10)
}

///|
test "Binom with location parameter" {
  let binom = Binom::new(2, 0.5, loc=5.0)

  // Test PMF with location shift
  assert_eq(binom.pmf(5.0), 0.25) // Should be C(2,0) * 0.5^2 = 0.25
  assert_eq(binom.pmf(6.0), 0.5) // Should be C(2,1) * 0.5^2 = 0.5
  assert_eq(binom.pmf(7.0), 0.25) // Should be C(2,2) * 0.5^2 = 0.25
  assert_eq(binom.pmf(4.0), 0.0) // Outside support
  assert_eq(binom.pmf(8.0), 0.0) // Outside support

  // Test mean with location
  let expected_mean = 5.0 + 2.0 * 0.5
  assert_true((binom.mean() - expected_mean).abs() < 1.0e-10)
}

///|
test "Binom PPF and quantiles" {
  let binom = Binom::new(4, 0.5, loc=0.0)

  // Test PPF at boundaries
  assert_eq(binom.ppf(0.0), -1.0)
  assert_eq(binom.ppf(1.0), 4.0)

  // Test PPF consistency with CDF
  let x = binom.ppf(0.5)
  assert_true(binom.cdf(x) >= 0.5)
  assert_true(binom.cdf(x - 1.0) < 0.5)
}

///|
test "Binom random sampling" {
  let binom = Binom::new(5, 0.4, loc=0.0)
  let seed = b"\x12\x34\x56\x78\x9a\xbc\xde\xf0\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb\xcc\xdd\xee\xff\x00\x01\x02\x03\x04\x05\x06\x07\x08"

  // Test sample generation
  let samples = binom.rvs(100, Some(seed))
  assert_eq(samples.length(), 100)

  // All samples should be within valid range
  let mut avg = 0.0
  for sample in samples {
    assert_true(sample >= 0.0)
    assert_true(sample <= 5.0)
    // Check if sample is integer (difference with rounded value should be 0)
    assert_true((sample - sample.round()).abs() < 1.0e-10)
    avg += sample
  }
  avg /= 100.0
  inspect(avg - binom.mean(), content="-0.07000000000000006")
}

///|
test "Binom entropy and interval" {
  let binom = Binom::new(3, 0.5, loc=0.0)

  // Test entropy is positive
  let entropy = binom.entropy()
  assert_true((entropy - 1.2554823251787535).abs() < 1.0e-10)

  // Test confidence interval
  let interval = binom.interval(0.05)
  let (lower, upper) = interval
  assert_true(lower <= upper)
  assert_eq(lower, 0.0)
  assert_eq(upper, 3.0)
}

///|
test "Binom edge cases" {
  // Test with p = 0 (degenerate case)
  let binom_zero = Binom::new(5, 0.0, loc=0.0)
  assert_true((binom_zero.pmf(0.0) - 1.0).abs() < 1.0e-10)
  assert_eq(binom_zero.pmf(1.0), 0.0)
  assert_eq(binom_zero.mean(), 0.0)

  // Test with p = 1 (degenerate case)
  let binom_one = Binom::new(3, 1.0, loc=0.0)
  assert_eq(binom_one.pmf(0.0), 0.0)
  assert_true((binom_one.pmf(3.0) - 1.0).abs() < 1.0e-10)
  assert_eq(binom_one.mean(), 3.0)
}
