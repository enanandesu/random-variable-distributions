///|
test "panic DLaplace::new/invalid_a" {
  ignore(DLaplace::new(0.0))
}

///|
test "DLaplace::pmf" {
  let a = 0.8
  let dist = DLaplace::new(a, loc=1.0)
  let pmf_loc = dist.pmf(1.0)
  let pmf_plus1 = dist.pmf(2.0)
  let pmf_minus1 = dist.pmf(0.0)
  assert_true((pmf_loc - 0.3799489622552249).abs() < 1.0e-12)
  assert_true((pmf_plus1 - 0.1707220736275535).abs() < 1.0e-12)
  assert_true((pmf_minus1 - 0.1707220736275535).abs() < 1.0e-12)
  assert_true((dist.pmf(1.5) - 0.0).abs() < 1.0e-12)
}

///|
test "DLaplace::cdf" {
  let a = 0.6
  let dist = DLaplace::new(a, loc=-2.0)
  assert_true((dist.cdf(-3.1) - 0.1944679423198219).abs() < 1.0e-12)
  assert_true((dist.cdf(-2.0) - 0.6456563062257954).abs() < 1.0e-12)
  let expected = 0.8932737304076198
  assert_true((dist.cdf(0.0) - expected).abs() < 1.0e-12)
}

///|
test "panic DLaplace::ppf/invalid_q_low" {
  let dist = DLaplace::new(0.9, loc=0.0)
  ignore(dist.ppf(0.0))
}

///|
test "panic DLaplace::ppf/invalid_q_high" {
  let dist = DLaplace::new(0.9, loc=0.0)
  ignore(dist.ppf(1.0))
}

///|
test "panic DLaplace::isf/invalid_q_low" {
  let dist = DLaplace::new(0.9, loc=0.0)
  ignore(dist.isf(0.0))
}

///|
test "panic DLaplace::isf/invalid_q_high" {
  let dist = DLaplace::new(0.9, loc=0.0)
  ignore(dist.isf(1.0))
}

///|
test "DLaplace::ppf" {
  let a = 0.9
  let dist = DLaplace::new(a, loc=0.0)
  assert_true(dist.ppf(0.1) == -2.0)
  assert_true(dist.ppf(0.3) == 0.0)
  assert_true(dist.ppf(0.5) == 0.0)
  assert_true(dist.ppf(0.7) == 0.0)
  assert_true(dist.ppf(0.9) == 2.0)
  assert_true(dist.ppf(0.99) == 4.0)
  assert_true(dist.ppf(0.999) == 7.0)
  assert_true(dist.ppf(1.0 - 1.0e-10) == 25.0)
}

///|
test "DLaplace::stats" {
  let a = 1.1
  let dist = DLaplace::new(a, loc=3.0)
  let expected_mean = 3.0
  let expected_var = 1.495844796352319
  let expected_std = 1.2230473401926514
  let expected_kurt = 3.668518553822257
  let expected_entropy = 1.5156779077596056
  assert_true((dist.mean() - expected_mean).abs() < 1.0e-12)
  assert_true((dist.varc() - expected_var).abs() < 1.0e-12)
  assert_true((dist.std() - expected_std).abs() < 1.0e-12)
  assert_true((dist.skew() - 0.0).abs() < 1.0e-12)
  assert_true((dist.kurt() - expected_kurt).abs() < 1.0e-12)
  assert_true((dist.entropy() - expected_entropy).abs() < 1.0e-12)
}

///|
test "DLaplace::log_functions" {
  let dist = DLaplace::new(0.7, loc=2.0)
  assert_true((dist.logpmf(2.0) - -1.0895270516938431).abs() < 1.0e-12)
  assert_true((dist.logpmf(3.0) - -1.789527051693843).abs() < 1.0e-12)
  assert_true((dist.logcdf(2.0) - -0.40318604888545784).abs() < 1.0e-12)
  assert_true(dist.logpmf(2.5) == @double.neg_infinity)
}

///|
test "DLaplace::rvs" {
  let dist = DLaplace::new(0.5, loc=0.0)
  let seed = b"12345678901234567890123456789012"
  let samples = dist.rvs(50000, rand=@random.Rand::chacha8(seed~))
  assert_eq(samples.length(), 50000)
  let max_k = 10
  let size = 2 * max_k + 1
  let counts = Array::make(size, 0)
  for sample in samples {
    assert_true(sample == sample.round())
    let k = sample.to_int()
    if k >= -max_k && k <= max_k {
      counts[k + max_k] += 1
    }
  }
  let total = 50000.0
  for idx = 0; idx < size; idx = idx + 1 {
    let k = idx - max_k
    let observed = counts[idx].to_double() / total
    let expected = dist.pmf(k.to_double())
    assert_true((observed - expected).abs() < 3.0e-2)
  }
}
