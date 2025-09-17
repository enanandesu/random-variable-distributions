///|
test "Boltzmann::pmf" {
  let dist = Boltzmann::new(0.5, 5)
  let lambda = 0.5
  let q = @math.exp(-lambda)
  let norm = -@math.expm1(-lambda) / -@math.expm1(-lambda * dist.n.to_double())
  let mut weight = 1.0
  for idx = 0; idx < dist.n; idx = idx + 1 {
    let expected = norm * weight
    let value = dist.pmf(idx.to_double())
    assert_true((value - expected).abs() < 1.0e-12)
    weight *= q
  }
  inspect(dist.pmf(-1.0), content="0")
  inspect(dist.pmf(5.0), content="0")
  inspect(dist.pmf(1.5), content="0")
}

///|
test "Boltzmann::pmf/location" {
  let dist = Boltzmann::new(0.7, 4, loc=2.0)
  assert_true((dist.pmf(2.0) - 0.5360094653651551).abs() < 1.0e-12)
  inspect(dist.pmf(1.0), content="0")
  inspect(dist.pmf(6.0), content="0")
}

///|
test "Boltzmann::cdf" {
  let dist = Boltzmann::new(1.2, 4, loc=3.0)
  assert_true((dist.cdf(2.9) - 0.0).abs() < 1.0e-9)
  let cdf_loc = dist.cdf(3.0)
  let pmf_loc = dist.pmf(3.0)
  assert_true((cdf_loc - pmf_loc).abs() < 1.0e-9)
  assert_true((dist.cdf(3.6) - cdf_loc).abs() < 1.0e-9)
  assert_true((dist.cdf(3.9) - cdf_loc).abs() < 1.0e-9)
  let pmf_next = dist.pmf(4.0)
  let cdf_next = dist.cdf(4.0)
  assert_true((cdf_next - (cdf_loc + pmf_next)).abs() < 1.0e-9)
  assert_true((dist.cdf(4.2) - cdf_next).abs() < 1.0e-9)
  let pmf_third = dist.pmf(5.0)
  let cdf_third = dist.cdf(5.0)
  assert_true((cdf_third - (cdf_next + pmf_third)).abs() < 1.0e-9)
  assert_true(dist.cdf(3.0) <= dist.cdf(4.0))
  assert_true(dist.cdf(4.0) <= dist.cdf(5.0))
  assert_true((dist.cdf(10.0) - 1.0).abs() < 1.0e-9)
}

///|
test "Boltzmann::ppf" {
  let dist = Boltzmann::new(0.9, 5, loc=1.0)
  inspect(dist.ppf(0.0), content="0")
  inspect(dist.ppf(1.0), content="5")
  let expected_0_2 = 1.0
  let expected_0_6 = 1.0
  let expected_0_9 = 3.0
  assert_true((dist.ppf(0.2) - expected_0_2).abs() < 1.0e-9)
  assert_true((dist.ppf(0.6) - expected_0_6).abs() < 1.0e-9)
  assert_true((dist.ppf(0.9) - expected_0_9).abs() < 1.0e-9)
}

///|
test "panic Boltzmann::ppf/invalid_low" {
  let dist = Boltzmann::new(1.0, 3)
  ignore(dist.ppf(-0.1))
}

///|
test "panic Boltzmann::ppf/invalid_high" {
  let dist = Boltzmann::new(1.0, 3)
  ignore(dist.ppf(1.1))
}

///|
test "Boltzmann::stats" {
  let dist = Boltzmann::new(0.8, 6, loc=2.0)
  let mean_expected = 2.766177994109335
  assert_true((dist.mean() - mean_expected).abs() < 1.0e-9)
  let var_expected = 1.180558866223071
  assert_true((dist.varc() - var_expected).abs() < 1.0e-9)
  assert_true((dist.std() - var_expected.sqrt()).abs() < 1.0e-9)
  let entropy_expected = 1.2012962761079906
  assert_true((dist.entropy() - entropy_expected).abs() < 1.0e-9)
  let skew_expected = 1.6080439600167251
  let kurt_expected = 2.344600506417025
  assert_true((dist.skew() - skew_expected).abs() < 1.0e-9)
  assert_true((dist.kurt() - kurt_expected).abs() < 1.0e-9)
}

///|
test "Boltzmann::log_functions" {
  let dist = Boltzmann::new(1.1, 5, loc=1.0)
  let lambda = 1.1
  let norm = -@math.expm1(-lambda) / -@math.expm1(-lambda * dist.n.to_double())
  let log_norm = @math.ln(norm)
  assert_true((dist.logpmf(1.0) - log_norm).abs() < 1.0e-9)
  assert_true((dist.logpmf(2.0) - (log_norm - lambda)).abs() < 1.0e-9)
  assert_true(dist.logpmf(0.5) == @double.neg_infinity)
  assert_true(dist.logpmf(10.0) == @double.neg_infinity)
  assert_true(dist.logcdf(0.0) == @double.neg_infinity)
  let cdf_value = dist.cdf(2.5)
  assert_true(cdf_value > 0.0)
  let expected_logcdf = @math.ln(cdf_value)
  assert_true((dist.logcdf(2.5) - expected_logcdf).abs() < 1.0e-9)
  let sf_value = dist.sf(1.5)
  let logsf = if sf_value <= 0.0 {
    @double.neg_infinity
  } else {
    @math.ln(sf_value)
  }
  assert_true((dist.logsf(1.5) - logsf).abs() < 1.0e-9)
  assert_true(dist.logsf(6.0) == @double.neg_infinity)
}

///|
test "Boltzmann::interval_and_isf" {
  let dist = Boltzmann::new(0.9, 6, loc=1.0)
  let (lower, upper) = dist.interval(0.2)
  assert_true(lower >= dist.loc - 1.0)
  assert_true(upper >= lower)
  assert_true((dist.isf(0.3) - dist.ppf(0.7)).abs() < 1.0e-9)
}

///|
test "Boltzmann::degenerate_case" {
  let dist = Boltzmann::new(0.7, 1, loc=4.0)
  assert_true((dist.pmf(4.0) - 1.0).abs() < 1.0e-9)
  inspect(dist.pmf(5.0), content="0")
  assert_true((dist.cdf(3.5) - 0.0).abs() < 1.0e-9)
  assert_true((dist.cdf(4.0) - 1.0).abs() < 1.0e-9)
  assert_true((dist.mean() - 4.0).abs() < 1.0e-9)
  assert_true((dist.varc() - 0.0).abs() < 1.0e-9)
  assert_true((dist.entropy() - 0.0).abs() < 1.0e-9)
  let samples = dist.rvs(10, None)
  for value in samples {
    assert_true((value - 4.0).abs() < 1.0e-9)
  }
}

///|
test "Boltzmann::rvs" {
  let dist = Boltzmann::new(0.6, 6)
  let seed = b"12345678901234567890123456789012"
  let samples = dist.rvs(50000, Some(seed))
  assert_eq(samples.length(), 50000)
  let counts = Array::make(dist.n, 0)
  for value in samples {
    let diff = value - dist.loc
    assert_true(diff >= 0.0)
    let rounded = diff.round()
    assert_true((diff - rounded).abs() < 1.0e-9)
    let idx = rounded.to_int()
    assert_true(idx >= 0)
    assert_true(idx < dist.n)
    counts[idx] += 1
  }
  let total = samples.length().to_double()
  for idx = 0; idx < dist.n; idx = idx + 1 {
    let observed = counts[idx].to_double() / total
    let expected = dist.pmf(dist.loc + idx.to_double())
    assert_true((observed - expected).abs() < 1.0e-2)
  }
}
