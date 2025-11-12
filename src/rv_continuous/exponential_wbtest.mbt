///|
test "panic Exponential::new/invalid_scale" {
  ignore(Exponential::new(scale=0.0))
}

///|
test "Exponential::pdf_cdf" {
  let dist = Exponential::new(loc=1.0, scale=2.0)
  assert_true((dist.pdf(0.5) - 0.0).abs() < 1.0e-12)
  assert_true((dist.pdf(1.0) - 0.5).abs() < 1.0e-12)
  assert_true((dist.pdf(3.0) - 0.18393972058572117).abs() < 1.0e-12)
  assert_true((dist.cdf(0.5) - 0.0).abs() < 1.0e-12)
  assert_true((dist.cdf(1.0) - 0.0).abs() < 1.0e-12)
  assert_true((dist.cdf(3.0) - 0.6321205588285577).abs() < 1.0e-12)
  assert_true((dist.sf(0.5) - 1.0).abs() < 1.0e-12)
  assert_true((dist.sf(3.0) - 0.36787944117144233).abs() < 1.0e-12)
}

///|
test "panic Exponential::ppf/invalid_q_low" {
  let dist = Exponential::new()
  ignore(dist.ppf(0.0))
}

///|
test "panic Exponential::ppf/invalid_q_high" {
  let dist = Exponential::new()
  ignore(dist.ppf(1.0))
}

///|
test "panic Exponential::isf/invalid_q_low" {
  let dist = Exponential::new()
  ignore(dist.isf(0.0))
}

///|
test "panic Exponential::isf/invalid_q_high" {
  let dist = Exponential::new()
  ignore(dist.isf(1.0))
}

///|
test "Exponential::ppf_isf" {
  let dist = Exponential::new(loc=1.5, scale=3.0)
  let q1 = 0.25
  let q2 = 0.75
  let expected_ppf1 = 2.363046217355343
  let expected_ppf2 = 5.6588830833596715
  assert_true((dist.ppf(q1) - expected_ppf1).abs() < 1.0e-12)
  assert_true((dist.ppf(q2) - expected_ppf2).abs() < 1.0e-12)
  assert_true((dist.isf(0.25) - dist.ppf(0.75)).abs() < 1.0e-12)
}

///|
test "Exponential::stats" {
  let dist = Exponential::new(loc=-2.0, scale=4.0)
  let expected_mean = 2.0
  let expected_median = 0.7725887222397811
  let expected_var = 16.0
  let expected_entropy = 2.386294361119891
  assert_true((dist.mean() - expected_mean).abs() < 1.0e-12)
  assert_true((dist.median() - expected_median).abs() < 1.0e-12)
  assert_true((dist.varc() - expected_var).abs() < 1.0e-12)
  assert_true((dist.std() - expected_var.sqrt()).abs() < 1.0e-12)
  assert_true((dist.skew() - 2.0).abs() < 1.0e-12)
  assert_true((dist.kurt() - 6.0).abs() < 1.0e-12)
  assert_true((dist.entropy() - expected_entropy).abs() < 1.0e-12)
}

///|
test "Exponential::moment" {
  let dist = Exponential::new(loc=1.5, scale=2.0)
  assert_true((dist.moment(0) - 1.0).abs() < 1.0e-12)
  assert_true((dist.moment(1) - dist.mean()).abs() < 1.0e-12)
  let expected_second = 16.25
  assert_true((dist.moment(2) - expected_second).abs() < 1.0e-12)
  let expected_fifth = 8128.21875
  assert_true((dist.moment(5) - expected_fifth).abs() < 1.0e-12)
}

///|
test "Exponential::log_functions" {
  let dist = Exponential::new(loc=0.5, scale=1.5)
  assert_eq(dist.logpdf(0.0), @double.neg_infinity)
  assert_true((dist.logpdf(1.5) - @math.ln(dist.pdf(1.5))).abs() < 1.0e-12)
  assert_eq(dist.logcdf(0.0), @double.neg_infinity)
  assert_true((dist.logcdf(2.0) - @math.ln(dist.cdf(2.0))).abs() < 1.0e-12)
  assert_true((dist.logsf(0.0) - 0.0).abs() < 1.0e-12)
  assert_true((dist.logsf(2.0) - @math.ln(dist.sf(2.0))).abs() < 1.0e-12)
}

///|
test "Exponential::interval" {
  let dist = Exponential::new(loc=-1.0, scale=0.8)
  let alpha = 0.2
  let (lower, upper) = dist.interval(alpha)
  assert_true((lower - -0.9157115874737389).abs() < 1.0e-12)
  assert_true((upper - 0.8420680743952369).abs() < 1.0e-12)
  assert_true(lower < upper)
}

///|
test "Exponential::rvs" {
  let dist = Exponential::new(loc=0.5, scale=1.5)
  let seed = b"12345678901234567890123456789012"
  let samples = dist.rvs(50000, rand=@random.Rand::chacha8(seed~))
  assert_eq(samples.length(), 50000)
  let mut sum = 0.0
  for sample in samples {
    assert_true(sample >= dist.loc)
    sum += sample
  }
  let mean_estimate = sum / 50000.0
  assert_true((mean_estimate - dist.mean()).abs() < 5.0e-2)
}
