///|
test "panic F::new/invalid_params" {
  ignore(F::new(-1.0, 5.0))
  ignore(F::new(3.0, 0.0))
  ignore(F::new(3.0, 5.0, scale=0.0))
}

///|
test "F::pdf_cdf" {
  let dist = F::new(5.0, 7.0, loc=1.0, scale=2.0)
  assert_true((dist.pdf(1.0) - 0.0).abs() < 1.0e-12)
  assert_true((dist.pdf(3.0) - 0.23073820875755294).abs() < 1.0e-12)
  assert_true((dist.pdf(5.0) - 0.0807343695696608).abs() < 1.0e-12)
  assert_true((dist.cdf(1.0) - 0.0).abs() < 1.0e-12)
  assert_true((dist.cdf(3.0) - 0.5187060870417945).abs() < 1.0e-12)
  assert_true((dist.cdf(5.0) - 0.8043267551331767).abs() < 1.0e-12)
  assert_true((dist.sf(5.0) - 0.1956732448668233).abs() < 1.0e-12)
}

///|
test "panic F::ppf/invalid_q" {
  let dist = F::new(5.0, 9.0)
  ignore(dist.ppf(0.0))
  ignore(dist.ppf(1.0))
}

///|
test "panic F::isf/invalid_q" {
  let dist = F::new(5.0, 9.0)
  ignore(dist.isf(0.0))
  ignore(dist.isf(1.0))
}

///|
test "F::ppf_isf" {
  let dist = F::new(5.0, 7.0, loc=1.0, scale=2.0)
  let ppf_025 = dist.ppf(0.25)
  let ppf_075 = dist.ppf(0.75)
  assert_true((ppf_025 - 2.0562379630764443).abs() < 1.0e-12)
  assert_true((ppf_075 - 4.422114006021369).abs() < 1.0e-12)
  assert_true((dist.cdf(ppf_025) - 0.25).abs() < 1.0e-12)
  assert_true((dist.cdf(ppf_075) - 0.75).abs() < 1.0e-12)
  assert_true((dist.isf(0.25) - dist.ppf(0.75)).abs() < 1.0e-12)
}

///|
test "F::stats" {
  let dist = F::new(5.0, 12.0, loc=-1.5, scale=3.0)
  assert_true((dist.mean() - 2.1).abs() < 1.0e-12)
  assert_true((dist.varc() - 9.72).abs() < 1.0e-12)
  assert_true((dist.skew() - 3.079201435678004).abs() < 1.0e-9)
  assert_true((dist.kurt() - 24.333333333333332).abs() < 1.0e-9)
  assert_true((dist.entropy() - 2.18020169141058).abs() < 1.0e-12)
}

///|
test "F::moment" {
  let dist = F::new(5.0, 12.0, loc=-1.5, scale=3.0)
  assert_true((dist.moment(0) - 1.0).abs() < 1.0e-12)
  assert_true((dist.moment(1) - dist.mean()).abs() < 1.0e-12)
  assert_true((dist.moment(2) - 14.129999999999999).abs() < 1.0e-10)
  assert_true((dist.moment(3) - 163.809).abs() < 1.0e-10)
}

///|
test "F::log_functions" {
  let dist = F::new(5.0, 7.0, loc=1.0, scale=2.0)
  assert_eq(dist.logpdf(1.0), @double.neg_infinity)
  assert_true((dist.logpdf(3.0) + 1.4664715065470475).abs() < 1.0e-12)
  assert_true((dist.logpdf(5.0) + 2.516590901316424).abs() < 1.0e-12)
  assert_eq(dist.logcdf(1.0), @double.neg_infinity)
  assert_true((dist.logcdf(5.0) + 0.2177496805089156).abs() < 1.0e-12)
  assert_true((dist.logsf(5.0) + 1.6313091289656863).abs() < 1.0e-12)
}

///|
test "F::interval" {
  let dist = F::new(5.0, 10.0, loc=0.5, scale=1.5)
  let alpha = 0.2
  let (lower, upper) = dist.interval(alpha)
  assert_true((lower - 0.9549036335316072).abs() < 1.0e-12)
  assert_true((upper - 4.282461029314436).abs() < 1.0e-12)
}

///|
test "F::rvs" {
  let dist = F::new(4.0, 9.0, loc=0.5, scale=2.5)
  let seed = b"12345678901234567890123456789012"
  let samples = dist.rvs(50000, Some(seed))
  assert_eq(samples.length(), 50000)
  let mut sum = 0.0
  let mut min_val = @double.infinity
  for sample in samples {
    assert_true(sample >= dist.loc)
    sum += sample
    if sample < min_val {
      min_val = sample
    }
  }
  let mean_estimate = sum / 50000.0
  assert_true((mean_estimate - dist.mean()).abs() < 5.0e-2)
  assert_true(min_val >= dist.loc - 1.0e-9)
}
