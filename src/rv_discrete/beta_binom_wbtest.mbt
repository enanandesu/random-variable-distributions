///|
test "panic BetaBinom::new/invalid_n_negative" {
  ignore(BetaBinom::new(-1, 1.0, 1.0))
}

///|
test "panic BetaBinom::new/invalid_a_nonpositive" {
  ignore(BetaBinom::new(1, 0.0, 1.0))
}

///|
test "panic BetaBinom::new/invalid_b_nonpositive" {
  ignore(BetaBinom::new(1, 1.0, 0.0))
}

///|
test "BetaBinom::new creates valid distribution" {
  let bb = BetaBinom::new(10, 2.0, 3.0, loc=0.0)
  assert_eq(bb.n, 10)
  assert_eq(bb.a, 2.0)
  assert_eq(bb.b, 3.0)
  assert_eq(bb.loc, 0.0)
  let bb_with_loc = BetaBinom::new(5, 1.0, 1.0, loc=2.0)
  assert_eq(bb_with_loc.loc, 2.0)
}

///|
test "BetaBinom pmf basic properties" {
  let bb = BetaBinom::new(3, 2.0, 1.0)

  // PMF should be 0 for negative values
  assert_eq(bb.pmf(-1.0), 0.0)

  // PMF should be 0 for values greater than n + loc
  assert_eq(bb.pmf(4.0), 0.0)

  // PMF should be 0 for non-integer values
  assert_eq(bb.pmf(1.5), 0.0)

  // Sum of all PMF values should be approximately 1
  let mut total = 0.0
  for k in 0..=3 {
    total += bb.pmf(k.to_double())
  }
  assert_eq((total - 1.0).abs() < 1.0e-10, true)
}

///|
test "BetaBinom pmf with location parameter" {
  let bb = BetaBinom::new(5, 2.0, 1.0, loc=5.0)
  assert_eq(bb.pmf(4.0), 0.0)
  assert_true((bb.pmf(5.0) - 0.047619047619047616).abs() < 1.0e-10)
  assert_true((bb.pmf(6.0) - 0.09523809523809529).abs() < 1.0e-10)
  assert_true((bb.pmf(7.0) - 0.1428571428571429).abs() < 1.0e-10)
  assert_true((bb.pmf(8.0) - 0.1904761904761905).abs() < 1.0e-10)
  assert_true((bb.pmf(9.0) - 0.2380952380952381).abs() < 1.0e-10)
  assert_true((bb.pmf(10.0) - 0.2857142857142857).abs() < 1.0e-10)
  assert_eq(bb.pmf(11.0), 0.0)
}

///|
test "BetaBinom cdf properties" {
  let bb = BetaBinom::new(5, 2.0, 3.0)

  // CDF should be 0 for values below support
  assert_eq(bb.cdf(-1.0), 0.0)

  // CDF should be 1 for values above support
  assert_eq(bb.cdf(10.0), 1.0)

  // CDF should be monotonically increasing
  let mut prev_cdf = bb.cdf(0.0)
  for k in 1..=5 {
    let curr_cdf = bb.cdf(k.to_double())
    assert_true(curr_cdf >= prev_cdf)
    prev_cdf = curr_cdf
  }

  // CDF at the end should be 1
  assert_eq((bb.cdf(5.0) - 1.0).abs() < 1.0e-10, true)
  assert_true((bb.cdf(4.1) - 0.9523809523809522).abs() < 1.0e-10)
}

///|
test "BetaBinom ppf properties" {
  let bb = BetaBinom::new(4, 1.5, 2.5)

  // PPF at 0 should give value below support
  assert_eq(bb.ppf(0.0), -1.0)

  // PPF at 1 should give maximum value
  assert_eq(bb.ppf(1.0), 4.0)

  // PPF should be increasing
  let q1 = bb.ppf(0.25)
  let q2 = bb.ppf(0.5)
  let q3 = bb.ppf(0.75)
  let q4 = bb.ppf(0.999)
  assert_true(q1 <= q2)
  assert_true(q2 <= q3)
  assert_eq(q1, 0)
  assert_eq(q2, 1)
  assert_eq(q3, 2)
  assert_eq(q4, 4)
}

///|
test "BetaBinom sf and isf consistency" {
  let bb = BetaBinom::new(6, 2.0, 2.0)
  for k in 0..=6 {
    let x = k.to_double()
    // sf(x) should equal 1 - cdf(x)
    let sf_val = bb.sf(x)
    let expected_sf = 1.0 - bb.cdf(x)
    assert_eq((sf_val - expected_sf).abs() < 1.0e-10, true)
  }

  // isf should be inverse of sf (approximately)
  assert_eq(bb.isf(0.3), 4)
  assert_eq(bb.isf(0.7), 2)
}

///|
test "BetaBinom log functions" {
  let bb = BetaBinom::new(4, 1.0, 2.0)
  for k in 0..=4 {
    let x = k.to_double()

    // logpmf should equal ln(pmf)
    let pmf_val = bb.pmf(x)
    if pmf_val > 0.0 {
      let logpmf_val = bb.logpmf(x)
      let expected_logpmf = @math.ln(pmf_val)
      assert_eq((logpmf_val - expected_logpmf).abs() < 1.0e-10, true)
    }

    // logcdf should equal ln(cdf)
    let cdf_val = bb.cdf(x)
    if cdf_val > 0.0 {
      let logcdf_val = bb.logcdf(x)
      let expected_logcdf = @math.ln(cdf_val)
      assert_eq((logcdf_val - expected_logcdf).abs() < 1.0e-10, true)
    }

    // logsf should equal ln(sf)
    let sf_val = bb.sf(x)
    if sf_val > 0.0 {
      let logsf_val = bb.logsf(x)
      let expected_logsf = @math.ln(sf_val)
      assert_eq((logsf_val - expected_logsf).abs() < 1.0e-10, true)
    }
  }
}

///|
test "BetaBinom moments" {
  let bb = BetaBinom::new(10, 3.0, 2.0)

  // Test mean
  let theoretical_mean = bb.loc + bb.n.to_double() * bb.a / (bb.a + bb.b)
  let computed_mean = bb.mean()
  assert_eq((computed_mean - theoretical_mean).abs() < 1.0e-10, true)

  // Test variance
  let n = bb.n.to_double()
  let a = bb.a
  let b = bb.b
  let theoretical_var = n *
    a *
    b *
    (n + a + b) /
    ((a + b) * (a + b) * (a + b + 1.0))
  let computed_var = bb.varc()
  assert_eq((computed_var - theoretical_var).abs() < 1.0e-10, true)

  // Test standard deviation
  let computed_std = bb.std()
  let expected_std = computed_var.sqrt()
  assert_eq((computed_std - expected_std).abs() < 1.0e-10, true)
}

///|
test "BetaBinom entropy" {
  let bb = BetaBinom::new(3, 2.0, 3.0)
  let entropy = bb.entropy()
  assert_true((entropy - 1.3220643652636423).abs() < 1.0e-10)
}

///|
test "BetaBinom interval" {
  let bb = BetaBinom::new(5, 2.0, 3.0)
  let (lower, upper) = bb.interval(0.44)

  // Lower should be less than or equal to upper
  assert_true(lower <= upper)
  assert_true((lower - 1).abs() < 1.0e-10)
  assert_true((upper - 3).abs() < 1.0e-10)
}

///|
test "BetaBinom rvs generates correct size" {
  let bb = BetaBinom::new(5, 2.0, 1.0)
  let samples = bb.rvs(
    100,
    rand=@random.Rand::chacha8(
      seed=b"\x12\x34\x56\x78\x9a\xbc\xde\xf0\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb\xcc\xdd\xee\xff\x00\x01\x02\x03\x04\x05\x06\x07\x08",
    ),
  )
  assert_eq(samples.length(), 100)

  // All samples should be within valid range
  for i in 0..<samples.length() {
    let sample = samples[i]
    assert_true(sample >= bb.loc)
    assert_true(sample <= bb.loc + bb.n.to_double())
    // Should be integer values (when loc is integer)
    assert_eq(sample, sample.round())
  }
}

///|
test "BetaBinom rvs with different seeds" {
  let bb = BetaBinom::new(3, 1.0, 1.0)
  let samples1 = bb.rvs(
    50,
    rand=@random.Rand::chacha8(
      seed=b"\x12\x34\x56\x78\x9a\xbc\xde\xf0\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb\xcc\xdd\xee\xff\x00\x01\x02\x03\x04\x05\x06\x07\x08",
    ),
  )
  let samples2 = bb.rvs(
    50,
    rand=@random.Rand::chacha8(
      seed=b"\x12\x34\x56\x78\x9a\xbc\xde\xf0\x11\x22\x33\x44\x55\x66\x77\x88\x99\xaa\xbb\xcc\xdd\xee\xff\x00\x01\x02\x03\x04\x05\x06\x07\x08",
    ),
  )
  let samples3 = bb.rvs(
    50,
    rand=@random.Rand::chacha8(seed=b"12345678901234567890987654321012"),
  )

  // Same seed should produce same results
  for i in 0..<50 {
    assert_eq(samples1[i], samples2[i])
  }

  // Different seeds should produce different results (with high probability)
  let mut differences = 0
  for i in 0..<50 {
    if samples1[i] != samples3[i] {
      differences += 1
    }
  }
  assert_true(differences > 10) // Should have many differences
  let samples4 = bb.rvs(
    1000,
    rand=@random.Rand::chacha8(seed=b"12345678901234567890987654321012"),
  )
  let mut avg = 0.0
  for sample in samples4 {
    avg += sample
  }
  avg /= 1000
  inspect(avg - bb.mean().abs(), content="-0.028000000000000025")
}

///|
test "BetaBinom edge cases" {
  // Test with very small parameters
  let bb_small = BetaBinom::new(1, 0.1, 0.1)
  assert_true(bb_small.pmf(0.0) > 0.0)
  assert_true(bb_small.pmf(1.0) > 0.0)

  // Test with n = 0 (degenerate case)
  let bb_zero = BetaBinom::new(0, 1.0, 1.0)
  assert_eq(bb_zero.pmf(0.0), 1.0)
  assert_eq(bb_zero.pmf(1.0), 0.0)

  // Test with large location parameter
  let bb_large_loc = BetaBinom::new(2, 1.0, 1.0, loc=100.0)
  assert_eq(bb_large_loc.pmf(99.0), 0.0)
  assert_true(bb_large_loc.pmf(100.0) > 0.0)
  assert_true(bb_large_loc.pmf(102.0) > 0.0)
  assert_eq(bb_large_loc.pmf(103.0), 0.0)
}

///|
test "BetaBinom skew and kurt" {
  let bb = BetaBinom::new(5, 1.2, 2.1, loc=3)
  assert_true((bb.skew() - 0.4579606781971745).abs() < 1.0e-10)
  assert_true((bb.kurt() - -0.8049953375303294).abs() < 1.0e-10)
}
