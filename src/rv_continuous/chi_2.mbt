///|
struct Chi2 {
  k : Double
}

///|
pub fn Chi2::new(k : Double) -> Chi2 {
  if k <= 0.0 {
    abort("k must be positive")
  }
  { k, }
}

///|
impl RVContinuous for Chi2 with rvs(self, size, seed) {
  let mut samples = Array::new(capacity=size)
  let normal = Normal::new()
  let mut seed_val = match seed {
    None => b"ABCDEFGHIJKLMNOPQRSTUVWXYZ123456"
    Some(val) => val
  }
  if self.k.floor() == self.k {
    for _ in 0..<size {
      let mut chi2_sample = 0.0
      let k_int = self.k.to_int()
      for _ in 0..<k_int {
        let z = normal.rvs(1, Some(seed_val))[0]
        seed_val = derive_seed(seed_val)
        chi2_sample += z * z
      }
      samples.push(chi2_sample)
    }
  } else {
    samples = Gamma::new(self.k / 2.0, theta=2.0).rvs(size, seed)
  }
  samples
}

///|
impl RVContinuous for Chi2 with pdf(self, x) {
  if x <= 0.0 {
    return 0.0
  }
  let k_half = self.k / 2.0
  let numerator = @math.pow(x, k_half - 1.0) * @math.exp(-x / 2.0)
  let denominator = @math.pow(2.0, k_half) * @kmath.gamma(k_half)
  numerator / denominator
}

///|
impl RVContinuous for Chi2 with cdf(self, x) {
  if x <= 0.0 {
    return 0.0
  }
  @sp.gammainc(self.k / 2.0, x / 2.0)
}

///|
impl RVContinuous for Chi2 with ppf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  let ans = @sp.func_inv(
    q,
    self.cdf(_),
    self.pdf(_),
    self.k,
    max_iterations=200,
  )
  match ans {
    None => abort("No solution found for ppf")
    Some(val) => val
  }
}

///|
impl RVContinuous for Chi2 with sf(self, x) {
  if x <= 0.0 {
    return 1.0
  }
  1.0 - self.cdf(x)
}

///|
impl RVContinuous for Chi2 with isf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  self.ppf(1.0 - q)
}

///|
impl RVContinuous for Chi2 with logpdf(self, x) {
  @math.ln(self.pdf(x))
}

///|
impl RVContinuous for Chi2 with logcdf(self, x) {
  @math.ln(self.cdf(x))
}

///|
impl RVContinuous for Chi2 with logsf(self, x) {
  @math.ln(self.sf(x))
}

///|
impl RVContinuous for Chi2 with moment(self, n) {
  if n < 0 {
    abort("moment order must be non-negative")
  }
  if n == 0 {
    return 1.0
  }
  @math.pow(2.0, n.to_double()) *
  @kmath.gamma(self.k / 2.0 + n.to_double()) /
  @kmath.gamma(self.k / 2.0)
}

///|
impl RVContinuous for Chi2 with mean(self) {
  self.k
}

///|
impl RVContinuous for Chi2 with median(self) {
  self.ppf(0.5)
}

///|
impl RVContinuous for Chi2 with varc(self) {
  2.0 * self.k
}

///|
impl RVContinuous for Chi2 with std(self) {
  (2.0 * self.k).sqrt()
}

///|
impl RVContinuous for Chi2 with entropy(self) {
  let k_half = self.k / 2.0
  k_half +
  @math.ln(2.0) +
  @kmath.ln_gamma(k_half) +
  (1.0 - k_half) * @kmath.digamma(k_half)
}

///|
impl RVContinuous for Chi2 with interval(self, alpha : Double) {
  if alpha <= 0.0 || alpha >= 1.0 {
    abort("alpha must be in (0, 1)")
  }
  let lower = self.ppf(alpha / 2.0)
  let upper = self.ppf(1.0 - alpha / 2.0)
  (lower, upper)
}

///|
impl RVContinuous for Chi2 with skew(self) {
  (8.0 / self.k).sqrt()
}

///|
impl RVContinuous for Chi2 with kurt(self) {
  12.0 / self.k
}
