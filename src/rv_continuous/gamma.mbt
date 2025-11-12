///|
struct Gamma {
  k : Double
  theta : Double
}

///|
pub fn Gamma::new(k : Double, theta? : Double = 1.0) -> Gamma {
  if k <= 0.0 {
    abort("k (shape parameter) must be positive")
  }
  if theta <= 0.0 {
    abort("theta (scale parameter) must be positive")
  }
  { k, theta }
}

///|
pub impl RVContinuous for Gamma with rvs(
  self,
  size,
  rand? = @random.Rand::new(),
) {
  Array::makei(
    size,
    if self.k >= 1.0 {
      _ => gamma_marsaglia_tsang(self.k, rand) * self.theta
    } else {
      _ => gamma_ahrens_dieter(self.k, rand) * self.theta
    },
  )
}

///|
fn gamma_marsaglia_tsang(k : Double, rand : @random.Rand) -> Double {
  let d = k - 1.0 / 3.0
  let c = 1.0 / (9.0 * d).sqrt()
  let mut ans = 0.0
  while true {
    let x = Normal::new().rvs(1, rand~)[0]
    let v = 1.0 + c * x
    if v > 0.0 {
      let v3 = v * v * v
      let u = rand.double()
      if u < 1.0 - 0.0331 * x * x * x * x {
        ans = d * v3
        break
      }
      if @math.ln(u) < 0.5 * x * x + d * (1.0 - v3 + @math.ln(v3)) {
        ans = d * v3
        break
      }
    }
  }
  ans
}

///|
fn gamma_ahrens_dieter(k : Double, rand : @random.Rand) -> Double {
  let e = @kmath.E
  while true {
    let u = rand.double()
    let v = rand.double()
    if u <= e / (e + k) {
      let x = @math.pow(v, 1.0 / k)
      let y = @math.exp(-x)
      if rand.double() <= y {
        break x
      }
    } else {
      let x = -@math.ln(v)
      let y = @math.pow(x, k - 1.0)
      if rand.double() <= y {
        break x
      }
    }
  } else {
    0.0
  }
}

///|
pub impl RVContinuous for Gamma with pdf(self, x) {
  if x <= 0.0 {
    return 0.0
  }
  let numerator = @math.pow(x, self.k - 1.0) * @math.exp(-x / self.theta)
  let denominator = @math.pow(self.theta, self.k) * @kmath.gamma(self.k)
  numerator / denominator
}

///|
pub impl RVContinuous for Gamma with cdf(self, x) {
  if x <= 0.0 {
    return 0.0
  }
  @sp.gammainc(self.k, x / self.theta)
}

///|
pub impl RVContinuous for Gamma with ppf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  let initial_guess = if self.k > 10.0 {
    let z = Normal::new().ppf(q)
    self.k + z * self.k.sqrt() + (z * z - 1.0) / 3.0
  } else {
    self.k * self.theta
  }
  let ans = @sp.func_inv(
    q,
    self.cdf(_),
    self.pdf(_),
    initial_guess,
    max_iterations=200,
  )
  match ans {
    None => abort("No solution found for ppf")
    Some(val) => val
  }
}

///|
pub impl RVContinuous for Gamma with sf(self, x) {
  if x <= 0.0 {
    return 1.0
  }
  1.0 - self.cdf(x)
}

///|
pub impl RVContinuous for Gamma with isf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  self.ppf(1.0 - q)
}

///|
pub impl RVContinuous for Gamma with logpdf(self, x) {
  @math.ln(self.pdf(x))
}

///|
pub impl RVContinuous for Gamma with logcdf(self, x) {
  @math.ln(self.cdf(x))
}

///|
pub impl RVContinuous for Gamma with logsf(self, x) {
  @math.ln(self.sf(x))
}

///|
pub impl RVContinuous for Gamma with moment(self, n) {
  if n < 0 {
    abort("moment order must be non-negative")
  }
  if n == 0 {
    return 1.0
  }
  @math.pow(self.theta, n.to_double()) *
  @kmath.gamma(self.k + n.to_double()) /
  @kmath.gamma(self.k)
}

///|
pub impl RVContinuous for Gamma with mean(self) {
  self.k * self.theta
}

///|
pub impl RVContinuous for Gamma with median(self) {
  self.ppf(0.5)
}

///|
pub impl RVContinuous for Gamma with varc(self) {
  self.k * self.theta * self.theta
}

///|
pub impl RVContinuous for Gamma with std(self) {
  (self.k * self.theta * self.theta).sqrt()
}

///|
pub impl RVContinuous for Gamma with entropy(self) {
  self.k +
  @math.ln(self.theta) +
  @kmath.ln_gamma(self.k) +
  (1.0 - self.k) * @kmath.digamma(self.k)
}

///|
pub impl RVContinuous for Gamma with interval(self, alpha : Double) {
  if alpha <= 0.0 || alpha >= 1.0 {
    abort("alpha must be in (0, 1)")
  }
  let lower = self.ppf(alpha / 2.0)
  let upper = self.ppf(1.0 - alpha / 2.0)
  (lower, upper)
}

///|
pub impl RVContinuous for Gamma with skew(self) {
  2.0 / self.k.sqrt()
}

///|
pub impl RVContinuous for Gamma with kurt(self) {
  6.0 / self.k
}
