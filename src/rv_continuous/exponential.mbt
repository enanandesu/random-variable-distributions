///|
struct Exponential {
  loc : Double
  scale : Double
}

///|
pub fn Exponential::new(
  loc? : Double = 0.0,
  scale? : Double = 1.0,
) -> Exponential {
  if scale <= 0.0 {
    abort("scale must be positive")
  }
  { loc, scale }
}

///|
fn exp_pow_int(base : Double, exponent : Int) -> Double {
  if exponent == 0 {
    1.0
  } else {
    let mut result = 1.0
    let mut b = base
    let mut e = exponent
    while e > 0 {
      if e % 2 == 1 {
        result *= b
      }
      b *= b
      e /= 2
    }
    result
  }
}

///|
fn exp_moment_term(loc : Double, scale : Double, n : Int) -> Double {
  let mut moment = 0.0
  for k in 0..=n {
    let binom = exp_binomial_coefficient(n, k).to_double()
    let loc_power = exp_pow_int(loc, n - k)
    let scale_power = exp_pow_int(scale, k)
    let factorial = @kmath.gamma(k.to_double() + 1.0)
    moment += binom * loc_power * scale_power * factorial
  }
  moment
}

///|
fn exp_binomial_coefficient(n : Int, k : Int) -> Int {
  if k < 0 || k > n {
    0
  } else if k == 0 || k == n {
    1
  } else {
    let mut result = 1
    let mut i = 1
    let mut j = n
    let kk = if k < n - k { k } else { n - k }
    while i <= kk {
      result = result * j / i
      i += 1
      j -= 1
    }
    result
  }
}

///|
impl RVContinuous for Exponential with rvs(self, size, seed) {
  let rng = match seed {
    Some(val) => @random.Rand::chacha8(seed=val)
    None => @random.Rand::new()
  }
  let samples = Array::new(capacity=size)
  for _i in 0..<size {
    let mut u = rng.double()
    if u <= 0.0 {
      u = 1.0e-12
    } else if u >= 1.0 {
      u = 1.0 - 1.0e-12
    }
    let value = self.loc - self.scale * @math.ln(1.0 - u)
    samples.push(value)
  }
  samples
}

///|
impl RVContinuous for Exponential with pdf(self, x) {
  if x < self.loc {
    0.0
  } else {
    @math.exp(-(x - self.loc) / self.scale) / self.scale
  }
}

///|
impl RVContinuous for Exponential with cdf(self, x) {
  if x < self.loc {
    0.0
  } else {
    1.0 - @math.exp(-(x - self.loc) / self.scale)
  }
}

///|
impl RVContinuous for Exponential with ppf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  self.loc - self.scale * @math.ln(1.0 - q)
}

///|
impl RVContinuous for Exponential with sf(self, x) {
  if x < self.loc {
    1.0
  } else {
    @math.exp(-(x - self.loc) / self.scale)
  }
}

///|
impl RVContinuous for Exponential with isf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  self.ppf(1.0 - q)
}

///|
impl RVContinuous for Exponential with logpdf(self, x) {
  if x < self.loc {
    @double.neg_infinity
  } else {
    -(x - self.loc) / self.scale - @math.ln(self.scale)
  }
}

///|
impl RVContinuous for Exponential with logcdf(self, x) {
  let value = self.cdf(x)
  if value <= 0.0 {
    @double.neg_infinity
  } else {
    @math.ln(value)
  }
}

///|
impl RVContinuous for Exponential with logsf(self, x) {
  let value = self.sf(x)
  if value <= 0.0 {
    @double.neg_infinity
  } else {
    @math.ln(value)
  }
}

///|
impl RVContinuous for Exponential with moment(self, n) {
  if n < 0 {
    abort("moment order must be non-negative")
  }
  if n == 0 {
    1.0
  } else {
    exp_moment_term(self.loc, self.scale, n)
  }
}

///|
impl RVContinuous for Exponential with entropy(self) {
  1.0 + @math.ln(self.scale)
}

///|
impl RVContinuous for Exponential with median(self) {
  self.loc + self.scale * @math.ln(2.0)
}

///|
impl RVContinuous for Exponential with mean(self) {
  self.loc + self.scale
}

///|
impl RVContinuous for Exponential with varc(self) {
  self.scale * self.scale
}

///|
impl RVContinuous for Exponential with std(self) {
  self.scale
}

///|
impl RVContinuous for Exponential with interval(self, alpha : Double) {
  if alpha <= 0.0 || alpha >= 1.0 {
    abort("alpha must be in (0, 1)")
  }
  let lower = self.ppf(alpha / 2.0)
  let upper = self.ppf(1.0 - alpha / 2.0)
  (lower, upper)
}

///|
impl RVContinuous for Exponential with skew(_) {
  2.0
}

///|
impl RVContinuous for Exponential with kurt(_) {
  6.0
}
