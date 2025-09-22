///|
struct F {
  df1 : Double
  df2 : Double
  loc : Double
  scale : Double
}

///|
pub fn F::new(
  df1 : Double,
  df2 : Double,
  loc? : Double = 0.0,
  scale? : Double = 1.0,
) -> F {
  if df1 <= 0.0 {
    abort("df1 must be positive")
  }
  if df2 <= 0.0 {
    abort("df2 must be positive")
  }
  if scale <= 0.0 {
    abort("scale must be positive")
  }
  { df1, df2, loc, scale }
}

///|
fn fisher_f_pow_int(base : Double, exponent : Int) -> Double {
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
fn fisher_f_binomial_coefficient(n : Int, k : Int) -> Int {
  if k < 0 || k > n {
    0
  } else if k == 0 || k == n {
    1
  } else {
    let mut result = 1
    let mut i = 1
    let mut j = n
    let limit = if k < n - k { k } else { n - k }
    while i <= limit {
      result = result * j / i
      i += 1
      j -= 1
    }
    result
  }
}

///|
fn fisher_f_standard_moment(df1 : Double, df2 : Double, order : Int) -> Double {
  if df2 <= 2.0 * order.to_double() {
    @double.infinity
  } else {
    let a = df1 / 2.0
    let b = df2 / 2.0
    let k = order.to_double()
    let ratio = @sp.beta(a + k, b - k) / @sp.beta(a, b)
    @math.pow(df2 / df1, k) * ratio
  }
}

///|
fn fisher_f_entropy(df1 : Double, df2 : Double, scale : Double) -> Double {
  let a = df1 / 2.0
  let b = df2 / 2.0
  let term = @sp.betaln(a, b) -
    (a - 1.0) * @kmath.digamma(a) -
    (b + 1.0) * @kmath.digamma(b) +
    (a + b) * @kmath.digamma(a + b) +
    @math.ln(df2 / df1) +
    @math.ln(scale)
  term
}

///|

///|
impl RVContinuous for F with rvs(self, size, seed) {
  let rng = match seed {
    Some(val) => @random.Rand::chacha8(seed=val)
    None => @random.Rand::new()
  }
  let samples = Array::new(capacity=size)
  let a = self.df1 / 2.0
  let b = self.df2 / 2.0
  for _i in 0..<size {
    let mut u = rng.double()
    if u <= 0.0 {
      u = 1.0e-12
    } else if u >= 1.0 {
      u = 1.0 - 1.0e-12
    }
    let mut z = @sp.betaincinv(a, b, u)
    if z <= 0.0 {
      z = 1.0e-12
    } else if z >= 1.0 {
      z = 1.0 - 1.0e-12
    }
    let y = (self.df2 * z) / (self.df1 * (1.0 - z))
    samples.push(self.loc + self.scale * y)
  }
  samples
}

///|
impl RVContinuous for F with pdf(self, x) {
  let y = (x - self.loc) / self.scale
  if y <= 0.0 {
    0.0
  } else {
    @math.exp(self.logpdf(x))
  }
}

///|
impl RVContinuous for F with cdf(self, x) {
  let y = (x - self.loc) / self.scale
  if y <= 0.0 {
    0.0
  } else {
    let numerator = self.df1 * y
    let denominator = numerator + self.df2
    let z = numerator / denominator
    @sp.betainc(self.df1 / 2.0, self.df2 / 2.0, z)
  }
}

///|
impl RVContinuous for F with ppf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  let a = self.df1 / 2.0
  let b = self.df2 / 2.0
  let mut z = @sp.betaincinv(a, b, q)
  if z <= 0.0 {
    z = 1.0e-12
  } else if z >= 1.0 {
    z = 1.0 - 1.0e-12
  }
  let y = (self.df2 * z) / (self.df1 * (1.0 - z))
  self.loc + self.scale * y
}

///|
impl RVContinuous for F with sf(self, x) {
  1.0 - self.cdf(x)
}

///|
impl RVContinuous for F with isf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  self.ppf(1.0 - q)
}

///|
impl RVContinuous for F with logpdf(self, x) {
  let y = (x - self.loc) / self.scale
  if y <= 0.0 {
    @double.neg_infinity
  } else {
    let ratio = self.df1 / self.df2
    let a = self.df1 / 2.0
    let b = self.df2 / 2.0
    a * @math.ln(ratio) +
    (a - 1.0) * @math.ln(y) -
    @sp.betaln(a, b) -
    (a + b) * @math.ln(1.0 + ratio * y) -
    @math.ln(self.scale)
  }
}

///|
impl RVContinuous for F with logcdf(self, x) {
  let value = self.cdf(x)
  if value <= 0.0 {
    @double.neg_infinity
  } else {
    @math.ln(value)
  }
}

///|
impl RVContinuous for F with logsf(self, x) {
  let value = self.sf(x)
  if value <= 0.0 {
    @double.neg_infinity
  } else {
    @math.ln(value)
  }
}

///|
impl RVContinuous for F with moment(self, n) {
  if n < 0 {
    abort("moment order must be non-negative")
  }
  if n == 0 {
    1.0
  } else {
    let mut total = 0.0
    for k in 0..=n {
      let binom = fisher_f_binomial_coefficient(n, k).to_double()
      let loc_power = fisher_f_pow_int(self.loc, n - k)
      let scale_power = fisher_f_pow_int(self.scale, k)
      let std_moment = fisher_f_standard_moment(self.df1, self.df2, k)
      if std_moment == @double.infinity {
        return @double.infinity
      }
      total += binom * loc_power * scale_power * std_moment
    }
    total
  }
}

///|
impl RVContinuous for F with entropy(self) {
  fisher_f_entropy(self.df1, self.df2, self.scale)
}

///|
impl RVContinuous for F with median(self) {
  self.ppf(0.5)
}

///|
impl RVContinuous for F with mean(self) {
  if self.df2 <= 2.0 {
    @double.infinity
  } else {
    self.loc + self.scale * self.df2 / (self.df2 - 2.0)
  }
}

///|
impl RVContinuous for F with varc(self) {
  if self.df2 <= 4.0 {
    @double.infinity
  } else {
    let numerator = 2.0 * self.df2 * self.df2 * (self.df1 + self.df2 - 2.0)
    let denominator = self.df1 *
      (self.df2 - 2.0) *
      (self.df2 - 2.0) *
      (self.df2 - 4.0)
    self.scale * self.scale * numerator / denominator
  }
}

///|
impl RVContinuous for F with std(self) {
  self.varc().sqrt()
}

///|
impl RVContinuous for F with interval(self, alpha : Double) {
  if alpha <= 0.0 || alpha >= 1.0 {
    abort("alpha must be in (0, 1)")
  }
  let lower = self.ppf(alpha / 2.0)
  let upper = self.ppf(1.0 - alpha / 2.0)
  (lower, upper)
}

///|
impl RVContinuous for F with skew(self) {
  if self.df2 <= 6.0 {
    @double.infinity
  } else {
    let numerator = (2.0 * self.df1 + self.df2 - 2.0) *
      (8.0 * (self.df2 - 4.0)).sqrt()
    let denominator = (self.df2 - 6.0) *
      (self.df1 * (self.df1 + self.df2 - 2.0)).sqrt()
    numerator / denominator
  }
}

///|
impl RVContinuous for F with kurt(self) {
  if self.df2 <= 8.0 {
    @double.infinity
  } else {
    let numerator = 12.0 *
      (
        self.df1 * (5.0 * self.df2 - 22.0) * (self.df1 + self.df2 - 2.0) +
        (self.df2 - 4.0) * (self.df2 - 2.0) * (self.df2 - 2.0)
      )
    let denominator = self.df1 *
      (self.df2 - 6.0) *
      (self.df2 - 8.0) *
      (self.df1 + self.df2 - 2.0)
    numerator / denominator
  }
}
