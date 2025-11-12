///|
struct DLaplace {
  a : Double
  loc : Double
}

///|
pub fn DLaplace::new(a : Double, loc? : Double = 0.0) -> DLaplace {
  if a <= 0.0 {
    abort("a must be positive")
  }
  { a, loc }
}

///|
fn dlaplace_q(a : Double) -> Double {
  @math.exp(-a)
}

///|
fn dlaplace_norm(a : Double) -> Double {
  @math.tanh(a / 2.0)
}

///|
fn dlaplace_one_minus_q(a : Double) -> Double {
  -@math.expm1(-a)
}

///|
pub impl RVDiscrete for DLaplace with rvs(self, size, seed) {
  let rng = match seed {
    Some(s) => @random.Rand::chacha8(seed=s)
    None => @random.Rand::new()
  }
  let c = dlaplace_norm(self.a)
  let result = Array::new(capacity=size)
  for _i in 0..<size {
    let u = rng.double()
    if u < c {
      result.push(self.loc)
    } else {
      let sign = if rng.double() < 0.5 { -1.0 } else { 1.0 }
      let v = rng.double()
      let mut geom = (@math.ln(1.0 - v) / -self.a).floor().to_int() + 1
      if geom < 1 {
        geom = 1
      }
      let offset = geom.to_double() * sign
      result.push(self.loc + offset)
    }
  }
  result
}

///|
pub impl RVDiscrete for DLaplace with pmf(self, k) {
  let diff = k - self.loc
  if diff != diff.round() {
    0.0
  } else {
    let abs_diff = diff.abs()
    let norm = dlaplace_norm(self.a)
    norm * @math.exp(-self.a * abs_diff)
  }
}

///|
pub impl RVDiscrete for DLaplace with logpmf(self, k) {
  let diff = k - self.loc
  if diff != diff.round() {
    @double.neg_infinity
  } else {
    @math.ln(dlaplace_norm(self.a)) - self.a * diff.abs()
  }
}

///|
pub impl RVDiscrete for DLaplace with cdf(self, x) {
  let q = dlaplace_q(self.a)
  let one_minus_q = dlaplace_one_minus_q(self.a)
  let norm = dlaplace_norm(self.a)
  let diff = (x - self.loc).floor().to_int()
  if diff < 0 {
    let power = (-diff).to_double()
    norm * q.pow(power) / one_minus_q
  } else {
    let power = diff.to_double() + 1.0
    1.0 - norm * q.pow(power) / one_minus_q
  }
}

///|
pub impl RVDiscrete for DLaplace with ppf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  let decay = dlaplace_q(self.a)
  let norm = dlaplace_norm(self.a)
  let one_minus_q = dlaplace_one_minus_q(self.a)
  let total_neg = norm * decay / one_minus_q
  if q <= total_neg {
    let ratio = q / total_neg
    let mut m = (@math.ln(ratio) / @math.ln(decay)).floor().to_int() + 1
    if m < 1 {
      m = 1
    }
    self.loc - m.to_double()
  } else {
    let pivot = total_neg + norm
    if q <= pivot {
      self.loc
    } else {
      let threshold = (1.0 - q) * one_minus_q / norm
      let mut m = (@math.ln(threshold) / @math.ln(decay)).ceil().to_int() - 1
      if m < 0 {
        m = 0
      }
      self.loc + m.to_double()
    }
  }
}

///|
pub impl RVDiscrete for DLaplace with sf(self, x) {
  1.0 - self.cdf(x)
}

///|
pub impl RVDiscrete for DLaplace with isf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  self.ppf(1.0 - q)
}

///|
pub impl RVDiscrete for DLaplace with logcdf(self, x) {
  let value = self.cdf(x)
  if value <= 0.0 {
    @double.neg_infinity
  } else {
    @math.ln(value)
  }
}

///|
pub impl RVDiscrete for DLaplace with logsf(self, x) {
  let value = self.sf(x)
  if value <= 0.0 {
    @double.neg_infinity
  } else {
    @math.ln(value)
  }
}

///|
pub impl RVDiscrete for DLaplace with entropy(self) {
  self.a / @math.sinh(self.a) - @math.ln(@math.tanh(self.a / 2.0))
}

///|
pub impl RVDiscrete for DLaplace with median(self) {
  self.loc
}

///|
pub impl RVDiscrete for DLaplace with mean(self) {
  self.loc
}

///|
pub impl RVDiscrete for DLaplace with varc(self) {
  let decay = dlaplace_q(self.a)
  let one_minus_q = dlaplace_one_minus_q(self.a)
  2.0 * decay / (one_minus_q * one_minus_q)
}

///|
pub impl RVDiscrete for DLaplace with std(self) {
  self.varc().sqrt()
}

///|
pub impl RVDiscrete for DLaplace with interval(self, alpha : Double) {
  if alpha < 0.0 || alpha > 1.0 {
    abort("alpha must be in [0, 1]")
  }
  let lower = self.ppf(alpha / 2.0)
  let upper = self.ppf(1.0 - alpha / 2.0)
  (lower, upper)
}

///|
pub impl RVDiscrete for DLaplace with skew(_) {
  0.0
}

///|
pub impl RVDiscrete for DLaplace with kurt(self) {
  let decay = dlaplace_q(self.a)
  let numerator = 1.0 + 10.0 * decay + decay * decay
  numerator / (2.0 * decay) - 3.0
}
