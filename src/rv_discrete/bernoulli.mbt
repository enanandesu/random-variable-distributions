///|
struct Bernoulli {
  p : Double
  loc : Double
}

///|
pub fn Bernoulli::new(p : Double, loc? : Double = 0.0) -> Bernoulli {
  { p, loc }
}

///|
impl RVDiscrete for Bernoulli with rvs(self, size, seed) {
  let rng = match seed {
    Some(s) => @random.Rand::chacha8(seed=s)
    None => @random.Rand::new()
  }
  let result = Array::new(capacity=size)
  for i = 0; i < size; i = i + 1 {
    let u = rng.double()
    if u < self.p {
      result.push(self.loc + 1.0)
    } else {
      result.push(self.loc)
    }
  }
  result
}

///|
impl RVDiscrete for Bernoulli with pmf(self, k) {
  if k == self.loc {
    1.0 - self.p
  } else if k == self.loc + 1.0 {
    self.p
  } else {
    0.0
  }
}

///|
impl RVDiscrete for Bernoulli with cdf(self, x) {
  if x < self.loc {
    0
  } else if x < self.loc + 1.0 {
    1.0 - self.p
  } else {
    1.0
  }
}

///|
impl RVDiscrete for Bernoulli with ppf(self, q) {
  if q < 0.0 || q > 1.0 {
    abort("q must be in [0, 1]")
  }
  if q == 0.0 {
    self.loc - 1.0
  } else if q <= 1.0 - self.p {
    self.loc
  } else {
    self.loc + 1.0
  }
}

///|
impl RVDiscrete for Bernoulli with sf(self, x) {
  1.0 - self.cdf(x)
}

///|
impl RVDiscrete for Bernoulli with isf(self, q) {
  if q < 0.0 || q > 1.0 {
    abort("q must be in [0, 1]")
  }
  self.ppf(1.0 - q)
}

///|
impl RVDiscrete for Bernoulli with logpmf(self, x) {
  @math.ln(self.pmf(x))
}

///|
impl RVDiscrete for Bernoulli with logcdf(self, x) {
  @math.ln(self.cdf(x))
}

///|
impl RVDiscrete for Bernoulli with logsf(self, x) {
  @math.ln(self.sf(x))
}

///|
impl RVDiscrete for Bernoulli with mean(self) {
  self.loc + self.p
}

///|
impl RVDiscrete for Bernoulli with median(self) {
  self.loc
}

///|
impl RVDiscrete for Bernoulli with varc(self) {
  self.p - self.p * self.p
}

///|
impl RVDiscrete for Bernoulli with std(self) {
  (self.p - self.p * self.p).sqrt()
}

///|
impl RVDiscrete for Bernoulli with entropy(self) {
  let q = 1.0 - self.p
  -q * @math.ln(q) - self.p * @math.ln(self.p)
}

///|
impl RVDiscrete for Bernoulli with interval(self, alpha : Double) {
  if alpha < 0.0 || alpha > 1.0 {
    abort("alpha must be in [0, 1]")
  }
  let lower = self.ppf(alpha / 2.0)
  let upper = self.ppf(1.0 - alpha / 2.0)
  (lower, upper)
}

///|
impl RVDiscrete for Bernoulli with skew(self) {
  let q = 1.0 - self.p
  let p = self.p
  (q - p) / (p * q).sqrt()
}

///|
impl RVDiscrete for Bernoulli with kurt(self) {
  let p = self.p
  (6.0 * p * p - 6.0 * p + 1.0) / (p * (1.0 - p))
}
