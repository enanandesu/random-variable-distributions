///|
struct BetaNegativeBinom {
  n : Int
  a : Double
  b : Double
  loc : Double
}

///|
pub fn BetaNegativeBinom::new(
  n : Int,
  a : Double,
  b : Double,
  loc? : Double = 0.0,
) -> BetaNegativeBinom {
  { n, a, b, loc }
}

///|
impl RVDiscrete for BetaNegativeBinom with rvs(self, size, seed) {
  let rng = match seed {
    Some(s) => @random.Rand::chacha8(seed=s)
    None => @random.Rand::new()
  }
  let result = Array::new(capacity=size)
  for i in 0..<size {
    let p = generate_beta(self.a, self.b, rng)
    let mut failures = 0.0
    let mut successes = 0
    while successes < self.n {
      let u = rng.double()
      if u < p {
        successes += 1
      } else {
        failures += 1.0
      }
    }
    result.push(self.loc + failures)
  }
  result
}

///|
impl RVDiscrete for BetaNegativeBinom with pmf(self, k) {
  @math.exp(self.logpmf(k))
}

///|
impl RVDiscrete for BetaNegativeBinom with cdf(self, x) {
  let k_max = (x - self.loc).floor().to_int()
  if k_max < 0 {
    0.0
  } else {
    let mut sum = 0.0
    for k = 0; k <= k_max; k = k + 1 {
      sum += self.pmf(self.loc + k.to_double())
    }
    sum
  }
}

///|
impl RVDiscrete for BetaNegativeBinom with ppf(self, q) {
  if q < 0.0 || q >= 1.0 {
    abort("q must be in [0, 1)")
  }
  if q == 0.0 {
    self.loc - 1.0
  } else {
    let mut cumsum = 0.0
    let mut k = 0
    while cumsum < q {
      cumsum += self.pmf(self.loc + k.to_double())
      if cumsum >= q {
        return self.loc + k.to_double()
      }
      k += 1
    }
    self.loc + k.to_double()
  }
}

///|
impl RVDiscrete for BetaNegativeBinom with sf(self, x) {
  1.0 - self.cdf(x)
}

///|
impl RVDiscrete for BetaNegativeBinom with isf(self, q) {
  if q <= 0.0 || q > 1.0 {
    abort("q must be in (0, 1]")
  }
  self.ppf(1.0 - q)
}

///|
fn log_combination(n : Int, k : Int) -> Double {
  if k < 0 || k > n {
    @double.neg_infinity
  } else if k == 0 || k == n {
    0.0
  } else {
    @kmath.lgamma((n + 1).to_double()) -
    @kmath.lgamma((k + 1).to_double()) -
    @kmath.lgamma((n - k + 1).to_double())
  }
}

///|
impl RVDiscrete for BetaNegativeBinom with logpmf(self, k) {
  let diff = k - self.loc
  if diff != diff.round() || diff < 0.0 {
    @double.neg_infinity
  } else {
    let k_int = diff.to_int()
    let lgcomb = log_combination(k_int + self.n - 1, k_int)
    lgcomb +
    @sp.betaln(self.n.to_double() + self.a, diff + self.b) -
    @sp.betaln(self.a, self.b)
  }
}

///|
impl RVDiscrete for BetaNegativeBinom with logcdf(self, x) {
  @math.ln(self.cdf(x))
}

///|
impl RVDiscrete for BetaNegativeBinom with logsf(self, x) {
  @math.ln(self.sf(x))
}

///|
impl RVDiscrete for BetaNegativeBinom with mean(self) {
  if self.a <= 1.0 {
    @double.infinity
  } else {
    self.loc + self.n.to_double() * self.b / (self.a - 1.0)
  }
}

///|
impl RVDiscrete for BetaNegativeBinom with median(self) {
  self.ppf(0.5)
}

///|
impl RVDiscrete for BetaNegativeBinom with varc(self) {
  if self.a <= 2.0 {
    @double.infinity
  } else {
    let a = self.a
    let b = self.b
    let r = self.n.to_double()
    r * b * (r + a - 1.0) * (b + a - 1.0) / (a - 2.0) / (a - 1.0) / (a - 1.0)
  }
}

///|
impl RVDiscrete for BetaNegativeBinom with std(self) {
  self.varc().sqrt()
}

///|
impl RVDiscrete for BetaNegativeBinom with entropy(self) {
  let mut entropy = 0.0
  let mut cumsum = 0.0
  let max_k = 1_000
  for k in 0..=max_k {
    let x = self.loc + k.to_double()
    let lgprob = self.logpmf(x)
    let prob = @math.exp(lgprob)
    if prob == 0.0 {
      break
    }
    entropy -= prob * lgprob
    cumsum += prob
    if 1.0 - cumsum < 1.0e-12 {
      break
    }
  }
  entropy
}

///|
impl RVDiscrete for BetaNegativeBinom with interval(self, alpha : Double) {
  if alpha < 0.0 || alpha > 1.0 {
    abort("alpha must be in [0, 1]")
  }
  let lower = self.ppf(alpha / 2.0)
  let upper = self.ppf(1.0 - alpha / 2.0)
  (lower, upper)
}

///|
impl RVDiscrete for BetaNegativeBinom with skew(self) {
  if self.a <= 3.0 {
    @double.infinity
  } else {
    let a = self.a
    let b = self.b
    let r = self.n.to_double()
    let numerator = (2.0 * r + a - 1.0) * (2.0 * b + a - 1.0) * (a - 2.0).sqrt()
    let denominator = (a - 3.0) * (r * b * (r + a - 1.0) * (b + a - 1.0)).sqrt()
    numerator / denominator
  }
}

///|
impl RVDiscrete for BetaNegativeBinom with kurt(self) {
  let n = self.n.to_double()
  let a = self.a
  let b = self.b
  let term = a - 2.0
  let term_2 = (a - 1.0) *
    (a - 1.0) *
    (a * a + a * (6.0 * b - 1.0) + 6.0 * (b - 1.0) * b) +
    3.0 *
    n *
    n *
    (
      (a + 5.0) * b * b +
      (a + 5.0) * (a - 1.0) * b +
      2.0 * (a - 1.0) * (a - 1.0)
    ) +
    3.0 *
    (a - 1.0) *
    n *
    (
      (a + 5.0) * b * b +
      (a + 5.0) * (a - 1.0) * b +
      2.0 * (a - 1.0) * (a - 1.0)
    )
  let denominator = (a - 4.0) *
    (a - 3.0) *
    b *
    n *
    (a + b - 1.0) *
    (a + n - 1.0)
  term * term_2 / denominator - 3.0
}
