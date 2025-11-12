///|
struct BetaBinom {
  n : Int
  a : Double
  b : Double
  loc : Double
}

///|
pub fn BetaBinom::new(
  n : Int,
  a : Double,
  b : Double,
  loc? : Double = 0.0,
) -> BetaBinom {
  if n < 0 {
    abort("n must be non-negative")
  }
  if a <= 0.0 {
    abort("a must be positive")
  }
  if b <= 0.0 {
    abort("b must be positive")
  }
  { n, a, b, loc }
}

///|
pub impl RVDiscrete for BetaBinom with rvs(
  self,
  size,
  rand? = @random.Rand::new(),
) {
  let result = Array::new(capacity=size)
  for i = 0; i < size; i = i + 1 {
    let p = generate_beta(self.a, self.b, rand)
    let mut count = 0.0
    for j = 0; j < self.n; j = j + 1 {
      let u = rand.double()
      if u < p {
        count += 1.0
      }
    }
    result.push(self.loc + count)
  }
  result
}

///|
fn generate_beta(a : Double, b : Double, rng : @random.Rand) -> Double {
  let mut x = 0.0
  let mut y = 0.0
  let mut sum = 0.0
  while {
          x = rng.double().pow(1.0 / a)
          y = rng.double().pow(1.0 / b)
          sum = x + y
          sum > 1.0
        } {
    continue
  }
  x / sum
}

///|
pub impl RVDiscrete for BetaBinom with pmf(self, k) {
  let diff = k - self.loc
  if diff != diff.round() {
    0.0
  } else {
    let k_int = diff.to_int()
    if k_int < 0 || k_int > self.n {
      0.0
    } else {
      let comb = combination(self.n, k_int).to_double()
      comb *
      @sp.beta(diff + self.a, self.n.to_double() - diff + self.b) /
      @sp.beta(self.a, self.b)
    }
  }
}

///|
pub impl RVDiscrete for BetaBinom with cdf(self, x) {
  let k_max = (x - self.loc).floor().to_int()
  if k_max < 0 {
    0.0
  } else if k_max >= self.n {
    1.0
  } else {
    let mut sum = 0.0
    for k = 0; k <= k_max; k = k + 1 {
      sum += self.pmf(self.loc + k.to_double())
    }
    sum
  }
}

///|
pub impl RVDiscrete for BetaBinom with ppf(self, q) {
  if q < 0.0 || q > 1.0 {
    abort("q must be in [0, 1]")
  }
  if q == 0.0 {
    self.loc - 1.0
  } else if q == 1.0 {
    self.loc + self.n.to_double()
  } else {
    let mut cumsum = 0.0
    for k = 0; k <= self.n; k = k + 1 {
      cumsum += self.pmf(self.loc + k.to_double())
      if cumsum >= q {
        return self.loc + k.to_double()
      }
    }
    self.loc + self.n.to_double()
  }
}

///|
pub impl RVDiscrete for BetaBinom with sf(self, x) {
  1.0 - self.cdf(x)
}

///|
pub impl RVDiscrete for BetaBinom with isf(self, q) {
  if q < 0.0 || q > 1.0 {
    abort("q must be in [0, 1]")
  }
  self.ppf(1.0 - q)
}

///|
pub impl RVDiscrete for BetaBinom with logpmf(self, x) {
  @math.ln(self.pmf(x))
}

///|
pub impl RVDiscrete for BetaBinom with logcdf(self, x) {
  @math.ln(self.cdf(x))
}

///|
pub impl RVDiscrete for BetaBinom with logsf(self, x) {
  @math.ln(self.sf(x))
}

///|
pub impl RVDiscrete for BetaBinom with mean(self) {
  self.loc + self.n.to_double() * self.a / (self.a + self.b)
}

///|
pub impl RVDiscrete for BetaBinom with median(self) {
  self.ppf(0.5)
}

///|
pub impl RVDiscrete for BetaBinom with varc(self) {
  self.n.to_double() *
  self.a *
  self.b *
  (self.n.to_double() + self.a + self.b) /
  ((self.a + self.b) * (self.a + self.b) * (self.a + self.b + 1.0))
}

///|
pub impl RVDiscrete for BetaBinom with std(self) {
  self.varc().sqrt()
}

///|
pub impl RVDiscrete for BetaBinom with entropy(self) {
  let mut entropy = 0.0
  for k in 0..=self.n {
    let prob = self.pmf(self.loc + k.to_double())
    if prob > 0.0 {
      entropy -= prob * @math.ln(prob)
    }
  }
  entropy
}

///|
pub impl RVDiscrete for BetaBinom with interval(self, alpha : Double) {
  if alpha < 0.0 || alpha > 1.0 {
    abort("alpha must be in [0, 1]")
  }
  let lower = self.ppf(alpha / 2.0)
  let upper = self.ppf(1.0 - alpha / 2.0)
  (lower, upper)
}

///|
pub impl RVDiscrete for BetaBinom with skew(self) {
  let a = self.a
  let b = self.b
  let n = self.n.to_double()
  let alpha = a + b
  let numerator = (alpha + 2 * n) * (b - a) * (1.0 + alpha).sqrt()
  let denominator = (alpha + 2.0) * (n * a * b * (n + alpha)).sqrt()
  numerator / denominator
}

///|
pub impl RVDiscrete for BetaBinom with kurt(self) {
  let a = self.a
  let b = self.b
  let n = self.n.to_double()
  let e_p = a / (a + b)
  let e_q = 1.0 - e_p
  let mut g2 = a + b
  g2 *= a + b - 1.0 + 6.0 * n
  g2 += 3.0 * a * b * (n - 2.0)
  g2 += 6.0 * n * n
  g2 -= 3.0 * e_p * b * n * (6.0 - n)
  g2 -= 18.0 * e_p * e_q * n * n
  g2 *= (a + b) * (a + b) * (1.0 + a + b)
  g2 /= n * a * b * (a + b + 2.0) * (a + b + 3.0) * (a + b + n)
  g2 -= 3.0
  g2
}
