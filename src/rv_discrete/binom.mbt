///|
struct Binom {
  n : Int
  p : Double
  loc : Double
}

///|
pub fn Binom::new(n : Int, p : Double, loc~ : Double = 0.0) -> Binom {
  { n, p, loc }
}

///|
impl RVDiscrete for Binom with rvs(self, size, seed) {
  let rng = match seed {
    Some(s) => @random.Rand::chacha8(seed=s)
    None => @random.Rand::new()
  }
  let result = Array::new(capacity=size)
  for i = 0; i < size; i = i + 1 {
    let mut count = 0.0
    for j = 0; j < self.n; j = j + 1 {
      let u = rng.double()
      if u < self.p {
        count += 1.0
      }
    }
    result.push(self.loc + count)
  }
  result
}

///|
fn combination(n : Int, k : Int) -> Int {
  if k < 0 || k > n {
    0
  } else if k == 0 || k == n {
    1
  } else if k == 1 {
    n
  } else {
    let k_opt = if k > n - k { n - k } else { k }
    let mut result = 1
    for i = 0; i < k_opt; i = i + 1 {
      result = result * (n - i) / (i + 1)
    }
    result
  }
}

///|
impl RVDiscrete for Binom with pmf(self, k) {
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
      self.p.pow(k_int.to_double()) *
      (1.0 - self.p).pow((self.n - k_int).to_double())
    }
  }
}

///|
impl RVDiscrete for Binom with cdf(self, x) {
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
impl RVDiscrete for Binom with ppf(self, q) {
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
impl RVDiscrete for Binom with sf(self, x) {
  1.0 - self.cdf(x)
}

///|
impl RVDiscrete for Binom with isf(self, q) {
  if q < 0.0 || q > 1.0 {
    abort("q must be in [0, 1]")
  }
  self.ppf(1.0 - q)
}

///|
impl RVDiscrete for Binom with logpmf(self, x) {
  @math.ln(self.pmf(x))
}

///|
impl RVDiscrete for Binom with logcdf(self, x) {
  @math.ln(self.cdf(x))
}

///|
impl RVDiscrete for Binom with logsf(self, x) {
  @math.ln(self.sf(x))
}

///|
impl RVDiscrete for Binom with moment(self, n) {
  if n < 0 {
    abort("moment order must be non-negative")
  }
  let mut sum = 0.0
  for k in 0..=self.n {
    let x = self.loc + k.to_double()
    sum += x.pow(n.to_double()) * self.pmf(x)
  }
  sum
}

///|
impl RVDiscrete for Binom with mean(self) {
  self.loc + self.n.to_double() * self.p
}

///|
impl RVDiscrete for Binom with median(self) {
  self.loc + (self.n.to_double() * self.p).round()
}

///|
impl RVDiscrete for Binom with varc(self) {
  self.n.to_double() * self.p * (1.0 - self.p)
}

///|
impl RVDiscrete for Binom with std(self) {
  (self.n.to_double() * self.p * (1.0 - self.p)).sqrt()
}

///|
impl RVDiscrete for Binom with entropy(self) {
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
impl RVDiscrete for Binom with interval(self, alpha : Double) {
  if alpha < 0.0 || alpha > 1.0 {
    abort("alpha must be in [0, 1]")
  }
  let lower = self.ppf(alpha / 2.0)
  let upper = self.ppf(1.0 - alpha / 2.0)
  (lower, upper)
}

///|
impl RVDiscrete for Binom with skew(self) {
  let p = self.p
  let q = 1.0 - p
  let n = self.n.to_double()
  (q - p) / (n * p * q).sqrt()
}

///|
impl RVDiscrete for Binom with kurt(self) {
  let p = self.p
  let q = 1.0 - p
  let n = self.n.to_double()
  (1.0 - 6.0 * p * q) / (n * p * q)
}
