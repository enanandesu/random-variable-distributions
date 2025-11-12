///|
struct Boltzmann {
  lambda : Double
  n : Int
  loc : Double
}

///|
pub fn Boltzmann::new(
  lambda : Double,
  n : Int,
  loc? : Double = 0.0,
) -> Boltzmann {
  if lambda <= 0.0 {
    abort("lambda must be positive")
  }
  if n <= 0 {
    abort("n must be positive")
  }
  { lambda, n, loc }
}

///|
fn boltzmann_q(lambda : Double) -> Double {
  @math.exp(-lambda)
}

///|
fn boltzmann_norm(lambda : Double, n : Int) -> Double {
  let numerator = -@math.expm1(-lambda)
  let denominator = -@math.expm1(-lambda * n.to_double())
  numerator / denominator
}

///|
fn boltzmann_moment_terms(
  lambda : Double,
  n : Double,
) -> (Double, Double, Double, Double, Double, Double) {
  if n <= 0.0 {
    (0.0, 0.0, 1.0, 0.0, 0.0, 0.0)
  } else {
    let z = boltzmann_q(lambda)
    if z == 0.0 {
      (0.0, 0.0, 1.0, 0.0, 0.0, 0.0)
    } else {
      let z_n = @math.exp(-lambda * n)
      let one_minus_z = -@math.expm1(-lambda)
      let one_minus_z_n = -@math.expm1(-lambda * n)
      let trm = one_minus_z_n / one_minus_z
      let mu = z / one_minus_z - n * z_n / one_minus_z_n
      let denom_z = one_minus_z * one_minus_z
      let denom_z_n = one_minus_z_n * one_minus_z_n
      let mut variance = z / denom_z - n * n * z_n / denom_z_n
      if variance < 0.0 && variance.abs() < 1.0e-12 {
        variance = 0.0
      }
      let trm_sq = trm * trm
      let mut trm2 = z * trm_sq - n * n * z_n
      if trm2 < 0.0 && trm2.abs() < 1.0e-12 {
        trm2 = 0.0
      }
      (mu, variance, trm, trm2, z, z_n)
    }
  }
}

///|
fn boltzmann_index(value : Double, loc : Double, n : Int) -> Int? {
  let diff = value - loc
  if diff != diff.round() {
    None
  } else {
    let idx = diff.to_int()
    if idx < 0 || idx >= n {
      None
    } else {
      Some(idx)
    }
  }
}

///|
pub impl RVDiscrete for Boltzmann with rvs(self, size, seed) {
  let rng = match seed {
    Some(s) => @random.Rand::chacha8(seed=s)
    None => @random.Rand::new()
  }
  let result = Array::new(capacity=size)
  let norm = boltzmann_norm(self.lambda, self.n)
  let q = boltzmann_q(self.lambda)
  let cdf = Array::new(capacity=self.n)
  let mut cumulative = 0.0
  let mut weight = 1.0
  for idx = 0; idx < self.n; idx = idx + 1 {
    let prob = norm * weight
    cumulative += prob
    cdf.push(cumulative)
    weight *= q
  }
  for _i in 0..<size {
    let u = rng.double()
    let mut value = self.loc + (self.n - 1).to_double()
    for idx = 0; idx < cdf.length(); idx = idx + 1 {
      if u <= cdf[idx] {
        value = self.loc + idx.to_double()
        break
      }
    }
    result.push(value)
  }
  result
}

///|
pub impl RVDiscrete for Boltzmann with pmf(self, k) {
  match boltzmann_index(k, self.loc, self.n) {
    Some(index) => {
      let norm = boltzmann_norm(self.lambda, self.n)
      norm * @math.exp(-self.lambda * index.to_double())
    }
    None => 0.0
  }
}

///|
pub impl RVDiscrete for Boltzmann with cdf(self, x) {
  let k_max = (x - self.loc).floor().to_int()
  if k_max < 0 {
    0.0
  } else if k_max >= self.n - 1 {
    1.0
  } else {
    let norm = boltzmann_norm(self.lambda, self.n)
    let ratio = boltzmann_q(self.lambda)
    let mut sum = 0.0
    let mut weight = 1.0
    for idx = 0; idx <= k_max; idx = idx + 1 {
      sum += norm * weight
      weight *= ratio
    }
    sum
  }
}

///|
pub impl RVDiscrete for Boltzmann with ppf(self, q) {
  if q < 0.0 || q > 1.0 {
    abort("q must be in [0, 1]")
  }
  if q == 0.0 {
    self.loc - 1.0
  } else if q == 1.0 {
    self.loc + (self.n - 1).to_double()
  } else {
    let norm = boltzmann_norm(self.lambda, self.n)
    let ratio = boltzmann_q(self.lambda)
    let mut cumulative = 0.0
    let mut weight = 1.0
    for idx = 0; idx < self.n; idx = idx + 1 {
      cumulative += norm * weight
      if cumulative >= q {
        return self.loc + idx.to_double()
      }
      weight *= ratio
    }
    self.loc + (self.n - 1).to_double()
  }
}

///|
pub impl RVDiscrete for Boltzmann with sf(self, x) {
  1.0 - self.cdf(x)
}

///|
pub impl RVDiscrete for Boltzmann with isf(self, q) {
  if q <= 0.0 || q > 1.0 {
    abort("q must be in (0, 1]")
  }
  self.ppf(1.0 - q)
}

///|
pub impl RVDiscrete for Boltzmann with logpmf(self, k) {
  match boltzmann_index(k, self.loc, self.n) {
    Some(index) => {
      let lambda = self.lambda
      let n = self.n.to_double()
      let log_norm = @math.ln(-@math.expm1(-lambda)) -
        @math.ln(-@math.expm1(-lambda * n))
      log_norm - lambda * index.to_double()
    }
    None => @double.neg_infinity
  }
}

///|
pub impl RVDiscrete for Boltzmann with logcdf(self, x) {
  let value = self.cdf(x)
  if value <= 0.0 {
    @double.neg_infinity
  } else {
    @math.ln(value)
  }
}

///|
pub impl RVDiscrete for Boltzmann with logsf(self, x) {
  let value = self.sf(x)
  if value <= 0.0 {
    @double.neg_infinity
  } else {
    @math.ln(value)
  }
}

///|
pub impl RVDiscrete for Boltzmann with entropy(self) {
  let norm = boltzmann_norm(self.lambda, self.n)
  let ratio = boltzmann_q(self.lambda)
  let mut weight = 1.0
  let mut entropy = 0.0
  for idx = 0; idx < self.n; idx = idx + 1 {
    let prob = norm * weight
    if prob > 0.0 {
      entropy -= prob * @math.ln(prob)
    }
    weight *= ratio
  }
  entropy
}

///|
pub impl RVDiscrete for Boltzmann with median(self) {
  self.ppf(0.5)
}

///|
pub impl RVDiscrete for Boltzmann with mean(self) {
  let n = self.n.to_double()
  let (mu, _, _, _, _, _) = boltzmann_moment_terms(self.lambda, n)
  self.loc + mu
}

///|
pub impl RVDiscrete for Boltzmann with varc(self) {
  let n = self.n.to_double()
  let (_, variance, _, _, _, _) = boltzmann_moment_terms(self.lambda, n)
  variance
}

///|
pub impl RVDiscrete for Boltzmann with std(self) {
  self.varc().sqrt()
}

///|
pub impl RVDiscrete for Boltzmann with interval(self, alpha : Double) {
  if alpha < 0.0 || alpha > 1.0 {
    abort("alpha must be in [0, 1]")
  }
  let lower = self.ppf(alpha / 2.0)
  let upper = self.ppf(1.0 - alpha / 2.0)
  (lower, upper)
}

///|
pub impl RVDiscrete for Boltzmann with skew(self) {
  let n = self.n.to_double()
  let (_, variance, trm, trm2, z, z_n) = boltzmann_moment_terms(self.lambda, n)
  if variance <= 0.0 || trm2 <= 0.0 {
    0.0
  } else {
    let n_cu = n * n * n
    let trm_cu = trm * trm * trm
    let numerator = z * (1.0 + z) * trm_cu - n_cu * z_n * (1.0 + z_n)
    let denom = trm2 * trm2.sqrt()
    if denom == 0.0 {
      0.0
    } else {
      numerator / denom
    }
  }
}

///|
pub impl RVDiscrete for Boltzmann with kurt(self) {
  let n = self.n.to_double()
  let (_, variance, trm, trm2, z, z_n) = boltzmann_moment_terms(self.lambda, n)
  if variance <= 0.0 || trm2 <= 0.0 {
    0.0
  } else {
    let n_sq = n * n
    let n_four = n_sq * n_sq
    let trm_sq = trm * trm
    let trm_four = trm_sq * trm_sq
    let numerator = z * (1.0 + 4.0 * z + z * z) * trm_four -
      n_four * z_n * (1.0 + 4.0 * z_n + z_n * z_n)
    let denom = trm2 * trm2
    if denom == 0.0 {
      0.0
    } else {
      numerator / denom
    }
  }
}
