///|
struct Normal {
  loc : Double
  scale : Double
}

///|
pub fn Normal::new(loc~ : Double = 0.0, scale~ : Double = 1.0) -> Normal {
  { loc, scale }
}

///|
impl RVContinuous for Normal with rvs(self, size, seed) {
  let samples = Array::new(capacity=size)
  let mut i = 0
  let rdm = match seed {
    None => @random.Rand::new()
    Some(val) => @random.Rand::chacha8(seed=val)
  }
  while i < size {
    let u1 = rdm.double()
    let u2 = rdm.double()
    let u1_safe = if u1 <= 0.0 { 1.0e-10 } else { u1 }
    let factor = (-2.0 * @math.ln(u1_safe)).sqrt()
    let theta = 2.0 * @math.PI * u2
    let z1 = factor * @math.cos(theta)
    let z2 = factor * @math.sin(theta)
    samples.push(self.loc + self.scale * z1)
    i += 1
    if i < size {
      samples.push(self.loc + self.scale * z2)
      i += 1
    }
  }
  samples
}

///|
impl RVContinuous for Normal with pdf(self, x) {
  let z = (x - self.loc) / self.scale
  let factor = 1.0 / (self.scale * (2.0 * @math.PI).sqrt())
  factor * @math.exp(-0.5 * z * z)
}

///|
impl RVContinuous for Normal with cdf(self, x) {
  let z = (x - self.loc) / self.scale
  0.5 * (1.0 + @kmath.erf(z / 2.0.sqrt()))
}

///|
impl RVContinuous for Normal with ppf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  self.loc + self.scale * 2.0.sqrt() * @kmath.erfinv(2.0 * q - 1.0)
}

///|
impl RVContinuous for Normal with sf(self, x) {
  1.0 - self.cdf(x)
}

///|
impl RVContinuous for Normal with isf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  self.ppf(1.0 - q)
}

///|
impl RVContinuous for Normal with logpdf(self, x) {
  @math.ln(self.pdf(x))
}

///|
impl RVContinuous for Normal with logcdf(self, x) {
  @math.ln(self.cdf(x))
}

///|
impl RVContinuous for Normal with logsf(self, x) {
  @math.ln(self.sf(x))
}

///|
impl RVContinuous for Normal with moment(self, n) {
  if n < 0 {
    abort("moment order must be non-negative")
  }
  if n == 0 {
    return 1.0
  }
  let mut moment = 0.0
  for k in 0..=n {
    let binom_coeff = binomial_coefficient(n, k).to_double()
    let mu_power = @math.pow(self.loc, (n - k).to_double())
    let sigma_power = @math.pow(self.scale, k.to_double())
    let std_moment = standard_normal_moment(k)
    moment += binom_coeff * mu_power * sigma_power * std_moment
  }
  moment
}

///|
impl RVContinuous for Normal with mean(self) {
  self.loc
}

///|
impl RVContinuous for Normal with median(self) {
  self.loc
}

///|
impl RVContinuous for Normal with varc(self) {
  self.scale * self.scale
}

///|
impl RVContinuous for Normal with std(self) {
  self.scale
}

///|
impl RVContinuous for Normal with entropy(self) {
  0.5 * (@math.ln(2.0 * @math.PI * self.scale * self.scale) + 1.0)
}

///|
impl RVContinuous for Normal with interval(self, alpha : Double) {
  if alpha <= 0.0 || alpha >= 1.0 {
    abort("alpha must be in (0, 1)")
  }
  let z_alpha_half = Normal::new().ppf(1.0 - alpha / 2.0)
  let margin = z_alpha_half * self.scale
  (self.loc - margin, self.loc + margin)
}

///|
impl RVContinuous for Normal with skew(_) {
  0.0
}

///|
impl RVContinuous for Normal with kurt(_) {
  0.0
}

///|
fn binomial_coefficient(n : Int, k : Int) -> Int {
  if k > n || k < 0 {
    0
  } else if k == 0 || k == n {
    1
  } else {
    let mut result = 1
    let k_min = if k < n - k { k } else { n - k }
    for i in 0..<k_min {
      result = result * (n - i) / (i + 1)
    }
    result
  }
}

///|
fn standard_normal_moment(n : Int) -> Double {
  match n {
    0 => 1.0
    1 => 0.0
    _ => if n % 2 == 1 { 0.0 } else { double_factorial(n - 1).to_double() }
  }
}

///|
fn double_factorial(n : Int) -> Int {
  if n <= 1 {
    1
  } else {
    n * double_factorial(n - 2)
  }
}

///|
fn derive_seed(seed : Bytes) -> Bytes {
  // 每次调用都会吐出一个新的 32 字节种子
  let r = @random.Rand::chacha8(seed~)
  Bytes::makei(32, _ => r.uint(limit=256U).to_byte())
}
