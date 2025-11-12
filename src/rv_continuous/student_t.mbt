///|
struct StudentT {
  df : Double
  loc : Double
  scale : Double
}

///|
pub fn StudentT::new(
  df : Double,
  loc? : Double = 0.0,
  scale? : Double = 1.0,
) -> StudentT {
  { df, loc, scale }
}

///|
pub impl RVContinuous for StudentT with rvs(
  self,
  size,
  rand? = @random.Rand::new(),
) {
  let normal = Normal::new()
  let chi2 = Chi2::new(self.df)
  Array::makei(size, _ => {
    let z = normal.rvs(1, rand~)[0]
    let v = chi2.rvs(1, rand~)[0]
    let t_standard = z / (v / self.df).sqrt()
    self.loc + self.scale * t_standard
  })
}

///|
pub impl RVContinuous for StudentT with pdf(self, x) {
  let z = (x - self.loc) / self.scale
  let k = (self.df + 1.0) / 2.0
  let factor = @kmath.gamma(k) /
    (self.scale * (self.df * @math.PI).sqrt() * @kmath.gamma(k - 0.5))
  factor * @math.pow(1 + z * z / self.df, -k)
}

///|
pub impl RVContinuous for StudentT with cdf(self, x) {
  let z = (x - self.loc) / self.scale
  let a = 0.5 * self.df
  let b = 0.5
  let x_beta = self.df / (self.df + z * z)
  if z >= 0.0 {
    1.0 - 0.5 * @sp.betainc(a, b, x_beta)
  } else {
    0.5 * @sp.betainc(a, b, x_beta)
  }
}

///|
pub impl RVContinuous for StudentT with ppf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  if q == 0.5 {
    return self.loc
  }
  let standard_t = StudentT::new(self.df)
  let initial_guess = if self.df > 30.0 {
    Normal::new().ppf(q)
  } else if self.df <= 2.0 {
    @math.tan(@math.PI * (q - 0.5))
  } else {
    let normal_approx = Normal::new().ppf(q)
    normal_approx * (self.df / (self.df - 2.0)).sqrt()
  }
  let ans_of_standard = @sp.func_inv(
    q,
    standard_t.cdf(_),
    standard_t.pdf(_),
    initial_guess,
    max_iterations=200,
  )
  let result_standard = match ans_of_standard {
    None => abort("No solution!")
    Some(val) => val
  }
  self.loc + self.scale * result_standard
}

///|
pub impl RVContinuous for StudentT with sf(self, x) {
  1.0 - self.cdf(x)
}

///|
pub impl RVContinuous for StudentT with isf(self, q) {
  if q <= 0.0 || q >= 1.0 {
    abort("q must be in (0, 1)")
  }
  self.ppf(1.0 - q)
}

///|
pub impl RVContinuous for StudentT with logpdf(self, x) {
  @math.ln(self.pdf(x))
}

///|
pub impl RVContinuous for StudentT with logcdf(self, x) {
  @math.ln(self.cdf(x))
}

///|
pub impl RVContinuous for StudentT with logsf(self, x) {
  @math.ln(self.sf(x))
}

///|
pub impl RVContinuous for StudentT with moment(self, n) {
  if n < 0 {
    abort("moment order must be non-negative")
  }
  if n == 0 {
    return 1.0
  }
  if n.to_double() >= self.df {
    abort(
      "Moment of order \{n} does not exist for t-distribution with df=\{self.df}",
    )
  }
  if self.scale == 1 && self.loc == 0 {
    return standard_t_moment(n, self.df)
  }
  let mut moment = 0.0
  for k in 0..=n {
    let binom_coeff = binomial_coefficient(n, k).to_double()
    let loc_power = @math.pow(self.loc, (n - k).to_double())
    let scale_power = @math.pow(self.scale, k.to_double())
    let std_moment = standard_t_moment(k, self.df)
    moment += binom_coeff * loc_power * scale_power * std_moment
  }
  moment
}

///|
pub impl RVContinuous for StudentT with mean(self) {
  if self.df <= 1.0 {
    abort("There is no mean when df<=1!")
  }
  self.loc
}

///|
pub impl RVContinuous for StudentT with median(self) {
  self.loc
}

///|
pub impl RVContinuous for StudentT with varc(self) {
  if self.df <= 2.0 {
    abort("There is no varc when df<=2!")
  }
  self.scale * self.scale * self.df / (self.df - 2.0)
}

///|
pub impl RVContinuous for StudentT with std(self) {
  if self.df <= 2.0 {
    abort("There is no std when df<=2!")
  }
  (self.scale * self.scale * self.df / (self.df - 2.0)).sqrt()
}

///|
pub impl RVContinuous for StudentT with entropy(self) {
  let nu = self.df
  let half_nu_plus_1 = (nu + 1.0) / 2.0
  let half_nu = nu / 2.0
  @math.ln(self.scale) +
  0.5 * @math.ln(nu * @math.PI) -
  @kmath.ln_gamma(half_nu_plus_1) +
  @kmath.ln_gamma(half_nu) +
  half_nu_plus_1 * (@kmath.digamma(half_nu_plus_1) - @kmath.digamma(half_nu))
}

///|
pub impl RVContinuous for StudentT with interval(self, alpha : Double) {
  if alpha <= 0.0 || alpha >= 1.0 {
    abort("alpha must be in (0, 1)")
  }
  (self.ppf(alpha / 2.0), self.ppf(1 - alpha / 2.0))
}

///|
pub impl RVContinuous for StudentT with skew(self) {
  if self.df <= 3.0 {
    abort("There is no skew when df<=3!")
  }
  0.0
}

///|
pub impl RVContinuous for StudentT with kurt(self) {
  if self.df <= 4.0 {
    abort("There is no kurt when df<=4!")
  }
  6.0 / (self.df - 4.0)
}

///|
fn standard_t_moment(n : Int, df : Double) -> Double {
  if n.to_double() >= df {
    abort(
      "Moment of order \{n} does not exist for t-distribution with df=\{df}",
    )
  }
  if n % 2 == 1 {
    0.0
  } else {
    let k = (n / 2).to_double()
    @math.pow(df, k) *
    @kmath.gamma(df / 2.0 - k) *
    @kmath.gamma(k + 0.5) /
    (@math.PI.sqrt() * @kmath.gamma(df / 2.0))
  }
}
