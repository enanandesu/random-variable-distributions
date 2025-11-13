///|
pub trait RVContinuous {
  rvs(Self, size : Int, rand? : @random.Rand) -> Array[Double]
  pdf(Self, x : Double) -> Double
  cdf(Self, x : Double) -> Double
  ppf(Self, q : Double) -> Double
  sf(Self, x : Double) -> Double
  isf(Self, q : Double) -> Double
  logpdf(Self, x : Double) -> Double = _
  logcdf(Self, x : Double) -> Double = _
  logsf(Self, x : Double) -> Double = _
  moment(Self, n : Int) -> Double
  entropy(Self) -> Double
  median(Self) -> Double
  mean(Self) -> Double
  varc(Self) -> Double
  std(Self) -> Double
  interval(Self, alpha : Double) -> (Double, Double)
  skew(Self) -> Double
  kurt(Self) -> Double
}

///|
impl RVContinuous with logpdf(self, x) {
  @math.ln(self.pdf(x))
}

///|
impl RVContinuous with logcdf(self, x) {
  @math.ln(self.cdf(x))
}

///|
impl RVContinuous with logsf(self, x) {
  @math.ln(self.sf(x))
}
