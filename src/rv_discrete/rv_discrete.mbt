///|
pub trait RVDiscrete {
  rvs(Self, size : Int, rand? : @random.Rand) -> Array[Double]
  pmf(Self, k : Double) -> Double
  cdf(Self, x : Double) -> Double
  ppf(Self, q : Double) -> Double
  sf(Self, x : Double) -> Double
  isf(Self, q : Double) -> Double
  logpmf(Self, k : Double) -> Double = _
  logcdf(Self, x : Double) -> Double = _
  logsf(Self, x : Double) -> Double = _
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
impl RVDiscrete with logpmf(self, x) {
  @math.ln(self.pmf(x))
}

///|
impl RVDiscrete with logcdf(self, x) {
  @math.ln(self.cdf(x))
}

///|
impl RVDiscrete with logsf(self, x) {
  @math.ln(self.sf(x))
}
