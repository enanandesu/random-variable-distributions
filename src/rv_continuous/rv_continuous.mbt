///|
pub trait RVContinuous {
  rvs(Self, size : Int, rand? : @random.Rand) -> Array[Double]
  pdf(Self, x : Double) -> Double
  cdf(Self, x : Double) -> Double
  ppf(Self, q : Double) -> Double
  sf(Self, x : Double) -> Double
  isf(Self, q : Double) -> Double
  logpdf(Self, x : Double) -> Double
  logcdf(Self, x : Double) -> Double
  logsf(Self, x : Double) -> Double
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
