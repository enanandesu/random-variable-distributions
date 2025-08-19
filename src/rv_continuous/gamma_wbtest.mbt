///|
test {
  let gamma = Gamma::new(3.2, theta=2.3)
  inspect(gamma.pdf(5.7), content="0.1108087756194874")
  inspect(gamma.cdf(5.7), content="0.4022957773845705")
  inspect(gamma.ppf(0.57), content="7.3138369928122735")
  inspect(gamma.sf(5.7), content="0.5977042226154294")
  inspect(gamma.isf(0.57), content="5.951733204299607")
  inspect(gamma.logpdf(5.7), content="-2.1999493054633494")
  inspect(gamma.logcdf(5.7), content="-0.9105676962669869")
  inspect(gamma.logsf(5.7), content="-0.5146592584028679")
  inspect(gamma.entropy(), content="2.7208683892594947")
  inspect(gamma.moment(4), content="12125.667240959998")
  inspect(gamma.median(), content="6.609167718893405")
  inspect(gamma.mean(), content="7.359999999999999")
  inspect(gamma.varc(), content="16.927999999999997")
  inspect(gamma.std(), content="4.114365078599612")
  inspect(
    gamma.interval(0.05),
    content="(1.6253681703591927, 17.34242834066752)",
  )
  inspect(gamma.skew(), content="1.118033988749895")
  inspect(gamma.kurt(), content="1.875")
}

///|
test {
  let gamma = Gamma::new(1.5)
  let arr = gamma.rvs(10, None)
  inspect(
    arr,
    content="[0.9729569921961162, 4.261539076438271, 0.3094879277576365, 1.1122774301356915, 0.7569410589802311, 2.480793710856602, 0.0866516154470591, 1.2610716714444075, 0.3589018878899852, 1.7586655383458065]",
  )
}
