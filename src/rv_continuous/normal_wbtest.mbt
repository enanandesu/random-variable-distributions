///|
test {
  let norm = Normal::new(loc=3.2, scale=2.3)
  inspect(norm.pdf(5.7), content="0.09607896282151948")
  inspect(norm.cdf(5.7), content="0.8614719874836898")
  inspect(norm.ppf(0.57), content="3.605660578995981")
  inspect(norm.sf(5.7), content="0.1385280125163102")
  inspect(norm.isf(0.57), content="2.7943394210040196")
  inspect(norm.logpdf(5.7), content="-2.3425848962153912")
  inspect(norm.logcdf(5.7), content="-0.14911273952227494")
  inspect(norm.logsf(5.7), content="-1.9766827173703665")
  inspect(norm.entropy(), content="2.251847656139777")
  inspect(norm.moment(4), content="513.8275000000001")
  inspect(norm.median(), content="3.2")
  inspect(norm.mean(), content="3.2")
  inspect(norm.varc(), content="5.289999999999999")
  inspect(norm.std(), content="2.3")
  inspect(
    norm.interval(0.05),
    content="(-1.307917164442122, 7.707917164442122)",
  )
  inspect(norm.skew(), content="0")
  inspect(norm.kurt(), content="0")
}

///|
test {
  let norm = Normal::new()
  let arr = norm.rvs(10, Some(b"12345678912345678998765432101234"))
  inspect(
    arr,
    content="[-0.38205584043417495, -0.23346926482400399, 0.13253431184974715, -1.6397320023872184, -0.37927401482268985, 0.497241294990723, -0.22747079179606697, 0.24827860372961044, -1.3123383074025374, -0.5220240942281283]",
  )
}
