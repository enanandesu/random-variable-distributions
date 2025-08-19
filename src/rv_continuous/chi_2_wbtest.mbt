///|
test {
  let c = Chi2::new(7.0)
  inspect(c.pdf(5.7), content="0.11933467990350571")
  inspect(c.cdf(5.7), content="0.42481764624201707")
  inspect(c.ppf(0.57), content="6.989250037623171")
  inspect(c.sf(5.7), content="0.575182353757983")
  inspect(c.isf(0.57), content="5.743486110561044")
  inspect(c.logpdf(5.7), content="-2.125823297205622")
  inspect(c.logcdf(5.7), content="-0.8560952697994322")
  inspect(c.logsf(5.7), content="-0.5530681514916809")
  inspect(c.entropy(), content="2.6362291812939116")
  inspect(c.moment(4), content="9009")
  inspect(c.median(), content="6.3458111955215255")
  inspect(c.mean(), content="7")
  inspect(c.varc(), content="14")
  inspect(c.std(), content="3.7416573867739413")
  inspect(c.interval(0.05), content="(1.6898691806773758, 16.01276427462951)")
  inspect(c.skew(), content="1.0690449676496976")
  inspect(c.kurt(), content="1.7142857142857142")
}

///|
test {
  let c = Chi2::new(7)
  let arr = c.rvs(10, None)
  inspect(
    arr,
    content="[6.442561285106298, 5.288446338554324, 6.611708627050616, 6.001664104247955, 9.923541034831635, 2.5574672818973974, 3.9695902906961487, 5.801548824773635, 4.268895608967385, 8.621065951317108]",
  )
}
