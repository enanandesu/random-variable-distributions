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
  let arr = c.rvs(10)
  inspect(
    arr,
    content="[10.264345884174421, 3.792653300580728, 5.955798100528671, 8.545474908392784, 5.917317049135864, 5.8773659867584245, 11.429732358768549, 1.2001530399360774, 6.982552028868126, 8.012671220936472]",
  )
}
