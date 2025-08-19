///|
test {
  let c = StudentT::new(4.3, loc=0.79, scale=1.53)
  inspect(c.pdf(5.7), content="0.00964821070268555")
  inspect(c.cdf(5.7), content="0.9852876771661471")
  inspect(c.ppf(0.57), content="1.0764340568299764")
  inspect(c.sf(5.7), content="0.014712322833852931")
  inspect(c.isf(0.57), content="0.5035659431700235")
  inspect(c.logpdf(5.7), content="-4.6409828002343705")
  inspect(c.logcdf(5.7), content="-0.014821622413996878")
  inspect(c.logsf(5.7), content="-4.219069848347401")
  inspect(c.entropy(), content="2.0880268458522355")
  inspect(c.moment(4), content="457.30691184086976")
  inspect(c.median(), content="0.79")
  inspect(c.mean(), content="0.79")
  inspect(c.varc(), content="4.376465217391305")
  inspect(c.std(), content="2.0920002909634845")
  inspect(c.interval(0.05), content="(-3.3435931943717563, 4.923593194371753)")
  inspect(c.skew(), content="0")
  inspect(c.kurt(), content="20.00000000000001")
}

///|
test {
  let c = StudentT::new(7.2)
  let arr = c.rvs(10, Some(b"11111ABCED33333abcde555556666678"))
  inspect(
    arr,
    content="[0.12467908609096608, 0.2768292605576927, 1.0504295090173157, 1.8827920640544755, -1.082633947187184, 0.04861538954791863, 1.0169856082229891, -0.6031720446586872, 1.6919164390125603, -1.7288510694355392]",
  )
}
