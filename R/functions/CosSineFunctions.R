# copied from: /home/yeung/projects/posttranscriptional_regulation/functions/CosSineFunctions.R
DelR.Dela <- function(a, b, n=4){
  # R = sqrt(a^2 + b^2)
  # derivative with respect to a
  # n is constant factor to scale amplitude. 
  # Mod of Fourier transform gives 1/2 mean to peak. Set n=2 for mean to peak. n=4 for peak to trough
  delRdela <- n * a / (sqrt(a ^ 2 + b ^ 2))
  return(delRdela)
}

DelR.Delb <- function(a, b, n=4){
  # R = sqrt(a^2 + b^2)
  # derivative with respect to b
  # Mod of Fourier transform gives 1/2 mean to peak. Set n=2 for mean to peak. n=4 for peak to trough
  delRdelb <- n * b / (sqrt(a ^ 2 + b ^ 2))
  return(delRdelb)
}

DelPhi.Dela <- function(a, b, omega = 2 * pi / 24){
  # phi = arctan(b / a)
  # derivative with respect to a
  # set omega to 1 if in radians, otherwise propagate omega
  delPhidela <- (-b / (a ^ 2 + b ^ 2)) / omega
  return(delPhidela)
}

DelPhi.Delb <- function(a, b, omega = 2 * pi / 24){
  # phi = arctan(b / a)
  # derivative with respect to b
  # set omega to 1 if in radians, otherwise propagate omega
  delPhidelb <- (a / (a ^ 2 + b ^ 2)) / omega
  return(delPhidelb)
}

Vectorize(GetAmp <- function(a, b, n=4){
  # Calculate amplitude and optionally propagate the errors
  # factor n gives ampltiude scaling
  R <- n * sqrt(a^2 + b^2)
  return(R)
}, vectorize.args = c("a", "b"))

Vectorize(GetAmp.se <- function(a, b, sig.a, sig.b, n=4){
  # print(paste("first term:", sig.a ^ 2 * DelR.Dela(a, b) ^ 2))
  # print(paste("second term:", sig.b ^ 2 + DelR.Delb(a, b) ^ 2))
  # factor n gives 
  sig.R <- sqrt(sig.a ^ 2 * DelR.Dela(a, b, n) ^ 2 + sig.b ^ 2 * DelR.Delb(a, b, n) ^ 2)
  return(sig.R)
}, vectorize.args = c("a", "b", "sig.a", "sig.b"))

Vectorize(AtanPositive <- function(b, a){
  # a: cos part
  # b: sin part
  # return atan2 but handle -pi to 0 to pi to 2pi
  phi <- atan2(b, a)
  phi <- sapply(phi, function(p){
    if (p < 0){
      p <- p + 2 * pi
    } 
    return(p)
  })
  return(phi)
}, vectorize.args = c("a", "b"))

Vectorize(GetPhi <- function(a, b, omega = 2 * pi / 24){
  # Calculate phi and optionally propagate the errors
  phi <- AtanPositive(b, a)
  phi <- phi / omega
  return(phi)
}, vectorize.args = c("a", "b"))

Vectorize(GetPhi.se <- function(a, b, sig.a, sig.b, omega = 2 * pi / 24){
  sig.phi <- sqrt(sig.a ^ 2 * DelPhi.Dela(a, b, omega) ^ 2 + sig.b ^ 2 * DelPhi.Delb(a, b, omega) ^ 2)
  return(sig.phi)
}, vectorize.args = c("a", "b", "sig.a", "sig.b"))

Vectorize(OscillateRelamp <- function(jmean, relamp, phase, jtime, w = 2 * pi / 24){
  y <- jmean * (1 + relamp * cos(w * jtime - phase))
  return(y)
})
