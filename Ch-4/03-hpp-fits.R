load("./data/whales.RData")
library(ascr)

cali.Capt <- list(bincapt = as.matrix(cali.CaptHist), toa = as.matrix(cali.ToA))

cali.hpp.hn <- fit.ascr(cali.Capt, cali.Traps, cali.Mask.Sea, detfn = "hn", sound.speed = 1500, hess = FALSE)
cali.hpp.hn.fixed <- fit.ascr(cali.Capt, cali.Traps, cali.Mask.Sea, detfn = "hn", sound.speed = 1500, hess = FALSE, fix = list(g0 = 1))
cali.hpp.hr <- fit.ascr(cali.Capt, cali.Traps, cali.Mask.Sea, detfn = "hr", sound.speed = 1500, sv = as.list(coef(cali.hpp.hn)), hess = FALSE) |>
  {\(x){
    fit.ascr(cali.Capt, cali.Traps, cali.Mask.Sea, detfn = "hr", sound.speed = 1500, sv = as.list(coef(x)), hess = FALSE)
  }}()
cali.hpp.th <- fit.ascr(cali.Capt, cali.Traps, cali.Mask.Sea, detfn = "th", sound.speed = 1500, hess = FALSE)
cali.hpp.hhn <- fit.ascr(cali.Capt, cali.Traps, cali.Mask.Sea, detfn = "hhn", sound.speed = 1500, hess = FALSE)

summary(cali.hpp.hn)
summary(cali.hpp.hn.fixed)
summary(cali.hpp.hr)
summary(cali.hpp.th)
summary(cali.hpp.hhn)

data.frame(detfn = c("hn", "hn.f", "hr", "th", "hhn"),
           Dc_km = c(coef(cali.hpp.hn)["D"] * 100 / 48, coef(cali.hpp.hn.fixed)["D"] * 100 / 48, coef(cali.hpp.hr)["D"] * 100 / 48, coef(cali.hpp.th)["D"] * 100 / 48, coef(cali.hpp.hhn)["D"] * 100 / 48),
           theta.1 = c(coef(cali.hpp.hn)["g0"], 1, coef(cali.hpp.hr)["g0"], coef(cali.hpp.th)["shape"], coef(cali.hpp.hhn)["lambda0"]),
           theta.2 = c(coef(cali.hpp.hn)["sigma"], coef(cali.hpp.hn.fixed)["sigma"], coef(cali.hpp.hr)["sigma"], coef(cali.hpp.th)["scale"], coef(cali.hpp.hhn)["sigma"]),
           theta.3 = c(NA, NA, coef(cali.hpp.hr)["z"], NA, NA),
           sigma.toa = c(coef(cali.hpp.hn)["sigma.toa"], coef(cali.hpp.hn.fixed)["sigma.toa"], coef(cali.hpp.hr)["sigma.toa"], coef(cali.hpp.th)["sigma.toa"], coef(cali.hpp.hhn)["sigma.toa"]),
           esa = c(coef(cali.hpp.hn, pars = "derived")["esa.1"], coef(cali.hpp.hn.fixed, pars = "derived")["esa.1"], coef(cali.hpp.hr, pars = "derived")["esa.1"], coef(cali.hpp.th, pars = "derived")["esa.1"], coef(cali.hpp.hhn, pars = "derived")["esa.1"]),
           AIC = c(AIC(cali.hpp.hn), AIC(cali.hpp.hn.fixed), AIC(cali.hpp.hr), AIC(cali.hpp.th), AIC(cali.hpp.hhn))) |>
  write.csv("code/whales/hpp-fits.csv")
