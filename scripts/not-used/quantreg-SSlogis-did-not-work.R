# build artificial data with multiplicative error
Dat <- NULL; Dat$x <- rep(1:25, 20)
set.seed(1)

Dat$y <- SSlogis(Dat$x, 10, 12, 2)*rnorm(500, 1, 0.1)

plot(prebinning_test$netDiv_CWM_std_tw,
     prebinning_test$NRI)

# fit first a nonlinear least-square regression

Dat.nls <- nls(NRI ~ poly(netDiv_CWM_std_tw, 4), data=prebinning_test); Dat.nls
lines(prebinning_test$netDiv_CWM_std_tw, 
      predict(Dat.nls, newdata=list(x=prebinning_test$netDiv_CWM_std_tw)), col=1)

# then fit the median using nlrq
Dat.nlrq <- nlrq(NRI ~ SSlogis(netDiv_CWM_std_tw, Asym, mid, scal), data=prebinning_test, tau=0.5, trace=TRUE)
lines(prebinning_test$netDiv_CWM_std_tw, 
      predict(Dat.nlrq, newdata=list(x=prebinning_test$netDiv_CWM_std_tw)), col=2)

# the 1st and 3rd quartiles regressions
Dat.nlrq <- nlrq(NRI ~ SSlogis(netDiv_CWM_std_tw, Asym, mid, scal), data=prebinning_test,  tau=0.25, trace=TRUE)
lines(prebinning_test$netDiv_CWM_std_tw, 
      predict(Dat.nlrq, newdata=list(x=prebinning_test$netDiv_CWM_std_tw)), col=3)

Dat.nlrq <- nlrq(NRI ~ SSlogis(netDiv_CWM_std_tw, Asym, mid, scal), data=prebinning_test,  tau=0.75, trace=TRUE)
lines(prebinning_test$netDiv_CWM_std_tw, 
      predict(Dat.nlrq, newdata=list(x=prebinning_test$netDiv_CWM_std_tw)), col=4)

# and finally "external envelopes" holding 95 percent of the data
# the 1st and 3rd quartiles regressions
Dat.nlrq <- nlrq(NRI ~ SSlogis(netDiv_CWM_std_tw, Asym, mid, scal), data=prebinning_test,  tau=0.025, trace=TRUE)
lines(prebinning_test$netDiv_CWM_std_tw, 
      predict(Dat.nlrq, newdata=list(x=prebinning_test$netDiv_CWM_std_tw)), col=5)

Dat.nlrq <- nlrq(NRI ~ SSlogis(netDiv_CWM_std_tw, Asym, mid, scal), data=prebinning_test,  tau=0.975, trace=TRUE)
lines(prebinning_test$netDiv_CWM_std_tw, 
      predict(Dat.nlrq, newdata=list(x=prebinning_test$netDiv_CWM_std_tw)), col=5)
leg <- c("least squares","median (0.5)","quartiles (0.25/0.75)",".95 band (0.025/0.975)")
legend(1, 12.5, legend=leg, lty=1, col=1:4)