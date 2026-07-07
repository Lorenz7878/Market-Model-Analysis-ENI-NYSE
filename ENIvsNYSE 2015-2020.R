ENI.P <- EniData2015.2020$Ultimo
NYSE.P <- NYSEnergy2015.2020$Ultimo
n <- length(ENI.P)

ENI.R <- ENI.P[-1]/ENI.P[-n]-1
NYSE.R <- NYSE.P[-1]/NYSE.P[-n]-1

plot(NYSE.R, ENI.R, pch=16)
abline(h=0,v=0)

# OLS
ols <- lm(ENI.R ~ NYSE.R)
hat.coef <- ols$coefficients
a <- -.5 ; int <- c(-a,a)
lines(int, hat.coef[1]+hat.coef[2]*int, col=2, lwd=2)

#kernel

ks <- ksmooth(NYSE.R, ENI.R, bandwidth = 0.05)  

lines(ks, col = "darkgreen", lwd = 2)
legend("topleft", legend = c("OLS", "Kernel"), col = c("red", "darkgreen"),
       lty = 1, lwd = 2, bty = "n")

confint(ols)
summary(ols)
anova(ols)

beta.hat <- summary(ols)$coefficients[2, 1]

residualerror <- resid(ols)
firm.specific.risk <- var(residualerror)
systematic.risk <- (beta.hat)^2 * var(NYSE.R)
var.eni.r <- systematic.risk + firm.specific.risk

residuals <- resid(ols)

#plot of density of residuals

plot(density(residuals), main="Density of residuals",col="red", lwd=2)

xseq <- seq(min(residuals), max(residuals), length = 100)
curve.norm <- dnorm(xseq, mean = mean(residuals), sd = sd(residuals))
lines(xseq, curve.norm, col = "blue", lwd = 2)


# Q-Q plot

qqnorm(residuals)
qqline(residuals, col="red")

#normality test bera and shapiro

shapiro.test(residuals)
library(tseries)
jarque.bera.test(residuals)

library(moments)

# Calcola skewness e kurtosis
skew <- skewness(residuals)
kurt <- kurtosis(residuals)

# Test di D'Agostino (skewness)
agostino.test(residuals)



# Neyman test
alpha <- 0.05
half.alpha <- alpha / 2
n1 <- n - 2
SE <- summary(ols)$coefficients[2, 2]
k.u <- qt(1 - half.alpha, n1) * SE
k.l <- qt(half.alpha, n1) * SE

if (beta.hat >= k.l && beta.hat <= k.u) {
  print("accept H0")
} else {
  print("reject H0")
}

effect.size <- beta.hat / SE



n <- n1 + 2  # number of observation

library(pwr)

power.result <- pwr.t.test(n = n, 
                           d = effect.size, 
                           sig.level = alpha, 
                           type = "one.sample", 
                           alternative = "two.sided")

print(power.result)

d.values <- seq(-1.5, 1.5, by = 0.01)
power.values <- sapply(d.values, function(d) {
  pwr.t.test(n = n, d = d, sig.level = alpha, 
             type = "one.sample", alternative = "two.sided")$power
})

plot(d.values, power.values, type = "l", lwd = 2, col = "blue",
     xlab = "beta",
     ylab = "Power",
     main = "power curve β̂")


abline(v = effect.size, col = "red", lty = 2) 

legend("bottomright",
       legend = c("Power curve"),
       col = c("blue"),
       lty = c(1, 2, 2), lwd = 2)



# Plot of Neyman Test
x_vals <- seq(-9 * SE, 9 * SE, length.out = 500)
t_vals <- (x_vals - 0) / SE
dens_vals <- dt(t_vals, df = n1)
dens_vals_rescaled <- dens_vals / max(dens_vals)

plot(x_vals, dens_vals_rescaled, type = "l", lwd = 2, col = "blue",
     xlab = expression(beta), ylab = "", yaxt = "n",
     main = "Neyman Test")

polygon(c(x_vals[x_vals <= k.l], k.l, k.l),
        c(dens_vals_rescaled[x_vals <= k.l], 0, 0),
        col = rgb(1, 0.6, 0.6, 0.6), border = NA)

polygon(c(k.u, x_vals[x_vals >= k.u], k.u),
        c(0, dens_vals_rescaled[x_vals >= k.u], 0),
        col = rgb(1, 0.6, 0.6, 0.6), border = NA)

abline(v = k.l, col = "darkred", lty = 2, lwd = 2)
abline(v = k.u, col = "darkred", lty = 2, lwd = 2)

abline(v = beta.hat, col = "red", lwd = 2)

legend("topleft",
       legend = c("Distribution t", "critical region", expression(hat(beta))),
       col = c("blue", rgb(1, 0.6, 0.6, 0.6), "red"),
       lty = c(1, NA, 1), lwd = c(2, NA, 2), pch = c(NA, 15, NA),
       pt.cex = 2, bg = "white")

#fisher test on beta.hat

test.s<-beta.hat/SE
p.value.beta <- 2* abs((1 - pt(test.s, n1)))

if (p.value.beta < alpha) {
  print("reject H0")
} else {
  print("accept H0")
}


#fisher test on alpha.hat

alpha.hat <- summary(ols)$coefficients[1, 1]
SE.alpha.hat<-summary(ols)$coefficients[1, 2]

t.statistic<- alpha.hat/SE.alpha.hat

p.value.alpha <- 2 * (1 - pt(t.statistic, n1))

if (p.value.alpha < alpha) {
  print("reject H0")
} else {
  print("accept H0")
}

#cyclical- non cyclical test  B0=1

B0<-1
K.alpha3<- qt(1-alpha,n1)*SE+B0

#set B0=1  H0:beta.hat=<B0 denotes a non cyclical-stock 
#          H1:beta.hat>=B0 denotes a cyclical stock

if (beta.hat < K.alpha3) {
  print("accept H0: non cyclical stock")
} else {
  print("accept H1: cyclical stock")
}

x_vals <- seq(B0 - 4 * SE, B0 + 4 * SE, length.out = 500)

# density function centered in B0
dens_vals <- dt((x_vals - B0) / SE, df = n1)

# Plot
plot(x_vals, dens_vals, type = "l", lwd = 2, col = "blue",
     xlab = expression(beta), ylab = "Density",
     main = "Cyclical vs Non-Cyclical Stock Test")

# critical region
polygon(c(K.alpha3, x_vals[x_vals >= K.alpha3], max(x_vals)),
        c(0, dens_vals[x_vals >= K.alpha3], 0),
        col = rgb(1, 0.6, 0.6, 0.6), border = NA)

# Line critical region
abline(v = K.alpha3, col = "red", lty = 2, lwd = 2)

# line of estimates of beta
abline(v = beta.hat, col = "darkgreen", lty = 2, lwd = 2)

# Legend
legend("topright",
       legend = c(expression(K[alpha]), expression(hat(beta))),
       col = c("red", "darkgreen"),
       lty = 2, lwd = 2, bty = "n")


