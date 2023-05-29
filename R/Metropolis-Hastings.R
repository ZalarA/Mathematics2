#Metropolis, bivariate normal, symmetric jumps

#density function is equal to
#p(x,y)=(2*pi*(1-rho^2))^(-1)*(x^2-2*rho*x*y+y^2)
#The other possibility:
#p(x,y)=p(y|x)*p(x)=p(x|y)*p(y),
#where 
#Y|X is distributed normally with mean rho*y and std sqrt(1-rho^2)
#X|Y is distributed normally with mean rho*x and std sqrt(1-rho^2)
#We will use p(y|x)*p(x) below.

set.seed(1234); 
m = 40000; rho = .8; sgm = sqrt(1 - rho^2)
xc = yc = numeric(m) # vectors of state components
xc[1] = -3; yc[1] = 3 # arbitrary starting values
jl = 1; jr = 1 # l and r limits of proposed jumps
for (i in 2:m)
{
  xc[i] = xc[i-1]; yc[i] = yc[i-1] # if jump rejected
  xp = runif(1, xc[i-1]-jl, xc[i-1]+jr) # proposed x coord
  yp = runif(1, yc[i-1]-jl, yc[i-1]+jr) # proposed y coord
  #
  nmtr = dnorm(xp)*dnorm(yp, rho*xp, sgm)
  dntr = dnorm(xc[i-1])*dnorm(yc[i-1], rho*xc[i-1], sgm)
  r = nmtr/dntr # density ratio
  acc = (min(r, 1) > runif(1)) # jump if acc == T
  if (acc) {xc[i] = xp; yc[i] = yp}
}
x = xc[(m/2+1):m]; y = yc[(m/2+1):m] # states after burn-in
round(c(mean(x), mean(y), sd(x), sd(y), cor(x,y)), 4)
mean(diff(x)==0) # proportion or proposals rejected
par(mfrow = c(1,2), pty="s")
plot(xc[1:100], yc[1:100], xlim=c(-4,4), ylim=c(-4,4), type="l")
plot(x, y, xlim=c(-4,4), ylim=c(-4,4), pch=".")
par(mfrow = c(1,1), pty="m")

#Metropolis, bivariate normal, asymmetric jumps
set.seed(1234); m = 40000; rho = .8; sgm = sqrt(1 - rho^2)
xc = yc = numeric(m) # vectors of state components
xc[1] = -3; yc[1] = 3 # arbitrary starting values
jl = 1.25; jr = 0.75 # l and r limits of proposed jumps
for (i in 2:m)
{
  xc[i] = xc[i-1]; yc[i] = yc[i-1] # if jump rejected
  xp = runif(1, xc[i-1]-jl, xc[i-1]+jr) # proposed x coord
  yp = runif(1, yc[i-1]-jl, yc[i-1]+jr) # proposed y coord
  nmtr = dnorm(xp)*dnorm(yp, rho*xp, sgm)
  dntr = dnorm(xc[i-1])*dnorm(yc[i-1], rho*xc[i-1], sgm)
  r = nmtr/dntr # density ratio
  acc = (min(r, 1) > runif(1)) # jump if acc == T
  if (acc) {xc[i] = xp; yc[i] = yp}
}
x = xc[(m/2+1):m]; y = yc[(m/2+1):m] # states after burn-in
round(c(mean(x), mean(y), sd(x), sd(y), cor(x,y)), 4)
mean(diff(x)==0) # proportion or proposals rejected
par(mfrow = c(1,2), pty="s")
plot(xc[1:100], yc[1:100], xlim=c(-4,4), ylim=c(-4,4), type="l")
plot(x, y, xlim=c(-4,4), ylim=c(-4,4), pch=".")
par(mfrow = c(1,1), pty="m")

#Metropolis-Hastings correction, bivariate normal, asymmetric
set.seed(2008)
m = 100000; xc = yc = numeric(m); xc[1] = 3; yc[1] = -3
rho = .8; sgm = sqrt(1 - rho^2); jl = 1.25; jr = .75
for (i in 2:m)
{
  xc[i] = xc[i-1]; yc[i] = yc[i-1] # if no jump
  xp = runif(1, xc[i-1]-jl, xc[i-1]+jr)
  yp = runif(1, yc[i-1]-jl, yc[i-1]+jr)
  nmtr.r = dnorm(xp)*dnorm(yp, rho*xp, sgm)
  dntr.r = dnorm(xc[i-1])*dnorm(yc[i-1], rho*xc[i-1], sgm)
  nmtr.adj = dunif(xc[i-1], xp-jl, xp+jr)*
    dunif(yc[i-1], yp-jl, yp+jr)
  dntr.adj = dunif(xp, xc[i-1]-jl, xc[i-1]+jr)*
    dunif(yp, yc[i-1]-jl, yc[i-1]+jr)
  r = nmtr.r/dntr.r; adj = nmtr.adj/dntr.adj
  acc = (min(r*adj, 1) > runif(1))
  if (acc) {xc[i] = xp; yc[i] = yp}
}
x = xc[(m/2+1):m]; y = yc[(m/2+1):m]
round(c(mean(x), mean(y), sd(x), sd(y), cor(x,y)) ,4)
mean(diff(xc)==0); mean(pmax(x, y) > 1.25)

par(mfrow=c(1,2), pty="s")
jump = diff(unique(x)); hist(jump, prob=T, col="wheat")
plot(x, y, xlim=c(-4,4), ylim=c(-4,4), pch=".")
par(mfrow=c(1,1), pty="m")
thinned= seq(1, m/2, by=100)
