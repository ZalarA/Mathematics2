#A Markov chain with 2 states and a transition matrix K1:

K1 = matrix(c( 0.5, 0.5,
              0.5, 0.5), nrow=2, ncol=2,byrow = T)
K1

#An arbitrary starting state 
v= runif(2)
#A starting state must be a distribution vector
v=v/sum(v)
#We do 30 steps of the Markov chain
for (i in 1:30)
  {v=v%*%K1
  print(v)}


#A Markov chain with 2 states and a transition matrix K2
#with huge tendency to remain in the present state:

K2 = matrix(c( 0.95, 0.05,
               0.05, 0.95), nrow=2, ncol=2,byrow = T)
K2

#An arbitrary starting state 
v= runif(2)
#A starting state must be a distribution vector
v=v/sum(v)
#We do 30 steps of the Markov chain
for (i in 1:30)
{v=v%*%K2
print(v)}

#A Markov chain with 2 states and a transition matrix K2
#with huge tendency to change the state:

K3 = matrix(c( 0.05, 0.95,
               0.95, 0.05), nrow=2, ncol=2,byrow = T)
K3

#An arbitrary starting state 
v= runif(2)
#A starting state must be a distribution vector
v=v/sum(v)
#We do 30 steps of the Markov chain
for (i in 1:30)
{v=v%*%K3
print(v)}


#Computation of E_Pi[f] for f(w1)=5,f(w2)=-1 using the 
#Markov chains above.

#first Markov chain 1
set.seed(1237)
m = 20000; n = 1:m; x = numeric(m); x[1] = 1
f=numeric(m)

# Simulation
for (i in 2:m)
{
  x[i]=rbinom(1,1,0.5)
  if (x[i-1]==0) f[i-1]=5
  else f[i-1]=-1
}
f
y = cumsum(x)/n #

# Results

plot(x[1:20], type="b", xlab="Step", ylab="State")
average=cumsum(f)/n
average[(m-10):m]
plot(cumsum(f)/n)
sd(average)
autocor=acf(f)
sum(autocor$acf[2:30])*2+var(f)

#Markov chain 2
set.seed(1237)
m = 20000; n = 1:m; x = numeric(m); x[1] = 1
f=numeric(m)

# Simulation
for (i in 2:m)
{
  r=rbinom(1,1,0.95)
  if (r==1) x[i]=x[i-1]
  else x[i]=1-x[i-1]
  if (x[i-1]==0) f[i-1]=5
  else f[i-1]=-1
}
f
y = cumsum(x)/n #

# Results

plot(x[1:20], type="b", xlab="Step", ylab="State")
average=cumsum(f)/n
average[(m-10):m]
plot(cumsum(f)/n)
sd(average)
autocor=acf(f)
sum(autocor$acf[2:30])*2+var(f)

#Markov chain 3
set.seed(1237)
m = 20000; n = 1:m; x = numeric(m); x[1] = 1
f=numeric(m)

# Simulation
for (i in 2:m)
{
  r=rbinom(1,1,0.05)
  if (r==1) x[i]=x[i-1]
  else x[i]=1-x[i-1]
  if (x[i-1]==0) f[i-1]=5
  else f[i-1]=-1
}
f
y = cumsum(x)/n #

# Results

plot(x[1:20], type="b", xlab="Step", ylab="State")
average=cumsum(f)/n
average[(m-10):m]
plot(cumsum(f)/n)
sd(average)
autocor=acf(f)
sum(autocor$acf[2:30])*2+var(f)