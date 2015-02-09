set.seed(55)

# ~~~
# ~~~ Hoff's code
# ~~~

# priors
mu0<-1.9  ; t20<-0.95^2
s20<-.01 ; nu0<-1

#data
y<-c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08)
n<-length(y) ; mean.y<-mean(y) ; var.y<-var(y)

S<-10000
PHI<-matrix(nrow=S,ncol=2)
PHI[1,]<-phi<-c( mean.y, 1/var.y)

### Gibbs sampling
for(s in 2:S) {

# generate a new theta value from its full conditional
mun<-  ( mu0/t20 + n*mean.y*phi[2] ) / ( 1/t20 + n*phi[2] )
t2n<- 1/( 1/t20 + n*phi[2] )
phi[1]<-rnorm(1, mun, sqrt(t2n) )

# generate a new sigma^2 value from its full conditional
nun<- nu0+n
s2n<- (nu0*s20 + (n-1)*var.y + n*(mean.y-phi[1])^2 ) /nun
phi[2]<- rgamma(1, nun/2, nun*s2n/2)

PHI[s,]<-phi         }


# ~~~
# ~~~ My code
# ~~~

PHI2 <- model.normal.semiconjugate.gibbs(
    phi.0 = c(mean.y, var.y)
  , gibbs.iters = 10000
  , mu.0 = mu0
  , tau2.0 = t20
  , nu.0 = nu0
  , sigma2.0 = s20
  , y.bar = mean.y
  , n.obs = n
  , s2 = var.y)


# ~~~
# ~~~ Compare output with Hoff's code 
# ~~~

(his <- quantile(PHI[,1],c(.025,.5,.975)))
(mine <- quantile(PHI2[,1],c(.025,.5,.975)))
for (i in 1:3) {
    expect_equal(his[i], mine[i], tolerance = 0.01)
}

(his <- quantile(PHI[,2],c(.025,.5, .975)))
(mine <- quantile(1/PHI2[,2],c(.025,.5, .975)))
for (i in 1:3) {
    expect_equal(his[i], mine[i], tolerance = 0.05)
}
