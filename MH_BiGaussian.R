library(VineCopula)
genY <- function (rho, N){ #Draw correlated normals
  X1 = rnorm(N)
  X2 = rnorm(N)
  X3 = rho*X1 + sqrt(1-rho^2)*X2
  Y1 = X1
  Y2 = X3
  Y = matrix(c(Y1,Y2),nrow = N, ncol=2)
  return (Y)
}


Co_PDF<- function(u1, u2, par){
  u0=BiCopPDF(u1, u2, family = 1, par = par)
  return(u0)
}

l_ratio <- function(Y,rho,rho_) #Likelihood ratio
  return (
    sum(log(Co_PDF(Y[,1],Y[,2],rho_)) - log(Co_PDF(Y[,1],Y[,2],rho)))
  )


  prior_ratio =function(rho,rho_)
    return(((1-rho_^2)^500)/((1-rho^2)^500))
 
posterior_ratio<- function(Y,rho,rho_){ #Use Bayes Formula
  return(l_ratio(Y,rho,rho_)* prior_ratio(rho,rho_))
}

rho_t=rho = 0.65
Y = genY(rho,1000)
Y=pobs(Y)
plot(Y, col='brown', cex=.5, pch=14)
#burn_in = 1000
#iterations = burn_in + 1000
iterations = 200
rho_0 = 0.5
rho = rho_0
s = c(0)
for (i in 1:iterations){
  rho_ = runif(1, min = rho -0.2, max = rho+0.2)
  alpha = min(1, 1/posterior_ratio(Y,rho,rho_))
  if (runif(1)<alpha){
    rho = rho_
    
  }
  # if (i >burn_in)
  #   s = c(s,rho)
  #print(s)
  s = c(s,rho)
  #print(alpha)
}
n = seq_along(s)
m = cumsum(s)/n
m2 = cumsum(s*s)/n
v = (m2 -m*m)*(n/(n-1))
#plot(m, col='red', pch=19, cex=.2 )



plot(m,  type = "b", pch = 19,
     col = "red", xlab = "MCMC iteration", ylab = "Rho", 
     lty = 3, lwd = 2, cex=.5, ylim=c(mean(m)-.3, mean(m)+.3))
abline(h=rho_t,col = "blue", 
       lty = 2, lwd = 2)
# # 3. Add a second line
# lines(rho_t, pch = 18, col = "blue", type = 'b', 
#       lty = 2, lwd = 2)


# 4. Add a legend to the plot and set legend lty
legend("bottomright", legend = c("Estimated Rho", "True Rho" ),
       col = c("red", "blue"), lty = 1:2, cex = 0.9 )

legend("topright", legend = "Increase MCMC iterations for a better result", col='lightgray', cex = 0.5 )
