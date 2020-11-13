

## ----echo=FALSE, include= TRUE---------------------------------------------------------------------------
lambda <- 2
a <- 2
b <- 2.2
plot({function (x) dbeta(x, a, b)}, main = "Loi beta de paramètres alpha=2, beta=2.2", ylab="f(y)")
plot({function (x) (dbeta(x, a, b)^lambda - 1)/lambda}, -0.2, 1.2, main = "La même loi beta après une transformation de Box et Cox avec lambda=2", ylab="h(f(y))")


## ----echo=FALSE------------------------------------------------------------------------------------------
# Initialisation vars
lambda_ <- 0.3
a <- 5
b <- 1
variance <- 2
n=50


## ----echo=FALSE, include= TRUE---------------------------------------------------------------------------
#Q1
set.seed(999)
X_obs <- rnorm(n) #loi gaussien centrée réduite
Epsilon <- rnorm(n=50, mean=0, sd = sqrt(variance))
# Par la définition
Z <- t(a + b%*%X_obs + Epsilon) #vector
# On considère que toutes les observations du jeu de données Y_i sont positives
Y <- (lambda_*Z +1)^(1/lambda_)
X <- as.matrix(cbind(1, X_obs))



## ----echo=FALSE, include= TRUE---------------------------------------------------------------------------
resY <- lm(Y~X)
resZ <- lm(Z~X)

{plot(resZ$fitted, Z, xlab = "vars ajustées", ylab = "vars observées",main="Pour Z")
abline(0, 1)}
{plot(resY$fitted, Y, xlab = "vars ajustées", ylab = "vars observées",main="Pour Y")
abline(0, 1)}


## ----echo=FALSE, include=TRUE----------------------------------------------------------------------------
#Etudier les résidus 
library(MASS)
#Z~X
{plot(resZ$fitted,studres(resZ),col=3,pch=3) #studentisés par validation croisée
points(resZ$fitted,stdres(resZ), col=2,pch=2 )# fitted to variance 1
points(resZ$fitted,resZ$residuals,main="différents résidus Z~X")
abline(h=2,lty=2)
abline(h=-2,lty=2)
legend("bottomright",c("resid", "stdres", "studres"),col=1:3, pch=1:3,cex=0.5)
}


## ----echo=FALSE, include=TRUE----------------------------------------------------------------------------
#Y~X
{plot(resY$fitted,studres(resY),col=3,pch=3) #studentisés par validation croisée
points(resY$fitted,stdres(resY), col=2,pch=2 )# fitted to variance 1
points(resY$fitted,resY$residuals,main="différents résidus Y~X")
abline(h=2,lty=2)
abline(h=-2,lty=2)
legend("bottomright",c("resid", "stdres", "studres"),col=1:3, pch=1:3,cex=0.5)
}


## ----echo=FALSE, include=TRUE----------------------------------------------------------------------------
qqnorm(studres(resZ),main="graphe quantile-quantile Z~X")
qqline(stdres(resZ)) # passe par le 1er et le 3ème quartile

qqnorm(studres(resY),main="graphe quantile-quantile Y~X")
qqline(stdres(resY)) # passe par le 1er et le 3ème quartile


## ----echo=FALSE, include=TRUE----------------------------------------------------------------------------
shapiro.test(studres(resZ))
shapiro.test(studres(resY))


## ----echo=FALSE------------------------------------------------------------------------------------------
#Q2
# Créer la matrice X du plan d'expérience
X <- as.matrix(cbind(1, X_obs))
# Interprétation les codes
Q = diag(1,n) - X%*%solve(t(X)%*%X)%*%t(X) 
Lmle = function(Z){
  n = length(Z)
  sig2 = (t(Z)%*%Q%*%Z) / n
  -n/2 * log(sig2)
}


## ----echo=FALSE, include=TRUE----------------------------------------------------------------------------
# Coder la fonction lmin(lambda, Y) qui calcule −Lmax.
vec1 = rep(1, 50)
Lmin = function(lambda, Y=Y) {
  Z_ <- ((Y^lambda)-1)/lambda
  -Lmle(Z_) + (1-lambda) * vec1%*%log(abs(Y)) + n/2 + (n/2)*log(2*pi)
}
# Tracer −Lmax
lambdas <- seq(0.01,2,0.01)
Vlmin = Vectorize(Lmin,"lambda")
plot(lambdas,Vlmin(lambdas, Y),
     main="-Lmax en fonction de lambda",
     ylab="",
     type="l",
     col="blue")



## ----echo=FALSE, warning=FALSE---------------------------------------------------------------------------
resopt = nlm(Lmin,Y=Y,p=2,hessian=TRUE)
lambda_est <- resopt$estimate
var_lambda_est <- 1/resopt$hessian


## ----message=FALSE, include=FALSE------------------------------------------------------------------------
alpha <- 0.05
q <- qnorm(1-alpha/2)
lambda_est - q*sqrt(var_lambda_est) #lower bound
lambda_est + q*sqrt(var_lambda_est) #upper bound


## ----warning=FALSE, include=FALSE------------------------------------------------------------------------
q_chi <- qchisq(1-alpha, df=1) #quantile
q_chi
lambda_0 <- c(1, 1/2, 0.3, 0.000001)
lambda_0
W_obs <-  (lambda_est - lambda_0)^2 / var_lambda_est
W_obs
p_value <- pchisq(W_obs, df=1, lower.tail = FALSE)
p_value
p_value-alpha>0 # conserve H_0?


## ----echo=FALSE, message=TRUE, include=TRUE--------------------------------------------------------------
lambda_0[4] <- 0.000001
for(i in 1:4){
  print(paste("lambda=", lambda_0[i]))
  TRV <- 2 *(Lmin(lambda_0[i], Y) - Lmin(lambda_est, Y))
  print(paste("TRV: ", TRV))
  p_value <- pchisq(TRV, df=1, lower.tail = FALSE)
  print(paste("p_value: ", p_value))
  print(paste("conserve H_0?", p_value-alpha>0)) # conserve H_0?
}



## ----include=FALSE---------------------------------------------------------------------------------------
################################# Merci de décommenter ci vous souhaitez l'exécuter ########################
#library('car')
#res <- powerTransform(Y~X, family="bcPower")
#summary(res)


## ----include=FALSE---------------------------------------------------------------------------------------
# Lire les données, vérfier que le data.frame obtenu comporte 27 observations.
df = read.csv2("NbCycleRupture.csv")
head(df)
dim(df)
str(df) #structure de data frame


## ----echo=FALSE, include= TRUE---------------------------------------------------------------------------
#M1
res1 <- lm(y~.,data=df) 
{plot(res1$fitted, df$y, xlab = "vars ajustées", ylab = "vars observées",main="M1:y~x1+x2+x3")
abline(0, 1)}

longueur <- 50*df$x1+300
amplitude <- df$x2+9
chargement <- df$x3*5 + 45
res2=lm(df$y~ longueur + amplitude +chargement)
{plot(res2$fitted, df$y, xlab = "vars ajustées", ylab = "vars observées",main="y~longueur+amplitude+chargement")
abline(0, 1)}


## ----echo=FALSE, include= TRUE---------------------------------------------------------------------------
summary(res1)
summary(res2)


## ----echo=FALSE, include=TRUE----------------------------------------------------------------------------
#M2
res3 <- lm(y~ df$x1 + df$x2 + df$x3 + I(df$x1^2) + I(df$x2^2) + I(df$x3^2) + I(df$x1*df$x2) + I(df$x1*df$x3) + I(df$x2*df$x3), data=df) 
{plot(res3$fitted, df$y, xlab = "vars ajustées", ylab = "vars observées",main="M2")
abline(0, 1)}


## ----echo=FALSE, include=TRUE----------------------------------------------------------------------------
summary(res3)


## ----echo=FALSE, include=FALSE---------------------------------------------------------------------------
anova(res1, res3)



