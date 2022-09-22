# libraries
library("quantmod")
library("pastecs")
library("ggplot2")
library("ISOweek")
library("lubridate")
library("car")
library("olsrr")


# data

getSymbols("JPM",from="2007-01-01",to="2022-9-20",src="yahoo",auto.assign=TRUE)
head(JPM)


# Return

Pvec <- as.vector(JPM$JPM.Adjusted)
Rvec <- (Pvec[-1] / Pvec[-length(Pvec)]-1)*100
# add as a column to our data set
Rvec = c(0,Rvec)
data = data.frame(cbind(Rvec,JPM))
# create TS df
ts = data.frame(data$Rvec)
index = as.numeric(strftime(as.Date('index(ts)', "%d-%m-%Y"), "%u"))
#compute the log return
rvec<-diff(log(Pvec))
head(ts)
head(data)


# weak form efficiency

# summary statistics
data.frame(stat.desc(data, norm = TRUE))
# distribution
# Overlaid histograms
ggplot(data, aes(x=Rvec)) +
  geom_histogram(binwidth=.5, alpha=.5, position="identity")+ ggtitle("Return Distribution Plot")
# Density plots with semi-transparent fill
ggplot(data, aes(x=Rvec)) + geom_density(alpha=.3, fill="red")+ ggtitle("Return Distribution Plot")
ggplot(data, aes(x=Rvec
)) + geom_boxplot() + 
  guides(fill=FALSE) + coord_flip() + ggtitle("Return Box Plot")
summary(data$date)

# return vector: time series
# make date a col
data <- cbind(date = rownames(data), data)
rownames(data) <- 1:nrow(data)
data$date <- as.POSIXct(data$date, format="%Y-%m-%d" )
#data$date <- strptime(data$date, format = "%Y-%m-%d")
p <- ggplot(data, aes(x=date, y=Rvec)) + geom_line() +  theme_minimal() + 
  labs(x="Date", y="Return", title="Return Time Series Plot")
p


#Fitting the AR Model to the time series

AR <- arima(Rvec, order = c(1,0,0)) # first-order autoregressive model:
print(AR)
#plot the series along with the fitted values
ts.plot(Rvec)
AR_fit <- Rvec - residuals(AR)
points(AR_fit, type = "l", col = 2, lty = 2)
#plotting the series plus the forecast and 95% prediction intervals
ts.plot(Rvec)
AR_forecast <- predict(AR, n.ahead = 500)$pred
AR_forecast_se <- predict(AR, n.ahead = 500)$se
points(AR_forecast, type = "l", col = 2)
points(AR_forecast - 2*AR_forecast_se, type = "l", col = 2, lty = 2)
points(AR_forecast + 2*AR_forecast_se, type = "l", col = 2, lty = 2)
#test the goodness of fit
# Find AIC of AR
AIC(AR)
# Find BIC of AR
BIC(AR)



# assumptions check

#(1)check if the residuals are stationary.
## timeseries plot of the residuals indicates stationarity
AR1_resid = resid(AR)
plot(AR1_resid)
#(2)check if the residuals are weak WN
acf(AR1_resid, main="Sample ACF for the residuals")
Box.test(AR1_resid, lag = 5, type = "Ljung-Box", fitdf = 1)
#(3)check if the volatility is constant
## plot timeseries of residual^2 and fit a scatter plot smoother to highlight changes
##reuslt:  The volatility seems to be changing.
par(mfrow=c(1,2));par(mar=c(3,3,3,3))
plot(resid(AR)**2, type="l", col=1, main = expression(residual^2))
smoother = loess((resid(AR)**2) ~ seq(1,length(resid(AR))), span=0.1)
lines(seq(1,length(resid(AR))),fitted(smoother),col=2)
#We also check for autocorrelation in the squared residuals by using an ACF plot of residual^2
##result: The squared residuals seemed to be slightly autocorrelated
acf((resid(AR)**2), main=expression("sample ACF of "~ residual^2))
#Ljung Box test result: no autocorrelation in the squared residuals
Box.test(resid(AR)**2, lag = 5, type = "Ljung-Box", fitdf = 1)
#(4)check normality of residual
qqnorm(AR1_resid, datax = TRUE,
       xlab = "normal quantile",
       ylab = "sample quantile of residuals",
       main = "normal probability plot for residuals")
qqline(AR1_resid, datax=TRUE, col = 2)


# Day of week effect

# convert from calendar date to week date and back to calendar date
weekday = wday(ymd(as.Date(data$date))) - 1
data=cbind(weekday, Rvec)
head(data)
D1=rep(0,length(Rvec))
D2=rep(0,length(Rvec))
D3=rep(0,length(Rvec))
D4=rep(0,length(Rvec))
D5=rep(0,length(Rvec))
for(i in c(1:length(Rvec))){
  if (weekday[i]==1){D1[i]=1}
  if (weekday[i]==2){D2[i]=1}
  if (weekday[i]==3){D3[i]=1}
  if (weekday[i]==4){D4[i]=1}
  if (weekday[i]==5){D5[i]=1}
}
#linear regression model
multi.fit = lm(Rvec~D1+D2+D3+D4+D5)
summary(multi.fit)



# model assumptions

#the residual plot
plot(multi.fit$residuals, pch = 16, ylab= "Residuals")
abline(h = 0, lty = 3)
#test residual homoscedasticity
ncvTest(multi.fit)
plot(multi.fit)
#test residual normality
r=multi.fit$residuals
shapiro.test(r)
hist(r,col="bisque", freq=FALSE, main=NA)
qqPlot(r)
#test if Monday return is statistically significant
multi.fit_exceptMonday = lm(Rvec~D2+D3+D4+D5)
anova(multi.fit_exceptMonday, multi.fit)
# non-zero residuals mean
print(mean(resid(multi.fit)))
# test for residuals autocorrelation
durbinWatsonTest(multi.fit)
# autocorrelation of predictors
df <- data.frame(D1, D2, D3, D4, D5)
cor(df)
ols_vif_tol(multi.fit)
# correlation and covariance between residuals and Xs
dd <- data.frame(D1, D2, D3, D4, D5, resid(multi.fit))
cor(dd)
cov(dd)