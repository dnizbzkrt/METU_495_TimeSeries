library(ggplot2)
library(fUnitRoots)
library(gridExtra)
library(fpp)
library(anomalize)
library(tibbletime)
library(tidyverse)
library(forecast)
library(tseries)
library(pdR)
library(uroot)
library(TSA)
#Taken from https://www.investing.com/rates-bonds/u.s.-10-year-bond-yield
df <- read.csv("10year.csv")

autoplot(ts(df$Close,start = 1998,frequency = 12))+
  labs(y = "Yield",title = "US 10-year Bond")+theme_minimal()
df$Date <- as.Date(df$Date)
df <- df[,c("Date","Close")]
test <- tail(df,36)
test <-ts(test$Close, start = c(2023, 1), frequency = 12)
train <- head(df, -36)
train <- as_tibble(train)
# Decompose time series data and detect anomalies
train %>% 
  time_decompose(Close, method = "stl", frequency = "auto", trend = "auto") %>%
  anomalize(remainder, method = "gesd", alpha = 0.05, max_anoms = 0.2) %>%
  plot_anomaly_decomposition()
train %>% 
  time_decompose(Close) %>%
  anomalize(remainder,alpha=0.05) %>%
  time_recompose() %>%
  plot_anomalies(time_recomposed = TRUE, ncol = 3, alpha_dots = 0.5)
train %>% 
  time_decompose(Close) %>%
  anomalize(remainder) %>%
  time_recompose() %>%
  filter(anomaly == 'Yes')
# Clean anomalies
train <- train %>% 
  time_decompose(Close) %>%
  anomalize(remainder,alpha=0.05) %>%
  clean_anomalies()

ts_data <- ts(train$observed_cleaned, start = 1998, frequency = 12)
train <- tsclean(ts_data)

p1<-ggAcf(train)+labs(title = "Autocorrelation Function Plot")+theme_minimal()
p2<-ggPacf(train)+labs(title = "Autocorrelation Function Plot")+theme_minimal()
grid.arrange(p1,p2,nrow=1)
autoplot(train)+labs(y="Yield",title = "train")+theme_minimal()
lambda <- BoxCox.lambda(train)
lambda
#The optimal lambda value is .79.
ts_data_t <-BoxCox(train,lambda) 
autoplot(ts_data_t)+labs(y="Log Yield",title = "train")+theme_minimal()
#since  therewas no change invariation we conti without transformation
kpss.test(train,c("Level")) #Not stationary
#Since p value is less than α we reject H0. Therefore, we don’t have enough evidence to claim that we have stationary data.
kpss.test(train,c("Trend")) #According to KPSS test the series has stoc trend.
mean(train)
#Since mean is not 0. we will use these adf tests
adfTest(train, lags=1, type="c")#lag choosen 1 because 1.lag is significant in  pacf
#Since p value is greater than α=0.05 , we fail to reject H0 The series is not stationary
#Now, we add trend term to the model
adfTest(train, lags=1, type="ct")
#Since p value is lower than α=0.05 , we reject H0. It means that
#we have non stationary system having determinictic trend. Therefore, we apply differencing using diff function.
out<-HEGY.test(wts=train, itsd=c(1,1,0), regvar=0, selectlags=list(mode="signf", Pmax=NULL))
out$stats
#tpi_1 is greater than 0.05. so, there is regular unit root
#fpi_11:12 is lower than 0.05.so there is no seasonal unit root
ndiffs(train)
#Thus, I need to take one regular difference.
d_ts_data_t <- diff(train)
autoplot(d_ts_data_t,main="Time Series Plot of the Differenced Data",ylab="ddata",xlab="Time")+theme_minimal()
ggPacf(d_ts_data_t,lag.max = 36)+labs(title = "Autocorrelation Function Plot")+theme_minimal()
ggAcf(d_ts_data_t,lag.max = 36)+labs(title = "Autocorrelation Function Plot")+theme_minimal()
grid.arrange(ggPacf(d_ts_data_t,lag.max = 36)+labs(title = "Autocorrelation Function Plot")+theme_minimal(),
             ggAcf(d_ts_data_t,lag.max = 36)+labs(title = "Autocorrelation Function Plot")+theme_minimal()
             ,nrow=1)
mean(d_ts_data_t)
#see the trend is removed. Then, check ACF and PACF.
adfTest(d_ts_data_t, lags=1, type="nc")#mean is very close to 0
#We have enough evidence to conclude that differenced series are stationary.
out2<-HEGY.test(wts=d_ts_data_t, itsd=c(1,0,0), regvar=0, selectlags=list(mode="signf", Pmax=NULL))
out2$stats
#there is regular unit root anymore
nsdiffs(d_ts_data_t)
#nsdiffs uses also ocsbtest and said no need to  seasonal diff
p3<-ggAcf(d_ts_data_t,lag.max=30)
p4<-ggPacf(d_ts_data_t,lag.max=30)
grid.arrange(p3,p4,nrow=1)
#ARIMA(2,1,2)
eacf(d_ts_data_t)
#ARIMA(1,1,1)
library(caschrono)
armaselect(d_ts_data_t)
#it is said that  ARIMA (0,1,0) is the suggested model.
# Assuming we have a list of candidate ARIMA models
candidate_models <- list(
  Arima(train, order = c(2, 1, 2)),
  Arima(train, order = c(1, 1, 1)),
  Arima(train, order = c(0, 1, 0)),
  Arima(train, order = c(2, 1, 0)),
  Arima(train, order = c(1, 1, 0)),
  Arima(train, order = c(0, 1, 1)),
  Arima(train, order = c(0, 1, 3)),
  Arima(train, order = c(3, 1, 0)),
  Arima(train, order = c(1, 1, 2)),
  Arima(train, order = c(2, 1, 1)),
  Arima(train, order = c(3, 1, 3)),
  Arima(train, order = c(4, 1, 2)),
  Arima(train, order = c(2, 1, 4))
  
  
)

sorted_indices <- order(sapply(candidate_models, function(model) model$aic))
# Print the three models with the lowest AIC values
for (i in seq_len(min(3, length(sorted_indices)))) {
  index <- sorted_indices[i]
  model <- candidate_models[[index]]
  order_values <- paste("order", model$arma[1], 1, model$arma[2], sep = " ")
  aic_value <- model$aic
  cat("Model", i, order_values, "AIC:", aic_value, "\n")
}
#Model 4 1 2 has smallest AIC with -14.93115 
Arima(train, order = c(3, 1, 0))
abs(0.1126  /0.0597   )#so it is not significant
Arima(train, order = c(2, 1, 1))
abs(-0.1298      /0.0692     )#so it is not significant
Arima(train, order = c(0, 1, 3))
abs(-0.1341/0.0596)#so it is not significant
Arima(train, order = c(2, 1, 2))#so it is not significant

Arima(train, order = c(4, 1, 2))
abs(-0.2392        /0.0589        )#so it is significant
abs(-0.7314      /0.1507     )#so it is significant
Arima(train, order = c(2, 1, 4))
abs(-0.6997          /0.1430           )#so it is significant
abs(-0.2388      /0.0586     )#so it is significant

fit1 <- Arima(train, order = c(4, 1, 2))
fit1
fit2 <- Arima(train, order = c(2, 1, 4))
fit2

#fit1 has lower AIC and BIC
r1 <- resid(fit1)
autoplot(r1)+geom_line(y=0)+theme_minimal()+ggtitle("Plot of The Residuals")
ggAcf(r1, main = "ACF of Residuals")
ggPacf(r1, main = "PACF of Residuals")
grid.arrange(ggAcf(r1, main = "ACF of Residuals")
,ggPacf(r1, main = "PACF of Residuals")
,nrow=1)

#are scattered around zero and it can be interpreted as zero mean.
ggplot(r1, aes(sample = r1)) +stat_qq()+geom_qq_line()+ggtitle("QQ Plot of the Residuals")+theme_minimal()
ggplot(r1,aes(x=r1))+geom_histogram(bins=20)+geom_density()+ggtitle("Histogram of Residuals")+theme_minimal()
shapiro.test(r1)
jarque.bera.test(r1)
library(vars)
#residuals are not distributed normal
Box.test(r1,lag=15,type = c("Ljung-Box"))
#Since p value is greater than 0.05, the residuals of the model are uncorrelated.
r2 <- resid(fit2)
autoplot(r2)+geom_line(y=0)+theme_minimal()+ggtitle("Plot of The Residuals")
ggAcf(r2, main = "ACF of Residuals")
ggPacf(r2, main = "PACF of Residuals")
grid.arrange(ggAcf(r2, main = "ACF of Residuals")
             ,ggPacf(r2, main = "PACF of Residuals")
             ,nrow=1)
#are scattered around zero and it can be interpreted as zero mean.
ggplot(r2, aes(sample = r2)) +stat_qq()+geom_qq_line()+ggtitle("QQ Plot of the Residuals")+theme_minimal()
ggplot(r2,aes(x=r2))+geom_histogram(bins=20)+geom_density()+ggtitle("Histogram of Residuals")+theme_minimal()
grid.arrange(ggplot(r2, aes(sample = r2)) +stat_qq()+geom_qq_line()+ggtitle("QQ Plot of the Residuals")+theme_minimal()

             ,ggplot(r2,aes(x=r2))+geom_histogram(bins=20)+geom_density()+ggtitle("Histogram of Residuals")+theme_minimal()
             ,nrow=1)
shapiro.test(r2)
jarque.bera.test(r2)
#residuals are not distributed normal
Box.test(r2,lag=15,type = c("Ljung-Box"))
#Since p value is greater than 0.05, the residuals of the model are uncorrelated.

library(rugarch)
library(rmgarch)
library(FinTS)
library(lmtest)
m = lm(r1 ~ 1+zlag(r1))
bgtest(m,order=15)
m = lm(r2 ~ 1+zlag(r2))
bgtest(m,order=15)

#Since p value is greater than 0.05, we have 95% confident that the residuals of the model are uncorrelated

p1<-ggAcf(as.vector(r1^2),main="ACF of Squared Residuals")+theme_minimal()
p2<-ggPacf(as.vector(r1^2),main="PACF of Squared Residuals")+theme_minimal() # homoscedasticity check
grid.arrange(p1,p2,ncol=2)
p1<-ggAcf(as.vector(r2^2),main="ACF of Squared Residuals")+theme_minimal()
p2<-ggPacf(as.vector(r2^2),main="PACF of Squared Residuals")+theme_minimal() # homoscedasticity check
grid.arrange(p1,p2,ncol=2)

ArchTest(r1)
ArchTest(r2)
#H0: Residuals exhibits no ARCH effects.
#H1: ARCH(lag) effects are present.
#Since p values is greather than 0.05 , we can not reject H0. Therefore, we can conclude there is no ARCH effects.


#now I will select arima(2,1,4)
#because it has low aic and parameter
test_arima <- accuracy(forecast(fit2,h=36)$mean,test)
summary(fit2)
train_arima <- accuracy(fit2)
accuracy(forecast(fit2,h=36)$mean,test)
train_arima <- accuracy(fit2)

autoplot(ts(df$Close, start = 1998, frequency = 12)) +
  autolayer(fitted(fit2), series = "Fitted", color = "red") +
  autolayer(forecast(fit2,h=36)$mean, series = "Forecast", color = "red") +
  autolayer(arima_forecast$upper[,2], series = "Upper PI", color = "red") +
  autolayer(arima_forecast$lower[,2], series = "Lower PI", color = "red") +
  geom_vline(xintercept = 2021, linetype = "dashed", color = "red") +
  ylab("Yield") + xlab("Year") +
  theme_minimal()
arima_forecast <- forecast(fit2,h=36)
arima_forecast$lower
#ets
library(TSA)
library(forecast)
ets_fit <- ets(train,model="ZZZ")
accuracy(forecast(ets_fit,h=36)$mean,test)
# identify optimal alpha parameter
beta <- seq(.0001, .5, by = .001)
RMSE <- NA
for(i in seq_along(beta)) {
  fit <- holt(train,
              beta = beta[i],
              h = 36)
  RMSE[i] <- accuracy(fit,
                      test)[2,2]
}
#ets not selected because holt perform better
# convert to a data frame and
# identify min alpha value
beta.fit <- data_frame(beta, RMSE)
beta.fit %>% filter(RMSE==min(RMSE))
holt.goog <- holt(train, h = 36,beta = 0.221)
summary(holt.goog)
test_tbats <- accuracy(holt.goog$mean, test)
train_tbats <- accuracy(holt.goog)
shapiro.test(resid(holt.goog))
# Create a time index for the entire dataset
autoplot(ts(df$Close, start = 1998, frequency = 12)) +
  autolayer(fitted(holt.goog), series = "Fitted", color = "red") +
  autolayer(holt.goog$mean, series = "Forecast", color = "red")+
  autolayer(holt.goog$upper[,2], series = "Upper PI",linetype = "dashed", color = "blue")+
  autolayer(holt.goog$lower[,2], series = "Lower PI",linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 2021, linetype = "dashed", color = "red") +
  ylab("Yield") + xlab("Year") +
  theme_minimal()
#tbats

tbatsmodel<-tbats(train)
shapiro.test(resid(tbatsmodel))
autoplot(train,main="TS plot of Train with TBATS Fitted") +autolayer(fitted(tbatsmodel), series="Fitted") +theme_minimal()
tbats_forecast<-forecast(tbatsmodel,h=36)

tbatsmodel
autoplot(tbats_forecast)+theme_minimal()
autoplot(tbats_forecast)+autolayer(test,series="actual",color="red")+theme_minimal()
t_t_tbats<- accuracy(tbats_forecast,test)
tbats_forecast$upper
autoplot(ts(df$Close,start = 1998,frequency = 12)) +
  autolayer(fitted(tbats_forecast), series = "Fitted", color = "red") +
  autolayer(tbats_forecast, series = "Forecast", color = "red")+
  autolayer(tbats_forecast$upper[,2], series = "Upper PI",linetype = "dashed", color = "blue")+
  autolayer(tbats_forecast$lower[,2], series = "Lower PI",linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 2021, linetype = "dashed", color = "red") +
  ylab("Yield") + xlab("Year") +
  theme_minimal()


#ntar
sizes <- c(1:20)
repeats <- c(1:20)
decays <- c(.0001:0.0005,0.001:0.005,0.01:0.05,0.1:0.5)
ş <- nnetar(train,size = 1,repeats = 2,decays=0)
ş$model
forecast(ş,h=36)
şş <- accuracy(forecast(ş,h=36),test)
şş[2,2]
sizes <- c(1,seq(2,10,by=2))
repeats <- c(1,seq(2,10,by=2))
decays <- c(0,0.001,0.01,0.1)
# Initialize variables to store best model and accuracy
best_model <- NULL
best_accuracy <- Inf
# Loop over parameters
for (size in sizes) {
  for (repeat_val in repeats) {
    for (decay_val in decays) {
      model <- nnetar(train, size = size, repeats = repeat_val, decays = decay_val,
                      lambda=NULL,PI=TRUE)
      cat(decay_val,repeat_val,size)
      forecast_result <- forecast(model, h = 36)
      
      acc <- accuracy(forecast_result, test)
      
      if (acc[2,2] < best_accuracy) {
        best_accuracy <- acc[2,2]
        best_model <- model
        print("best")
      }
    }
  }
}
bmn <- forecast(best_model, h = 36)
best_model
nnetar(train, size = 1, repeats = 1, decays = 0, lambda=NULL,PI=TRUE)
fitted(forecast(best_model, h = 36))
accuracy(forecast(best_model, h = 36),test)
shapiro.test(resid(best_model))
tss <- ts(df$Close,start = 1998,frequency = 12)
autoplot(tss) +
  autolayer(fitted(best_model), series = "Fitted", color = "red") +
  autolayer(forecast(best_model,PI=TRUE, h = 36)$mean, series = "Forecast", color = "red")+
  autolayer(forecast(best_model,PI=TRUE, h = 36)$upper[,2], series = "Upper PI",linetype = "dashed", color = "blue")+
  autolayer(forecast(best_model,PI=TRUE, h = 36)$lower[,2], series = "Lower PI",linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 2021, linetype = "dashed", color = "red") +
  ylab("Yield") + xlab("Year") +
  theme_minimal()
#prophet
library(prophet)

ds<-c(seq(as.Date("1998/01/01"),as.Date("2020/12/01"),by="month"))
dk<-data.frame(ds,y=head(df$Close,-36))
seasonality.prior.scales <- seq(5, 50, by = 5)
changepoint.prior.scales <- seq(0.1, 1, by = 0.1)
changepoint.range <- seq(.25,1,.25)
best_acc <- Inf
best_params <- c()

# Nested loops for parameter tuning
for (seasonality_scale in seasonality.prior.scales) {
  for (changepoint_scale in changepoint.prior.scales) {
    for (changepoint.rangef in changepoint.range) {
    # Fit prophet model
    pro <- prophet(dk,
                   seasonality.prior.scale = seasonality_scale,
                   changepoint.prior.scale = changepoint_scale,
                   changepoint.range=changepoint.rangef)
    future <- make_future_dataframe(pro, periods = 36, freq = 'month')
    forecast <- predict(pro, future)
    # Calculate accuracy
    acc <- accuracy(tail(forecast$yhat, 36), test)
    # Check if the current result is better than the best
    if (acc[2] < best_acc) {  # Assuming acc[2] is RMSE
      best_acc <- acc[2]
      best_params <- c(seasonality_scale, changepoint_scale,changepoint.range)
    }
  }
  }
}

# Print the best result
cat("Best RMSE:", best_acc, "\n")
cat(best_params)

pro <- prophet(dk,seasonality.prior.scale=15,changepoint.prior.scale=0.1,changepoint.range=.75)
summary(pro)
future<-make_future_dataframe(pro,periods = 36,freq='month')
forecast <- predict(pro, future)
accuracy(tail(forecast$yhat,36),test)
shapiro.test(pro$history$y_scaled)
df <- read.csv("10year.csv")
testt <- tail(df,36)
trainn <- head(df,-36)

test_data_prophet <- testt %>%
  rename(ds = Date, y = Close) %>%
  mutate(ds = as.Date(ds))

forecast <- forecast %>%
  dplyr::select(ds, yhat, yhat_lower, yhat_upper) %>%
  dplyr::left_join(test_data_prophet, by = "ds")

ggplot(forecast, aes(x = ds)) +
  geom_line(aes(y = yhat), color = "black", size = 1) +  
  geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper), fill = "gray", alpha = 0.5) +
  geom_line(aes(y = y), color = "red", size = 1) +  
  labs(title = "Prophet vs Actual Data", x = "Date", y = "Yield") +
  theme_minimal()


test_arima <- accuracy(forecast(fit3,h=36)$mean,test)
train_arima <- accuracy(fit3)

test_ets <- accuracy(holt.goog$mean, test)
train_ets <- accuracy(holt.goog)
t_t_tbats<- accuracy(tbats_forecast,test)

#arima
round(train_arima,3)
round(test_arima,3)
#ets
round(accuracy(holt.goog$mean, test),3)
round(accuracy(holt.goog),3)
#tbats
round(accuracy(tbats_forecast,test),3)
#nnetar
round(accuracy(forecast(best_model, h = 36),test),3)
#prophet
round(accuracy(tail(forecast$yhat,36),test),3)
round(accuracy(forecast$yhat,train),3)
