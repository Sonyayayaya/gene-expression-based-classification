library(tseries)
library(forecast)
require(caret)
#install.packages("MLmetrics")
library(MLmetrics)
library(boot)

luad_expMat <- read.csv('C:.../genes.csv')

# Check the type of data
# str(luad_expMat)


# Selecting the line for the Dickey-Fowler test
# Select rows of the transposed dataset - patient's gene sets
new_vector <- unlist(luad_expMat[12, 2:35], use.names = FALSE)

# Create a raw data plot
plot.ts(new_vector, xlab = "Index", ylab = "Gene expression")

# There are two hypotheses for each vector:
# H0: A time series is non-stationary. 
# In other words, it has some time-dependent structure 
# and does not have a constant dispersion in time.
# H1: The time series is stationary.

# If the p-value from the test is less than some 
# significance level (e.g., α = 0.05), 
# then we can reject the null hypothesis and conclude 
# that the time series is stationary.
adf.test(new_vector)

# Make the data into a stationary series
stationary_vector = diff(log(new_vector))

# Diff series plot
plot.ts(stationary_vector, xlab = "Index", ylab = "Gene expression")

# Checking the stationarity of the data
adf.test(stationary_vector)

# Create a function to calculate MAPE
mape_func <- function(data, indices) {
  
  # The factor of dividing the dataset into training and test data
  smp_size <- floor(0.75 * length(data))
  
  # Initialize the pseudo-random number generator
  set.seed(12)
  
  # Create a range of numbers to divide the dataset
  indices <- sample(seq_len(length(data)), size = smp_size)
  
  # Split the dataset
  train_data <- data[indices]
  test_data <- data[-indices]
  
  # Create a SARIMAX model
  model_sarimax <- auto.arima(train_data, seasonal = TRUE, xreg = NULL)
  
  # Creating forecasts
  forecast_sarimax <- forecast(model_sarimax, h = length(test_data))
  
  # MAPE metric calculation
  mape <- MAPE(forecast_sarimax$mean, test_data)
  return(mape)
}

# Perform bootstrapping
bootstrapped_mape <- boot(data = stationary_vector, statistic = mape_func, R = 100)

# Confidence intervals for the MAPE
ci <- boot.ci(boot.out = bootstrapped_mape, type = "basic", conf = 0.95)
print(ci)

# Print the model results
print(model_sarimax)