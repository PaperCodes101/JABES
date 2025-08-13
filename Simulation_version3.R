#install.packages("INLA")

library(MASS)
library(spdep)
library(INLA)
library(dplyr)
library(ggplot2)
#library(glmmTMB)

#-----------------------------
# 1. Simulation parameters
#-----------------------------
#set.seed(123)

merged_data = read.csv("~/Research Spring 2025/Research Paper 01 Practice codes/January to March Data/merged_data.csv")
colnames(merged_data)

merged_data = merged_data %>%
  mutate(County_Name.x = ifelse(FIP == "48033", "Borden", County_Name.x))

merged_data = merged_data %>%
  mutate(County_Name.x = ifelse(FIP == "48123", "De Witt", County_Name.x))

merged_data = merged_data[order(merged_data$FIP),]

merged_data <- merged_data[ , -c(8, 9, 10, 12, 13, 14, 16)]

n_counties <- 254  # number of counties

# Updating outliers
quantile(merged_data$Measles_Count, 0.997)

q_high <- quantile(merged_data$Measles_Count, 0.997)
q_high_new <- as.numeric(ceiling(q_high))

merged_data$Measles_Count[merged_data$Measles_Count > q_high] <- q_high_new
print(merged_data$Measles_Count)

summary(merged_data$MMR_Avg)
any(merged_data$MMR_Avg == 0)

# Updating outliers
quantile(merged_data$MMR_Avg, 0.5)

q_high_MMR <- quantile(merged_data$MMR_Avg, 0.5)

merged_data$MMR_Avg[merged_data$MMR_Avg ==0] <- q_high_MMR
print(merged_data$MMR_Avg)
any(merged_data$MMR_Avg == 0)

county   = merged_data$County_Name.x
temp_obs = merged_data$Average_Temperature  
prec_obs = merged_data$Average_Precipitation    
vacc_obs = merged_data$MMR_Avg
expected_obs = merged_data$Expected_Cases
Y_obs = merged_data$Measles_Count

#----------------------------------------------#

M1_waic<-c()
M2_waic<-c()
M3_waic<-c()
M4_waic<-c()
M5_waic<-c()
M6_waic<-c()

B = 2000

for (i in 1:B) {  #i = 1
  
  simulate_dataset <- function() {
    
    # Covariates
    X1 <- rnorm(n_counties, mean = vacc_obs, sd = 0.05)       # MMR
    X2 <- rnorm(n_counties, mean = temp_obs, sd = 0.05)       # Temp
    X3 <- rnorm(n_counties, mean = prec_obs, sd = 0.05)       # Prec
    X4 <- expected_obs   #expected cases
    
    lambda <- Y_obs
    
    gamma0 <- 0.730
    
    pi <- exp(gamma0) / (1 + exp(gamma0))
    
    # Zero-inflated Poisson counts
    z <- rbinom(n_counties, 1, pi)                              # Bernoulli: 1 = forced zero
    Y <- ifelse(z == 1, 0, rpois(n_counties, lambda))           # ZIP outcome
    
    data.frame(Y, X1, X2, X3, X4, FIPS = merged_data$FIP , county = 1:n_counties)
  }
  
  
  sim_data <- simulate_dataset()
  summary(sim_data)
  
  #hist(sim_data$Y, breaks = 50, main = "Simulated Y")
  
  shape_data = st_read("~/Research Spring 2025/Research Paper 01 Practice codes/January to March Data/Texas_County_Boundaries_Detailed.shp")
  shape_data = shape_data[order(shape_data$CNTY_FIPS), ]
  
  sim_data$FIPS <- as.character(sim_data$FIPS)
  shape_data$CNTY_FIPS <- as.character(shape_data$CNTY_FIPS)
  tx_map <- left_join(shape_data, sim_data, by = c("CNTY_FIPS" = "FIPS"))
  
  # ggplot(tx_map) +
  #  geom_sf(aes(fill = Y), color = "gray30", size = 0.2) +
  #  scale_fill_viridis_c(option = "plasma", name = "Simulated\nCount (Y)") +
  #  theme_minimal() +
  #  labs(title = "Simulated ZIP Disease Counts by Texas County",
  #       caption = "Counts simulated from ZIP + BYM model") +
  #  theme(axis.text = element_blank(),
  #        axis.ticks = element_blank(),
  #        panel.grid = element_blank())
  
  nb_1 <- poly2nb(shape_data)
  head(nb_1)
  
  nb2INLA("map_1.adj", nb_1)
  g_1 <- inla.read.graph(filename = "map_1.adj")
  
  # Model Fitting
  # Zip model with spatial effect bym_Model 04
  
  # ZIP with BYM spatial effect - Model 04
  formula1 <- Y ~ offset(log(X4)) + X1 + X2 + X3 + 
    f(county, model = "bym", graph = g_1)
  
  # ZIP with BYM spatial effect - Model 04
  res1 <- tryCatch({
    inla(
      formula1,
      family = "zeroinflatedpoisson1",
      data = sim_data,
      control.predictor = list(compute = TRUE, link = 1),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)
    )
  }, error = function(e) {
    cat("Model 04 (ZIP + BYM) failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  M4_waic[i] <- if (!is.null(res1)) res1$waic$waic else NA
  
  
  # Poisson with BYM spatial effect - Model 01
  formula2 <- Y ~ offset(log(X4)) + X1 + X2 + X3 + 
    f(county, model = "bym", graph = g_1)
  
  res2 <- tryCatch({
    inla(
      formula2,
      family = "poisson",
      data = sim_data,
      control.predictor = list(compute = TRUE, link = 1),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)
    )
  }, error = function(e) {
    cat("Model 01 (Poisson + BYM) failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  M1_waic[i] <- if (!is.null(res2)) res2$waic$waic else NA
  
  
  # ZIP with IID random effect - Model 06
  formula3 <- Y ~ offset(log(X4)) + X1 + X2 + X3 + 
    f(county, model = "iid")
  
  res3 <- tryCatch({
    inla(
      formula3,
      family = "zeroinflatedpoisson1",
      data = sim_data,
      control.predictor = list(compute = TRUE, link = 1),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)
    )
  }, error = function(e) {
    cat("Model 06 (ZIP + IID) failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  M6_waic[i] <- if (!is.null(res3)) res3$waic$waic else NA
  
  
  # Poisson with IID random effect - Model 03
  formula4 <- Y ~ offset(log(X4)) + X1 + X2 + X3 + 
    f(county, model = "iid")
  
  res4 <- tryCatch({
    inla(
      formula4,
      family = "poisson",
      data = sim_data,
      control.predictor = list(compute = TRUE, link = 1),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)
    )
  }, error = function(e) {
    cat("Model 03 (Poisson + IID) failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  M3_waic[i] <- if (!is.null(res4)) res4$waic$waic else NA
  
  
  
  # DEfine leroux model
  
  # Adjacency matrix W (symmetric, binary)
  W <- as(nb2mat(nb_1, style="B"), "Matrix")
  
  # Diagonal matrix D: number of neighbors per region
  D <- Diagonal(x = rowSums(W))
  
  # Intrinsic CAR precision matrix
  Q <- D - W
  
  # Identity matrix
  n <- nrow(Q)
  I <- Diagonal(n)
  
  # Matrix C = I - Q
  C <- I - Q
  
  
  # Poisson with spatial effect - Leroux (generic1) - Model 02
  formula5 <- Y ~ offset(log(X4)) + X1 + X2 + X3 + 
    f(county, model = "generic1", Cmatrix = C)
  
  res5 <- tryCatch({
    inla(
      formula5,
      family = "poisson",
      data = sim_data,
      control.predictor = list(compute = TRUE, link = 1),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)
    )
  }, error = function(e) {
    cat("Model 02 (Poisson + Leroux) failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  M2_waic[i] <- if (!is.null(res5)) res5$waic$waic else NA
  
  
  # ZIP with spatial effect - Leroux (generic1) - Model 05
  formula6 <- Y ~ offset(log(X4)) + X1 + X2 + X3 + 
    f(county, model = "generic1", Cmatrix = C)
  
  res6 <- tryCatch({
    inla(
      formula6,
      family = "zeroinflatedpoisson1",
      data = sim_data,
      control.predictor = list(compute = TRUE, link = 1),
      control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE)
    )
  }, error = function(e) {
    cat("Model 05 (ZIP + Leroux) failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  M5_waic[i] <- if (!is.null(res6)) res6$waic$waic else NA
  
  
  print(i)
  
}

print(M1_waic)
print(M4_waic)
print(M2_waic)
print(M3_waic)
print(M5_waic)
print(M6_waic)

M1B_waic<-na.omit(M1_waic)
M2B_waic<-na.omit(M2_waic)
M3B_waic<-na.omit(M3_waic)
M4B_waic<-na.omit(M4_waic)
M5B_waic<-na.omit(M5_waic)
M6B_waic<-na.omit(M6_waic)

summary(M1B_waic)
summary(M2B_waic)
summary(M3B_waic)
summary(M4B_waic)
summary(M5B_waic)
summary(M6B_waic)

quantile(M1B_waic, 0.99)
quantile(M2B_waic, 0.99)
quantile(M3B_waic, 0.99)
quantile(M4B_waic, 0.99)
quantile(M5B_waic, 0.99)
quantile(M6B_waic, 0.99)

q_high_M1 <- quantile(M1B_waic, 0.99)
M1B_waic[M1B_waic > q_high_M1] <- q_high_M1


q_high_M2<- quantile(M2B_waic, 0.99)
M2B_waic[M2B_waic > q_high_M2] <- q_high_M2

q_high_M3<- quantile(M3B_waic, 0.99)
M3B_waic[M3B_waic > q_high_M3] <- q_high_M3

q_high_M4<- quantile(M4B_waic, 0.99)
M4B_waic[M4B_waic > q_high_M4] <- q_high_M4

q_high_M5<- quantile(M5B_waic, 0.99)
M5B_waic[M5B_waic > q_high_M5] <- q_high_M5

q_high_M6<- quantile(M6B_waic, 0.99)
M6B_waic[M6B_waic > q_high_M6] <- q_high_M6

AM1_waic<-mean(M1B_waic)
AM2_waic<-mean(M2B_waic)
AM3_waic<-mean(M3B_waic)
AM4_waic<-mean(M4B_waic)
AM5_waic<-mean(M5B_waic)
AM6_waic<-mean(M6B_waic)

SM1_waic<-sd(M1B_waic)
SM2_waic<-sd(M2B_waic)
SM3_waic<-sd(M3B_waic)
SM4_waic<-sd(M4B_waic)
SM5_waic<-sd(M5B_waic)
SM6_waic<-sd(M6B_waic)

data.frame(AM1_waic, AM2_waic, AM3_waic, AM4_waic, AM5_waic, AM6_waic)
data.frame(SM1_waic, SM2_waic, SM3_waic, SM4_waic, SM5_waic, SM6_waic)
