library(INLA)
library(caret)
library(dplyr)
library(sf)
library(leaflet)
library(spdep)
library(sn)
library(coda)
library(ggplot2)



merged_data = read.csv("~/Research Spring 2025/Research Paper 01 Practice codes/January to March Data/merged_data.csv")

merged_data = merged_data %>%
  mutate(County_Name.x = ifelse(FIP == "48033", "Borden", County_Name.x))

merged_data = merged_data %>%
  mutate(County_Name.x = ifelse(FIP == "48123", "De Witt", County_Name.x))

# Updating outliers
quantile(merged_data$Measles_Count, 0.997)

q_high <- quantile(merged_data$Measles_Count, 0.997)
q_high_new <- as.numeric(ceiling(q_high))

merged_data$Measles_Count[merged_data$Measles_Count > q_high] <- q_high_new
print(merged_data$Measles_Count)


merged_data = merged_data %>%
  mutate(MMR_Avg = ifelse(County_Name.x == "Borde", 0.8750, MMR_Avg))

#merged_data = merged_data %>%
#  mutate(MMR_Avg = ifelse(County_Name.x == "Crane", 0.9767, MMR_Avg))

#merged_data = merged_data %>%
#  mutate(MMR_Avg = ifelse(County_Name.x == "Stonewall", 1.0, MMR_Avg))

summary(merged_data$MMR_Avg)
any(merged_data$MMR_Avg == 0)

# Updating outliers
quantile(merged_data$MMR_Avg, 0.5)

q_high_MMR <- quantile(merged_data$MMR_Avg, 0.5)

merged_data$MMR_Avg[merged_data$MMR_Avg ==0] <- q_high_MMR
print(merged_data$MMR_Avg)
any(merged_data$MMR_Avg == 0)

class(merged_data)
shape_data = st_read("~/Research Spring 2025/Research Paper 01 Practice codes/January to March Data/Texas_County_Boundaries_Detailed.shp")

shape_data = shape_data[order(shape_data$CNTY_FIPS), ]
merged_data = merged_data[order(merged_data$FIP),]

#To check the raws are same
names1 <- merged_data$County_Name.x
names2 <- shape_data$CNTY_NM

identical(names1, names2) # [1] TRUE

# Calculate centroids
centroids_1 <- st_centroid(shape_data)
geo_1 <- data.frame(
  county.names = shape_data$CNTY_NM,  # replace NAME with your county name column
  x = st_coordinates(centroids_1)[, 1],
  y = st_coordinates(centroids_1)[, 2]
)

#To check the raws are same
names3 <- merged_data$County_Name.x
names4 <- geo_1$county.names

identical(names3, names4) # TRUE


colnames(merged_data)

merged_data <- merged_data %>%
  rename(county.names = County_Name.x)

#To check the raws are same
names5 <- merged_data$county.names
names6 <- geo_1$county.names

identical(names5, names6) # TRUE

data_1 <- merge(geo_1, merged_data, by = "county.names", sort = FALSE)

names7 <- merged_data$county.names
names8 <- data_1$county.names

identical(names7, names8) #TRUE


spatial_polygon_1 <- as(shape_data, "Spatial")
class(spatial_polygon_1)

polygon_1 <- shape_data

#To check the raws are same
names9 <- data_1$county.names
names10 <- spatial_polygon_1@data$CNTY_NM
names11 <- polygon_1$CNTY_NM

identical(names9, names10) # TRUE
identical(names9, names11) #TRUE

New_data_1 <- list(
  geo = geo_1,
  data = data_1,
  spatial.polygon = spatial_polygon_1,
  polygon = polygon_1
)

class(New_data_1)
names(New_data_1)

head(New_data_1$data)

map_1 <- New_data_1$spatial.polygon
plot(map_1)

map_1 <- st_as_sf(map_1)
map_1 <- st_transform(map_1, crs = 4326)

# Check New_data1 and map_1 raws are in the same order to check counties are not messed up
names12 <- New_data_1$data$county.names
names13 <- map_1$CNTY_NM

identical(names12, names13) # TRUE


d_1 <- New_data_1$data[, c("county.names", "Measles_Count", "Expected_Cases", "MMR_Avg", "Average_Temperature", "Average_Precipitation")]

d_1$SIR <- d_1$Measles_Count / d_1$Expected_Cases

head(d_1)

names14 <- d_1$county.names
names15 <- map_1$CNTY_NM

identical(names14, names15) # TRUE

max(d_1$SIR)
min(d_1$SIR)

# Check if the county names are exactly equal and in the same order before the combination
all(map_1$CNTY_NM == d_1$county.names) #TRUE

map_1 <- cbind(map_1, d_1)

names16 <- d_1$county.names
names17 <- map_1$CNTY_NM

identical(names16, names17) #TRUE


head(map_1)

#### Neighbours #####

map_1$MMR_Avg = map_1$MMR_Avg * 100


nb_1 <- poly2nb(map_1)
head(nb_1)

nb2INLA("map_1.adj", nb_1)
g_1 <- inla.read.graph(filename = "map_1.adj")

# Define spatial structured Random effect
map_1$idareau_1 <- 1:nrow(map_1)
data_1$idareau_1 <- 1:nrow(data_1)

names18 = map_1$CNTY_NM
names19 = data_1$county.names

identical(names18, names19) # TRUE

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


# ---- Step 1: Create 5 folds ----
k <- 5
folds <- createFolds(data_1$Measles_Count, k = k, list = TRUE)

# ---- Step 2: Initialize storage for results ----

rmse_poisson_spatial_bym <- numeric(k)
mae_poisson_spatial_bym <- numeric(k)

rmse_poisson_spatial_le <- numeric(k)
mae_poisson_spatial_le <- numeric(k)

rmse_poisson_iid <- numeric(k)
mae_poisson_iid <- numeric(k)

rmse_zip_spatial_bym <- numeric(k)
mae_zip_spatial_bym <- numeric(k)

rmse_zip_spatial_le <- numeric(k)
mae_zip_spatial_le <- numeric(k)

rmse_zip_iid <- numeric(k)
mae_zip_iid <- numeric(k)


# ---- Step 3: Run cross-validation ----
for (i in 1:k) {
  cat("\nFold", i, "\n")
  
  test_idx <- folds[[i]]
  train_idx <- setdiff(seq_len(nrow(data_1)), test_idx)
  
  train_data <- data_1[train_idx, ]
  test_data <- data_1[test_idx, ]
  
  # --------------------------
  # Poisson Model
  # --------------------------
  formula_poisson_spatial_bym <- Measles_Count ~ offset(log(Expected_Cases)) + MMR_Avg + Average_Temperature + Average_Precipitation +
    f(idareau_1, model = "bym", graph = g_1)
  
  result_poisson_spatial_bym <- inla(
    formula_poisson_spatial_bym,
    data = data_1,
    family = "poisson",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE)
  )
  
  preds_poisson_bym <- result_poisson_spatial_bym$summary.fitted.values$mean[test_idx]
  obs_poisson_bym <- data_1$Measles_Count[test_idx]
  
  rmse_poisson_spatial_bym[i] <- sqrt(mean((obs_poisson_bym - preds_poisson_bym)^2))
  mae_poisson_spatial_bym[i] <- mean(abs(obs_poisson_bym - preds_poisson_bym))
  
  
  # --------------------------
  # Poisson Model
  # --------------------------
  formula_poisson_spatial_le <- Measles_Count ~ offset(log(Expected_Cases)) + MMR_Avg + Average_Temperature + Average_Precipitation +
    f(idareau_1,  model = "generic1", Cmatrix = C)
  
  result_poisson_spatial_le <- inla(
    formula_poisson_spatial_le,
    data = data_1,
    family = "poisson",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE)
  )
  
  preds_poisson_spatial_le <- result_poisson_spatial_le$summary.fitted.values$mean[test_idx]
  obs_poisson_spatial_le <- data_1$Measles_Count[test_idx]
  
  rmse_poisson_spatial_le[i] <- sqrt(mean((obs_poisson_spatial_le - preds_poisson_spatial_le)^2))
  mae_poisson_spatial_le[i] <- mean(abs(obs_poisson_spatial_le - preds_poisson_spatial_le))
  
  
  # --------------------------
  # Poisson Model iid
  # --------------------------
  formula_poisson_iid <- Measles_Count ~ offset(log(Expected_Cases)) + MMR_Avg + Average_Temperature + Average_Precipitation +
    f(idareau_1, model = "iid")
  
  result_poisson_iid <- inla(
    formula_poisson_iid,
    data = data_1,
    family = "poisson",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE)
  )
  
  preds_poisson_iid <- result_poisson_iid$summary.fitted.values$mean[test_idx]
  obs_poisson_iid <- data_1$Measles_Count[test_idx]
  
  rmse_poisson_iid[i] <- sqrt(mean((obs_poisson_iid - preds_poisson_iid)^2))
  mae_poisson_iid[i] <- mean(abs(obs_poisson_iid - preds_poisson_iid))
  
  
  # --------------------------
  # zip Model spatial
  # --------------------------
  formula_zip_spatial_bym <- Measles_Count ~ offset(log(Expected_Cases)) + MMR_Avg + Average_Temperature + Average_Precipitation +
    f(idareau_1, model = "bym", graph = g_1)
  
  result_zip_spatial_bym <- inla(
    formula_zip_spatial_bym,
    data = data_1,
    family = "zeroinflatedpoisson1",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE)
  )
  
  preds_zip_bym <- result_zip_spatial_bym$summary.fitted.values$mean[test_idx]
  obs_zip_bym <- data_1$Measles_Count[test_idx]
  
  rmse_zip_spatial_bym[i] <- sqrt(mean((obs_zip_bym - preds_zip_bym)^2))
  mae_zip_spatial_bym[i] <- mean(abs(obs_zip_bym - preds_zip_bym))
  
  
  # --------------------------
  # zip Model spatial
  # --------------------------
  formula_zip_spatial_le <- Measles_Count ~ offset(log(Expected_Cases)) + MMR_Avg + Average_Temperature + Average_Precipitation +
    f(idareau_1,  model = "generic1", Cmatrix = C)
  
  result_zip_spatial_le <- inla(
    formula_zip_spatial_le,
    data = data_1,
    family = "zeroinflatedpoisson1",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE)
  )
  
  preds_zip_le <- result_zip_spatial_le$summary.fitted.values$mean[test_idx]
  obs_zip_le <- data_1$Measles_Count[test_idx]
  
  rmse_zip_spatial_le[i] <- sqrt(mean((obs_zip_le - preds_zip_le)^2))
  mae_zip_spatial_le[i] <- mean(abs(obs_zip_le - preds_zip_le))
  
  
  
  # --------------------------
  # zip Model iid
  # --------------------------
  formula_zip_iid <- Measles_Count ~ offset(log(Expected_Cases)) + MMR_Avg + Average_Temperature + Average_Precipitation +
    f(idareau_1, model = "iid")
  
  result_zip_iid <- inla(
    formula_zip_iid,
    data = data_1,
    family = "zeroinflatedpoisson1",
    control.predictor = list(compute = TRUE),
    control.compute = list(dic = TRUE)
  )
  
  preds_zip_iid <- result_zip_iid$summary.fitted.values$mean[test_idx]
  obs_zip_iid <- data_1$Measles_Count[test_idx]
  
  rmse_zip_iid[i] <- sqrt(mean((obs_zip_iid - preds_zip_iid)^2))
  mae_zip_iid[i] <- mean(abs(obs_zip_iid - preds_zip_iid))
}
# ---- Step 4: Aggregate results ----
cat("\n--- Cross-Validation Results ---\n")
cat("Poisson Model with spatial effect:\n")
cat("  Mean RMSE:", mean(rmse_poisson_spatial_bym), "\n")
cat("  Mean MAE :", mean(mae_poisson_spatial_bym), "\n\n")

cat("\n--- Cross-Validation Results ---\n")
cat("Poisson Model with leroux:\n")
cat("  Mean RMSE:", mean(rmse_poisson_spatial_le), "\n")
cat("  Mean MAE :", mean(mae_poisson_spatial_le), "\n\n")


cat("\n--- Cross-Validation Results ---\n")
cat("Poisson Model iid:\n")
cat("  Mean RMSE:", mean(rmse_poisson_iid), "\n")
cat("  Mean MAE :", mean(mae_poisson_iid), "\n\n")

cat("\n--- Cross-Validation Results ---\n")
cat("zip Model with spatial:\n")
cat("  Mean RMSE:", mean(rmse_zip_spatial_bym), "\n")
cat("  Mean MAE :", mean(mae_zip_spatial_bym), "\n\n")


cat("zip Model with leroux:\n")
cat("  Mean RMSE:", mean(rmse_zip_spatial_le), "\n")
cat("  Mean MAE :", mean(mae_zip_spatial_le), "\n\n")

cat("\n--- Cross-Validation Results ---\n")
cat("zip Model with iid:\n")
cat("  Mean RMSE:", mean(rmse_zip_iid), "\n")
cat("  Mean MAE :", mean(mae_zip_iid), "\n\n")

