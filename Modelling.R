#install.packages("sn")

library(INLA)
library(dplyr)
library(sf)
library(leaflet)
library(spdep)
library(sn)
library(coda)
library(ggplot2)
library(scales)

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
colnames(shape_data)

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

####----- plot SIR values in texas map------ ###

# Create custom breaks focusing on low-end values
breaks_SIR <- c(0, 1, 2.5, 5, 10, 25, 50, 100, max(map_1$SIR, na.rm = TRUE))

# Corresponding colors (you can tweak intensity to match your aesthetic)
color_scale_SIR <- scale_fill_gradientn(
  colors = c(
    "#FFFFCC", # pale yellow
    "#FFE0B2", # light orange
    "#FDBE85", # soft orange
    "#FD8D3C", # medium orange
    "#FC4E2A", # orange-red
    "#E31A1C", # red
    "#BD0026", # dark red
    "#800026", # deep red
    "#4D004B", # purplish
    "#2A0044"  # dark purple
  ),
  values = rescale(breaks_SIR),  # rescale to 0–1
  breaks = breaks_SIR,
  name = "SIR"
)

# Final plot
ggplot(map_1) +
  geom_sf(aes(fill = SIR), color = "black", size = 0.1) +
  color_scale_SIR +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(2.0, "cm"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

# Choose only key breakpoints for legend
legend_breaks <- c(0, 50, 100, 150, 200, 250, 300)
legend_labels <- c("0", "50", "100", "150", "200", "250", "300")

color_scale_SIR <- scale_fill_gradientn(
  colors = c(
    "#FFFFCC", "#FFE0B2", "#FDBE85", "#FD8D3C",
    "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#4D004B"
  ),
  values = rescale(breaks_SIR),
  breaks = legend_breaks,
  labels = legend_labels,
  name = "SIR"
)

theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 10, face = "bold"),
  legend.key.width = unit(1.5, "cm")
)

ggplot(map_1) +
  geom_sf(aes(fill = SIR), color = "black", size = 0.1) +
  scale_fill_gradientn(
    colors = c(
      "#FFFFCC", "#FFE0B2", "#FDBE85", "#FD8D3C",
      "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#4D004B"
    ),
    values = rescale(breaks_SIR),
    breaks = legend_breaks,
    labels = legend_labels,
    name = "SIR"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.75, "cm"),
    legend.key.height = unit(1.0, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )
#==================================================================#



#### Neighbours #####

map_1$MMR_Avg = map_1$MMR_Avg * 100


nb_1 <- poly2nb(map_1)
head(nb_1)

nb2INLA("map_1.adj", nb_1)
g_1 <- inla.read.graph(filename = "map_1.adj")

# Define spatial structured Random effect
map_1$idareau_1 <- 1:nrow(map_1)

stopifnot(nrow(map_1) == g_1$n)   # same number of areas

#### Poisson with spatial effect ####

formula33 <- Measles_Count ~  offset(log(Expected_Cases)) + MMR_Avg + Average_Temperature + Average_Precipitation +
  f(idareau_1, model = "bym", graph = g_1)

res33 <- inla(formula33,
              family = "poisson", data = map_1,
              control.predictor = list(compute = TRUE),
              control.compute = list(return.marginals.predictor = TRUE, dic = TRUE, waic = TRUE)
)

summary(res33)

#### Poisson with IID ####

formula44 <- Measles_Count ~ offset(log(Expected_Cases)) + MMR_Avg + Average_Temperature + Average_Precipitation +
  f(idareau_1, model = "iid")

res44 <- inla(formula44,
              family = "poisson", data = map_1,
              control.predictor = list(compute = TRUE),
              control.compute = list(return.marginals.predictor = TRUE, dic = TRUE, waic = TRUE)
)

summary(res44)


### Zip1 with spatial effect ####

formula55 <- Measles_Count ~ offset(log(Expected_Cases)) + MMR_Avg + Average_Temperature + Average_Precipitation +
  f(idareau_1, model = "bym", graph = g_1)

res55 <- inla(formula55,
              family = "zeroinflatedpoisson1", data = map_1,
              control.predictor = list(compute = TRUE),
              control.compute = list(return.marginals.predictor = TRUE, config = TRUE, dic = TRUE, waic = TRUE)
)

summary(res55)

### zip1 with IID ###

formula66 <- Measles_Count ~ offset(log(Expected_Cases)) + MMR_Avg + Average_Temperature + Average_Precipitation +
  f(idareau_1, model = "iid")

res66 <- inla(formula66,
              family = "zeroinflatedpoisson1", data = map_1,
              control.predictor = list(compute = TRUE),
              control.compute = list(return.marginals.predictor = TRUE, cpo = TRUE, config = TRUE, dic = TRUE, waic = TRUE)
)

summary(res66)
