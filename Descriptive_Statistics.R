library(readr)
library(GGally)


#install.packages("MazamaSpatialPlots")

library(MazamaSpatialPlots) # to easily create thematic maps, particularly for US states and counties
library(usmap) # to easily create visually appealing maps of the United States
library(ggplot2) # use ggplot2 to add layer for visualization

merged_data = read.csv('~/Research Spring 2025/Research Paper 01 Practice codes/January to March Data/merged_data.csv')
colnames(merged_data)

merged_data = merged_data %>%
  mutate(County_Name.x = ifelse(FIP == "48033", "Borden", County_Name.x))

merged_data = merged_data %>%
  mutate(County_Name.x = ifelse(FIP == "48123", "De Witt", County_Name.x))


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

# Updating outliers
quantile(merged_data$Measles_Count, 0.997)

q_high <- quantile(merged_data$Measles_Count, 0.997)
q_high_new <- as.numeric(ceiling(q_high))

merged_data$Measles_Count[merged_data$Measles_Count > q_high] <- q_high_new
print(merged_data$Measles_Count)

# Read Texas counties shape file
texas_counties = st_read("~/Research Spring 2025/Research Paper 01 Practice codes/January to March Data/Texas_County_Boundaries_Detailed.shp")
colnames(texas_counties)
colnames(merged_data)

merged_data = merged_data %>%rename(CNTY_FIPS = FIP)
merged_data = merged_data %>%mutate(CNTY_FIPS = as.character(CNTY_FIPS))
texas_measles = texas_counties %>%left_join(merged_data, by = "CNTY_FIPS")
colnames(texas_measles)


#### draw the scatter plot and correlation matrix as in Figure 01 ####

ggpairs(data = texas_measles, columns = c(15, 17, 25, 21))




#### Visualize the choropleth Maps ####
head(texas_measles)

texas_measles = texas_measles %>%rename(countyFIPS = CNTY_FIPS) # CountyMap function need "countyFIPS" command

texas_measles$stateCode = "TX" # CountyMap function need StateCode
texas_measles$stateName = "Texas" # CountyMap function need StateName

class(texas_measles)

# Convert sf object to data frame
texas_measles = as.data.frame(texas_measles)
class(texas_measles)

library(sf)
library(ggplot2)

# Convert the data frame to sf object
texas_measles = st_as_sf(texas_measles)
class(texas_measles)


par(mfrow = c(1,1), mar = c(5.1,5,1.5,2.1))

# Define color scale for Measles Count
measles_color_scale = scale_fill_gradientn(
  colors = c("#FFFFCC", "#FD8D3C", "#FC4E2A", 
             "#E31A1C", "#BD0026", "#800026"), 
  values = scales::rescale(c(0, 1, 5, 10, 30, 50,
                             70, 90, 100, 200, 250, 300)),
  name = "Measles Count"
)

ggplot(data = texas_measles) +
  geom_sf(aes(fill = Measles_Count), color = "black") +
  measles_color_scale +
  theme_minimal()+
  # labs(title = "Measles Count in Texas Counties") +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        axis.text = element_blank(),      # Remove latitude & longitude labels
        axis.ticks = element_blank(),     # Remove axis ticks
        panel.grid = element_blank()     # Remove grid lines)
  )


#Define the color scale to population size
population_color_scale = scale_fill_gradientn(
  colors = c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", 
             "#E31A1C", "#BD0026", "#800026"), 
  values = scales::rescale(c(0, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 1e6, 5e6)),
  name = "Population Size")


# Create the Texas for population size
ggplot(data = texas_measles) +
  geom_sf(aes(fill = Population_Size), color = "black") +
  population_color_scale +
  theme_minimal() +
  #labs(title = "Population size in Texas counties") +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        axis.text = element_blank(),      # Remove latitude & longitude labels
        axis.ticks = element_blank(),     # Remove axis ticks
        panel.grid = element_blank()     # Remove grid lines)
  )


# Define color scale for MMR Vaccination Rate
mmr_color_scale = scale_fill_gradientn(
  colors = c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", 
             "#800026"), 
  values = scales::rescale(c(0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1)),
  name = "MMR Vaccination Rate"
)


# Plot MMR Vaccination Rate Map
ggplot(data = texas_measles) +
  geom_sf(aes(fill = MMR_Avg), color = "black") +
  mmr_color_scale +
  theme_minimal() +
  #labs(title = "MMR Vaccination Rate in U.S. States") +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        axis.text = element_blank(),      # Remove latitude & longitude labels
        axis.ticks = element_blank(),     # Remove axis ticks
        panel.grid = element_blank()     # Remove grid lines)
  )

# Define color scale for Average Temperature
tmp_color_scale = scale_fill_gradientn(
  colors = c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", 
             "#800026"), 
  values = scales::rescale(c(55, 60, 65, 70, 75, 80)),
  name = "Average Temperature"
)


# Plot MMR Vaccination Rate Map
ggplot(data = texas_measles) +
  geom_sf(aes(fill = Average_Temperature), color = "black") +
  tmp_color_scale +
  theme_minimal() +
  #labs(title = "MMR Vaccination Rate in U.S. States") +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        axis.text = element_blank(),      # Remove latitude & longitude labels
        axis.ticks = element_blank(),     # Remove axis ticks
        panel.grid = element_blank()     # Remove grid lines)
  )



# Define color scale for Average Temperature
prcp_color_scale = scale_fill_gradientn(
  colors = c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", 
             "#800026"), 
  values = scales::rescale(c(0, 1, 2, 3, 4, 5, 6, 7)),
  name = "Average Precipitation"
)


# Plot MMR Vaccination Rate Map
ggplot(data = texas_measles) +
  geom_sf(aes(fill = Average_Precipitation), color = "black") +
  prcp_color_scale +
  theme_minimal() +
  #labs(title = "MMR Vaccination Rate in U.S. States") +
  theme(legend.position = "bottom",
        legend.key.size = unit(1, "cm"),
        axis.text = element_blank(),      # Remove latitude & longitude labels
        axis.ticks = element_blank(),     # Remove axis ticks
        panel.grid = element_blank()     # Remove grid lines)
  )
