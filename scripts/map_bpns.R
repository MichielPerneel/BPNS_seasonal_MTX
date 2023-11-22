#remotes::install_github("EMODnet/EMODnetWCS")
#remotes::install_github("EMODnet/EMODnetWFS")

library(EMODnetWFS)
library(EMODnetWCS)
library(ows4R)
library(sf)
library(ggsn)
library(ggplot2)
library(terra)
library(tidyterra)
library(svglite)

# Set stations
stations <- data.frame(
  label = c("ZGO2", "120", "330", "130", "780", "700"),
  longitude = c(2.500000, 2.701167, 2.808333, 2.905000, 3.058000, 3.220000),
  latitude = c(51.333333, 51.185000, 51.433333, 51.270833, 51.471167, 51.376667)
  ) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Set coastal towns
towns <- data.frame(
  label = c("Ostend"),
  latitude = c(51.21404),
  longitude = c(2.92196)
  ) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Define area bounding box
bbox <- st_as_sfc("MULTIPOLYGON (((1.8 51, 1.8 52, 3.8 52, 3.8 51, 1.8 51)))")
st_crs(bbox) <- 4326

# Get coastline
wfs <- WFSClient$new("https://geo.vliz.be/geoserver/MarineRegions/wfs", "2.0.0")

countries <- wfs$getFeatures("MarineRegions:worldcountries_esri_2014",
       cql_filter = "territory IN ('Belgium', 'France', 'Netherlands')") %>%
       st_cast(to = "GEOMETRYCOLLECTION")

# Turn GEOMETRYCOLLECTION into multipolygons
geoms <- lapply(countries$the_geom, `[`)
mp <- lapply(geoms, function(x) sf::st_multipolygon(x = x))
sfc_mp <- sf::st_sfc(mp)
countries$mp <- sfc_mp
countries <- sf::st_set_geometry(countries, sfc_mp)
countries$mp <- NULL
rm(geoms);rm(mp);rm(sfc_mp)

# Crop countries to bounding box
st_crs(countries) <- 4326
countries <- countries %>% st_intersection(bbox)

# Get belgian part of the north sea
be_eez <- wfs$getFeatures("MarineRegions:eez_boundaries",
  cql_filter = "line_id IN (3730, 3729, 128, 127, 126, 125, 124)") %>%
  st_cast(to = "MULTILINESTRING") %>%
  st_cast(to = "LINESTRING")

st_crs(be_eez) <- 4326

# Get bathymetry
wcs <- emdn_init_wcs_client(service = "bathymetry")

cov <- emdn_get_coverage(wcs,
                         coverage_id = "emodnet__mean",
                         bbox = c(xmin = 1.8,
                                  ymin = 51,
                                  xmax = 3.8,
                                  ymax = 52)
)

# Mask bathymetry
cov_masked <- terra::mask(cov, countries, inverse = TRUE)

# Remove bathymetry values above 0
cov_masked[cov_masked > 0] <- NA

# Plot
plot.new()
map <- ggplot() +
  # Add Bathymetry
  geom_spatraster(data = cov_masked) +
  scale_fill_hypso_c(
    palette = "wiki-2.0_bathy",
    name = "Depth (m)",
    breaks = c(0, -20, -40, -60, -80),
    labels = c("0", "20", "40", "60", "80")
    ) +

  # Add countries and boundaries background
  geom_sf(mapping = aes(), data = countries) +

  # Add stations and local towns
  geom_sf(mapping = aes(), data = stations, size = 0.5) +
  geom_sf_text(aes(label = stations$label, geometry = stations$geometry),
              size = 2.5,
              nudge_x = 0,
              nudge_y = 0.05) +
  geom_sf(mapping = aes(), data = towns, size = 0.5) +
  geom_sf_text(aes(label = towns$label, geometry = towns$geometry),
              size = 2.5,
              nudge_x = 0.2,
              nudge_y = 0) +

  # Add north arrow and scalebar
  north(location = "topleft",
        #scale = 0.1,
        symbol = 10, # https://oswaldosantos.github.io/ggsn/
        x.min = 1.8,
        y.min = 51,
        x.max = 3.8,
        y.max = 52) +
  scalebar(location = "bottomright",
           dist = 10,
           height = 0.03,
           dist_unit = "km",
           border.size = 0.03,
           transform = TRUE,
           model = "WGS84",
           x.min = 1.9,
           y.min = 51.05,
           x.max = 3.7,
           y.max = 51.9,
           st.size = 1.8,
           st.dist = 0.03) +

  # Set blank theme
  theme_minimal() +
  theme(text = element_text(size = 8)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5)
        ) +
  labs(title = "Sampling locations in the BPNS")
map
dev.off()

# Save map as 90 mm wide and 74 mm high at 600 dpi as svg
ggsave("figures/environmental/BPNS_map.svg", map,
        width = 3.5, height = 3, unit = "in", dpi = 600)