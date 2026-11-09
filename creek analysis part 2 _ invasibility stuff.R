# ============================================================
# Libraries
# ============================================================
library(tidyverse)
library(sf)
library(geosphere)
library(purrr)

# ============================================================
# Step 1: Load Data
# ============================================================
nas <- read_csv("nas.csv")
sites <- read_csv("aremp_sites_attribute_table.csv")
watersheds <- read_csv("aremp_watersheds_attribute_table.csv")
target_ws <- read_csv("ALL_YEARS_watershed_richness.csv")  # 139 watersheds

# Fix typo in watershed names
watersheds <- watersheds %>%
  mutate(CREEK_CODE = if_else(CREEK_CODE == "OREIF", "OREFI", CREEK_CODE))

# Add missing WADNR manually (provided coordinates)
if (!"WADNR" %in% watersheds$CREEK_CODE) {
  watersheds <- watersheds %>%
    add_row(CREEK_CODE = "WADNR",
            LAT_DD = 47.80673077,
            LONG_DD = -124.0082626)
}

# Restrict to target 139 watersheds
target_list <- unique(target_ws[[1]])
watersheds <- watersheds %>% filter(CREEK_CODE %in% target_list)

# ============================================================
# Step 2: Simplify NAS Taxonomic Groups
# ============================================================
macroinv_groups <- c(
  "Mollusks-Bivalves", "Mollusks-Gastropods",
  "Crustaceans-Amphipods", "Crustaceans-Crayfish",
  "Crustaceans-Copepods", "Crustaceans-Isopods",
  "Crustaceans-Crabs", "Crustaceans-Shrimp",
  "Bryozoans", "Coelenterates-Hydrozoans", "Tunicates"
)
amphibian_groups <- c("Amphibians-Frogs", "Amphibians-Salamanders")
fish_groups      <- c("Fishes", "Marine Fishes")
mammal_groups    <- c("Mammals")

nas <- nas %>%
  mutate(Group_Simplified = case_when(
    Group_ %in% fish_groups      ~ "Fishes",
    Group_ %in% amphibian_groups ~ "Amphibians",
    Group_ %in% macroinv_groups  ~ "Macroinvertebrates",
    Group_ %in% mammal_groups    ~ "Mammals",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Group_Simplified)) %>%
  select(Scientific_Name, Group_Simplified, Latitude, Longitude, HUC_8_Number)

# ============================================================
# Step 3: Convert to sf objects
# ============================================================
sites_sf <- st_as_sf(sites, coords = c("Longitude", "Latitude"), crs = 4326)
nas_sf   <- st_as_sf(nas,   coords = c("Longitude", "Latitude"), crs = 4326)

# Load HUC8 shapefile and fix field name
huc8_sf <- st_read("huc8.shape.shp") %>%
  st_transform(4326) %>%
  rename(HUC8 = huc8)

# ============================================================
# Step 4: Spatially join sites to HUC8s
# ============================================================
sites_sf <- st_join(sites_sf, huc8_sf %>% select(HUC8), join = st_within)
sites_sf <- sites_sf %>% filter(!is.na(HUC8))

# ============================================================
# Step 5: Calculate Site-Level Invasibility
# ============================================================
calc_invasibility <- function(group_name) {
  nas_grp <- nas_sf %>% filter(Group_Simplified == group_name)
  
  sites_sf %>%
    rowwise() %>%
    mutate(
      inv_richness = {
        site_huc <- HUC8
        spp_in_huc <- nas_grp %>% filter(HUC_8_Number == site_huc)
        n_distinct(spp_in_huc$Scientific_Name)
      },
      inv_proximity = {
        site_huc <- HUC8
        site_latlon <- st_coordinates(geometry)
        spp_in_huc <- nas_grp %>% filter(HUC_8_Number == site_huc)
        
        if (nrow(spp_in_huc) == 0) {
          0
        } else {
          spp_in_huc %>%
            group_by(Scientific_Name) %>%
            summarise(min_dist = min(distHaversine(site_latlon, st_coordinates(geometry))/1000),
                      .groups = "drop") %>%   # distance in km
            summarise(sum_inv = sum(1 / (min_dist + 1)), .groups = "drop") %>%
            pull(sum_inv)
        }
      },
      inv_index = inv_richness + inv_proximity
    ) %>%
    ungroup() %>%
    st_drop_geometry() %>%
    select(SITE_ID, inv_richness, inv_proximity, inv_index) %>%
    rename_with(~paste0(c("rich_", "prox_", "inv_"), tolower(group_name)),
                c(inv_richness, inv_proximity, inv_index))
}

groups <- unique(nas$Group_Simplified)

site_scores <- reduce(
  map(groups, calc_invasibility),
  full_join,
  by = "SITE_ID"
)

site_scores <- site_scores %>%
  mutate(total_invasibility = rowSums(select(., starts_with("inv_")), na.rm = TRUE)) %>%
  left_join(sites %>% select(SITE_ID, CREEK_CODE), by = "SITE_ID")

# ============================================================
# Step 6: Aggregate to Watersheds → Keep per-taxon rich/prox
# ============================================================
baseline_scores <- site_scores %>%
  group_by(CREEK_CODE) %>%
  summarise(
    fish_rich       = mean(rich_fishes, na.rm = TRUE),
    fish_prox       = mean(prox_fishes, na.rm = TRUE),
    amphib_rich     = mean(rich_amphibians, na.rm = TRUE),
    amphib_prox     = mean(prox_amphibians, na.rm = TRUE),
    mammals_rich    = mean(rich_mammals, na.rm = TRUE),
    mammals_prox    = mean(prox_mammals, na.rm = TRUE),
    macro_rich      = mean(rich_macroinvertebrates, na.rm = TRUE),
    macro_prox      = mean(prox_macroinvertebrates, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  right_join(watersheds %>% select(CREEK_CODE, LAT_DD, LONG_DD), by = "CREEK_CODE") %>%
  mutate(across(everything(), ~replace_na(., 0)))

# ============================================================
# Step 7: Neighbor Smoothing → Create per-taxon spillover
# ============================================================
ws_sf <- st_as_sf(baseline_scores, coords = c("LONG_DD", "LAT_DD"), crs = 4326)
dist_matrix <- st_distance(ws_sf)
closest_neighbors <- apply(dist_matrix, 1, function(x) order(x)[2:3])

neighbor_tbl <- map_dfr(1:nrow(ws_sf), function(i) {
  tibble(
    CREEK_CODE = ws_sf$CREEK_CODE[i],
    N1 = ws_sf$CREEK_CODE[closest_neighbors[1, i]],
    N2 = ws_sf$CREEK_CODE[closest_neighbors[2, i]]
  )
})

neighbor_tbl <- neighbor_tbl %>%
  left_join(baseline_scores, by = c("N1" = "CREEK_CODE")) %>%
  rename_with(~paste0("N1_", .),
              c(fish_rich, fish_prox, amphib_rich, amphib_prox,
                mammals_rich, mammals_prox, macro_rich, macro_prox)) %>%
  left_join(baseline_scores, by = c("N2" = "CREEK_CODE")) %>%
  rename_with(~paste0("N2_", .),
              c(fish_rich, fish_prox, amphib_rich, amphib_prox,
                mammals_rich, mammals_prox, macro_rich, macro_prox))

final_scores <- baseline_scores %>%
  left_join(neighbor_tbl, by = "CREEK_CODE") %>%
  mutate(
    fish_spill    = rowMeans(cbind(N1_fish_rich + N1_fish_prox,
                                   N2_fish_rich + N2_fish_prox), na.rm = TRUE),
    amphib_spill  = rowMeans(cbind(N1_amphib_rich + N1_amphib_prox,
                                   N2_amphib_rich + N2_amphib_prox), na.rm = TRUE),
    mammals_spill = rowMeans(cbind(N1_mammals_rich + N1_mammals_prox,
                                   N2_mammals_rich + N2_mammals_prox), na.rm = TRUE),
    macro_spill   = rowMeans(cbind(N1_macro_rich + N1_macro_prox,
                                   N2_macro_rich + N2_macro_prox), na.rm = TRUE)
  )

# ============================================================
# Step 8: Rescale per-taxon components and build scores
# ============================================================
rescale100 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (rng[1] == rng[2]) rep(50, length(x)) else 1 + (x - rng[1]) * (99 / (rng[2] - rng[1]))
}

# --- No spillover ---
baseline_scores <- baseline_scores %>%
  mutate(
    fish_inv       = rescale100(fish_rich)   + rescale100(fish_prox),
    amphibians_inv = rescale100(amphib_rich) + rescale100(amphib_prox),
    mammals_inv    = rescale100(mammals_rich)+ rescale100(mammals_prox),
    macro_inv      = rescale100(macro_rich)  + rescale100(macro_prox),
    total_inv      = fish_inv + amphibians_inv + mammals_inv + macro_inv
  )

# --- With spillover ---
final_scores <- final_scores %>%
  mutate(
    fish_inv       = rescale100(fish_rich)   + rescale100(fish_prox)   + rescale100(fish_spill),
    amphibians_inv = rescale100(amphib_rich) + rescale100(amphib_prox) + rescale100(amphib_spill),
    mammals_inv    = rescale100(mammals_rich)+ rescale100(mammals_prox)+ rescale100(mammals_spill),
    macro_inv      = rescale100(macro_rich)  + rescale100(macro_prox)  + rescale100(macro_spill),
    total_inv      = fish_inv + amphibians_inv + mammals_inv + macro_inv
  )

# ============================================================
# Step 9: Save Results (only relevant columns)
# ============================================================
baseline_out <- baseline_scores %>%
  rename(Watershed = CREEK_CODE) %>%
  select(Watershed, fish_inv, amphibians_inv, mammals_inv, macro_inv, total_inv)

final_out <- final_scores %>%
  rename(Watershed = CREEK_CODE) %>%
  select(Watershed, fish_inv, amphibians_inv, mammals_inv, macro_inv, total_inv)

write_csv(baseline_out, "C:/Users/hp/Box/eDNA postdoc/Datasets/AREMP/CREEK ANALYSIS/Invasibility_NoSpillover.csv")
write_csv(final_out,    "C:/Users/hp/Box/eDNA postdoc/Datasets/AREMP/CREEK ANALYSIS/Invasibility_WithSpillover.csv")
