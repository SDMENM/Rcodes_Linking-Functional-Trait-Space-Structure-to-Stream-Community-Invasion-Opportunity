# ============================================================
# STEP 1: Setup and file paths
# ============================================================

library(tidyverse)
library(readxl)
library(stringr)

files <- c(
  "AREMP2016-Cutoff468_deduplicated_20250821_Final.xlsx",
  "AREMP2017-Cuttoff520_deduplicated_20250821_Final.xlsx",
  "AREMP2018_Cutoff337_deduplicated_20250821_Final.xlsx",
  "AREMP2019-Cuttoff342_deduplicated_20250904_FinalFixed.xlsx",
  "AREMP2021_Cutoff314_deduplicated_20250821_Final.xlsx"
)

out_root <- "out_watershed"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# STEP 2: Helper functions
# ============================================================

normalize_placeholder_na <- function(x) {
  s <- str_trim(as.character(x))
  s[s %in% c("", " ", NA)] <- NA_character_
  s_low <- str_to_lower(s)
  mask <- s_low %in% c("na","n/a","none","unknown","unidentified","unk","null")
  s[mask] <- NA_character_
  s
}

is_genus_only <- function(species, genus = NULL) {
  x <- normalize_placeholder_na(species)
  g <- normalize_placeholder_na(genus)
  x_low <- str_to_lower(x)
  pat <- paste(c("^sp\\.?$", "^spp\\.?$", "^cf\\.?$", "^aff\\.?$"), collapse = "|")
  is_blank   <- is.na(x_low)
  is_marker  <- !is.na(x_low) & str_detect(x_low, pat)
  is_same_g  <- !is.na(g) & !is.na(x_low) & (x_low == str_to_lower(g))
  is_blank | is_marker | is_same_g
}

is_weird_species <- function(species) {
  x <- normalize_placeholder_na(species)
  ifelse(is.na(x), FALSE,
         str_detect(x, "[0-9]") | str_detect(str_to_lower(x), "\\b[x×]\\b"))
}

exclude_taxa <- c(
  "Homo sapiens",
  "Genus unk_g",
  "unk_g uncultured fungus",
  "Salvelinus environmental sample",
  "Cottus sp. 'Nippenose Valley'"
)

is_excluded_taxon <- function(taxon) {
  taxon %in% exclude_taxa
}

# ============================================================
# STEP 3: Read and clean one year file (example: 2016)
# ============================================================

infile <- files[1]
raw <- read_excel(infile)

dat <- raw
if ("Site#" %in% names(dat) && !"Site" %in% names(dat)) dat <- dat %>% rename(Site = `Site#`)
if ("Watershed" %in% names(dat) && !"Watershed" %in% names(dat)) dat <- dat %>% rename(Watershed = Watershed)

needed <- c("Watershed","Site","Species","Genus","Class","Phylum","Order","Depth")
dat <- dat %>% select(any_of(needed))

# ============================================================
# STEP 4: Standardize and initial filtering
# ============================================================

edna_clean <- dat %>%
  mutate(
    Species = normalize_placeholder_na(Species),
    Genus   = normalize_placeholder_na(Genus),
    Depth   = suppressWarnings(as.numeric(Depth)),
    Watershed = str_to_upper(Watershed)   # FORCE UPPERCASE
  ) %>%
  filter(!is.na(Genus), !is.na(Watershed))

# ============================================================
# STEP 5: Remove controls
# ============================================================

edna_filtered <- edna_clean %>%
  filter(!str_detect(Site, regex("plate", ignore_case = TRUE)),
         str_detect(Site, "\\d|integrated")) %>%
  filter(!is.na(Watershed))

# ============================================================
# STEP 6: Assign taxonomic groups
# ============================================================

edna_filtered <- edna_filtered %>%
  mutate(
    TaxaGroup = case_when(
      Phylum == "Chordata" & Class %in% c(
        "Actinopteri","Actinopterygii","Hyperoartia","Petromyzonti","Petromyzontida",
        "Elasmobranchii","Holocephali","Chondrichthyes","Myxini"
      ) ~ "Fish",
      Phylum == "Chordata" & Class == "Amphibia"  ~ "Amphibian",
      Phylum == "Chordata" & Class == "Mammalia" ~ "Mammal",
      Class %in% c("Insecta","Malacostraca","Branchiopoda",
                   "Bivalvia","Gastropoda",
                   "Arachnida","Maxillopoda","Ostracoda") ~ "Macroinvertebrate",
      Phylum == "Annelida" ~ "Macroinvertebrate",
      Phylum == "Chordata" & is.na(Class) & Order %in% c(
        "Salmoniformes","Cypriniformes","Perciformes","Scorpaeniformes",
        "Acipenseriformes","Siluriformes","Gadiformes","Clupeiformes",
        "Esociformes","Petromyzontiformes"
      ) ~ "Fish",
      TRUE ~ "Other"
    )
  ) %>%
  filter(TaxaGroup %in% c("Fish","Amphibian","Mammal","Macroinvertebrate"))

# ============================================================
# STEP 7: Split into species-level and genus-only
# ============================================================

species_level <- edna_filtered %>%
  filter(!is_genus_only(Species, Genus), !is_weird_species(Species)) %>%
  mutate(Taxon = if_else(str_detect(Species, paste0("^", Genus, "\\b")),
                         Species, paste(Genus, Species)),
         Rank = "species")

genus_level <- edna_filtered %>%
  filter(is_genus_only(Species, Genus)) %>%
  mutate(Taxon = paste("Genus", Genus), Rank = "genus")

edna_kept <- bind_rows(species_level, genus_level)

# ============================================================
# STEP 8: Remove redundant genus rows and excluded taxa
# ============================================================

genus_species_overlap <- species_level %>% distinct(Watershed, Genus)

edna_kept <- edna_kept %>%
  mutate(Genus = if_else(Rank == "species", word(Taxon, 1), str_remove(Taxon, "^Genus "))) %>%
  filter(!(Rank == "genus" & Genus %in% genus_species_overlap$Genus & 
             Watershed %in% genus_species_overlap$Watershed)) %>%
  filter(!is_excluded_taxon(Taxon))

# ============================================================
# STEP 9: Summarize watershed × taxon depth
# ============================================================

species_depth_summary <- edna_kept %>%
  group_by(Watershed, Taxon, Rank, TaxaGroup) %>%
  summarise(TotalDepth = sum(Depth, na.rm = TRUE), .groups = "drop") %>%
  arrange(TaxaGroup, Watershed, Taxon)

# ============================================================
# STEP 10: Build unique taxa table
# ============================================================

unique_taxa <- edna_kept %>%
  distinct(Taxon, Rank, TaxaGroup) %>%
  mutate(
    Genus = if_else(Rank == "species", word(Taxon, 1), str_remove(Taxon, "^Genus ")),
    Species = if_else(Rank == "species", str_remove(Taxon, paste0("^", Genus, "\\s+")), NA_character_),
    TraitKey = if_else(Rank == "species", Taxon, paste("Genus", Genus))
  ) %>%
  arrange(TaxaGroup, Rank, Genus, Species)

# ============================================================
# STEP 11: Write out results for test year
# ============================================================

year_label <- "2016"
year_dir <- file.path(out_root, year_label)
dir.create(year_dir, showWarnings = FALSE, recursive = TRUE)

write_csv(species_depth_summary, file.path(year_dir, paste0(year_label, "_watershed_taxon_depth.csv")))
write_csv(unique_taxa, file.path(year_dir, "unique_taxa.csv"))

# ============================================================
# STEP 12: Function to process all years
# ============================================================

process_one_year <- function(infile) {
  m <- str_match(basename(infile), "([12][0-9]{3})")
  year_label <- if (!is.na(m[,2])) m[,2] else tools::file_path_sans_ext(basename(infile))
  
  raw <- read_excel(infile)
  dat <- raw
  if ("Site#" %in% names(dat) && !"Site" %in% names(dat)) dat <- dat %>% rename(Site = `Site#`)
  if ("Watershed" %in% names(dat) && !"Watershed" %in% names(dat)) dat <- dat %>% rename(Watershed = Watershed)
  
  needed <- c("Watershed","Site","Species","Genus","Class","Phylum","Order","Depth")
  dat <- dat %>% select(any_of(needed))
  
  edna_clean <- dat %>%
    mutate(
      Species   = normalize_placeholder_na(Species),
      Genus     = normalize_placeholder_na(Genus),
      Depth     = suppressWarnings(as.numeric(Depth)),
      Watershed = str_to_upper(Watershed),
      Watershed = if_else(Watershed == "OREIF", "OREFI", Watershed)   
    ) %>%
    filter(!is.na(Genus), !is.na(Watershed))
  
  edna_filtered <- edna_clean %>%
    filter(!str_detect(Site, regex("plate", ignore_case = TRUE)),
           str_detect(Site, "\\d|integrated")) %>%
    filter(!is.na(Watershed)) %>%
    mutate(
      TaxaGroup = case_when(
        Phylum == "Chordata" & Class %in% c("Actinopteri","Actinopterygii","Hyperoartia","Petromyzonti",
                                            "Petromyzontida","Elasmobranchii","Holocephali","Chondrichthyes","Myxini") ~ "Fish",
        Phylum == "Chordata" & Class == "Amphibia"  ~ "Amphibian",
        Phylum == "Chordata" & Class == "Mammalia" ~ "Mammal",
        Class %in% c("Insecta","Malacostraca","Branchiopoda","Bivalvia","Gastropoda",
                     "Arachnida","Maxillopoda","Ostracoda") ~ "Macroinvertebrate",
        Phylum == "Annelida" ~ "Macroinvertebrate",
        Phylum == "Chordata" & is.na(Class) & Order %in% c(
          "Salmoniformes","Cypriniformes","Perciformes","Scorpaeniformes",
          "Acipenseriformes","Siluriformes","Gadiformes","Clupeiformes",
          "Esociformes","Petromyzontiformes"
        ) ~ "Fish",
        TRUE ~ "Other"
      )
    ) %>%
    filter(TaxaGroup %in% c("Fish","Amphibian","Mammal","Macroinvertebrate"))
  
  species_level <- edna_filtered %>%
    filter(!is_genus_only(Species, Genus), !is_weird_species(Species)) %>%
    mutate(Taxon = if_else(str_detect(Species, paste0("^", Genus, "\\b")),
                           Species, paste(Genus, Species)),
           Rank = "species")
  
  genus_level <- edna_filtered %>%
    filter(is_genus_only(Species, Genus)) %>%
    mutate(Taxon = paste("Genus", Genus), Rank = "genus")
  
  edna_kept <- bind_rows(species_level, genus_level)
  
  genus_species_overlap <- species_level %>% distinct(Watershed, Genus)
  
  edna_kept <- edna_kept %>%
    mutate(Genus = if_else(Rank == "species", word(Taxon, 1), str_remove(Taxon, "^Genus "))) %>%
    filter(!(Rank == "genus" & Genus %in% genus_species_overlap$Genus & 
               Watershed %in% genus_species_overlap$Watershed)) %>%
    filter(!is_excluded_taxon(Taxon))
  
  species_depth_summary <- edna_kept %>%
    group_by(Watershed, Taxon, Rank, TaxaGroup) %>%
    summarise(TotalDepth = sum(Depth, na.rm = TRUE), .groups = "drop") %>%
    mutate(Year = year_label) %>%
    arrange(TaxaGroup, Watershed, Taxon)
  
  unique_taxa <- edna_kept %>%
    distinct(Taxon, Rank, TaxaGroup) %>%
    mutate(
      Genus = if_else(Rank == "species", word(Taxon, 1), str_remove(Taxon, "^Genus ")),
      Species = if_else(Rank == "species", str_remove(Taxon, paste0("^", Genus, "\\s+")), NA_character_),
      TraitKey = if_else(Rank == "species", Taxon, paste("Genus", Genus)),
      Year = year_label
    ) %>%
    arrange(TaxaGroup, Rank, Genus, Species)
  
  # --- NEW: write per-year outputs ---
  year_dir <- file.path(out_root, year_label)
  dir.create(year_dir, showWarnings = FALSE, recursive = TRUE)
  
  write_csv(species_depth_summary, file.path(year_dir, paste0(year_label, "_watershed_taxon_depth.csv")))
  write_csv(unique_taxa, file.path(year_dir, "unique_taxa.csv"))
  
  list(year = year_label, watershed_depth = species_depth_summary, unique_taxa = unique_taxa)
}


# ============================================================
# STEP 13: Apply to all years
# ============================================================

results <- lapply(files, process_one_year)

# ============================================================
# STEP 14: Combine results
# ============================================================

all_depth <- bind_rows(lapply(results, `[[`, "watershed_depth")) %>%
  arrange(TaxaGroup, Year, Watershed, Taxon)

all_unique_taxa <- bind_rows(lapply(results, `[[`, "unique_taxa")) %>%
  arrange(TaxaGroup, Year, Genus, Species)

write_csv(all_depth, file.path(out_root, "ALL_YEARS_watershed_taxon_depth.csv"))
write_csv(all_unique_taxa, file.path(out_root, "ALL_YEARS_unique_taxa_species_plus_genus.csv"))

# ============================================================
# STEP 14b: Watershed richness across all years
# ============================================================

watershed_richness <- all_depth %>%
  group_by(Watershed) %>%
  summarise(WatershedRichness = n_distinct(Taxon), .groups = "drop") %>%
  arrange(Watershed)

write_csv(watershed_richness, file.path(out_root, "ALL_YEARS_watershed_richness.csv"))

# ============================================================
# STEP 15: Watershed × Taxon matrices
# ============================================================

watershed_abund <- all_depth %>%
  group_by(Year, Watershed, Taxon) %>%
  summarise(TotalDepth = sum(TotalDepth, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = c(Year, Watershed), names_from = Taxon, values_from = TotalDepth, values_fill = 0)
write_csv(watershed_abund, file.path(out_root, "ALL_YEARS_watershed_taxon_matrix_depth.csv"))

watershed_pa <- watershed_abund %>%
  mutate(across(-c(Year, Watershed), ~ as.integer(. > 0)))
write_csv(watershed_pa, file.path(out_root, "ALL_YEARS_watershed_taxon_matrix_pa.csv"))

# ============================================================
# STEP 15b: Pooled Watershed × Taxon PA matrix (across all years)
# ============================================================

watershed_pooled_pa <- all_depth %>%
  group_by(Watershed, Taxon) %>%
  summarise(Present = as.integer(sum(TotalDepth, na.rm = TRUE) > 0), .groups = "drop") %>%
  pivot_wider(id_cols = Watershed, names_from = Taxon, values_from = Present, values_fill = 0)

write_csv(watershed_pooled_pa, file.path(out_root, "ALL_YEARS_watershed_taxon_matrix_pa_POOLED.csv"))

# ============================================================
# STEP 16: Regional pool
# ============================================================

regional_pool <- all_unique_taxa %>%
  distinct(Taxon, Rank, TaxaGroup, Genus, Species, TraitKey) %>%
  arrange(TaxaGroup, Rank, Genus, Species)
write_csv(regional_pool, file.path(out_root, "ALL_YEARS_unique_taxa_species_plus_genus.csv"))

# ============================================================
# STEP 17: Taxon–trait map
# ============================================================

taxon_map <- regional_pool %>%
  select(Taxon, Rank, Genus, Species, TraitKey, TaxaGroup) %>%
  arrange(TaxaGroup, Rank, Genus, Species)
write_csv(taxon_map, file.path(out_root, "ALL_YEARS_Taxon_TraitKey_Map.csv"))
