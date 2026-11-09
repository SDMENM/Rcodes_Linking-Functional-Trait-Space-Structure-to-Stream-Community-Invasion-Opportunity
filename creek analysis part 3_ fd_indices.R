# ============================================================
# Functional & Taxonomic Diversity Indices at Watershed Level
# ============================================================

library(tidyverse)
library(mFD)
library(funrar)
library(betapart)

# ----------------------------
# Input files
# ----------------------------
traits_all       <- read_csv("traits_all.csv")   # harmonized traits
watershed_taxa   <- read_csv("out_watershed/ALL_YEARS_watershed_taxon_matrix_pa_POOLED.csv")
invasibility_df  <- read_csv("Invasibility_WithSpillover.csv")

# ----------------------------
# Prepare watershed × taxon PA matrix
# ----------------------------
watershed_taxa <- watershed_taxa %>%
  column_to_rownames("Watershed") %>%
  mutate(across(everything(), as.numeric))

# ----------------------------
# Prepare traits
# ----------------------------
traits_comb <- traits_all %>%
  select(TraitKey, Body_Size, Reproductive_Output, Trophic_Role,
         Dispersal_Capacity, Thermal_Preference, Habitat_Preference, Life_Span) %>%
  column_to_rownames("TraitKey") %>%
  mutate(
    Body_Size            = as.numeric(Body_Size),
    Life_Span            = as.numeric(Life_Span),
    Reproductive_Output  = as.ordered(Reproductive_Output),
    Dispersal_Capacity   = as.ordered(Dispersal_Capacity),
    Thermal_Preference   = as.ordered(Thermal_Preference),
    Trophic_Role         = as.factor(Trophic_Role),
    Habitat_Preference   = as.factor(Habitat_Preference)
  )

# ----------------------------
# Align taxa (shared only)
# ----------------------------
shared_taxa <- intersect(colnames(watershed_taxa), rownames(traits_comb))
watershed_taxa <- watershed_taxa[, shared_taxa, drop = FALSE]
traits_comb    <- traits_comb[shared_taxa, , drop = FALSE]

# ----------------------------
# Trait categories
# ----------------------------
trait_types_shared <- tibble(
  trait_name = colnames(traits_comb),
  trait_type = c("Q", "O", "N", "O", "O", "N", "Q")
)

# ----------------------------
# Build functional space (3 axes)
# ----------------------------
dist_comb <- mFD::funct.dist(traits_comb, trait_types_shared, metric = "gower")
fspaces_comb <- mFD::quality.fspaces(dist_comb, maxdim_pcoa = 3,
                                     deviation_weighting = c("absolute", "squared"))

sp_coords <- fspaces_comb$details_fspaces$sp_pc_coord
write_csv(fspaces_comb$details_fspaces$pc_eigenvalues, "FD_eigenvalues.csv")

# ----------------------------
# Richness filter thresholds
# ----------------------------
watershed_richness <- rowSums(watershed_taxa > 0)
watersheds_ok <- list(
  fdis_fori = names(watershed_richness[watershed_richness >= 2]),
  feve      = names(watershed_richness[watershed_richness >= 3]),
  hull      = names(watershed_richness[watershed_richness >= 4])
)

# ----------------------------
# Functional Diversity indices (keep all 139)
# ----------------------------
get_fd <- function(sites, ind_vect) {
  if (length(sites) == 0) return(NULL)
  mFD::alpha.fd.multidim(
    sp_faxes_coord = sp_coords,
    asb_sp_w       = as.matrix(watershed_taxa[sites, , drop = FALSE]),
    ind_vect       = ind_vect,
    scaling        = FALSE
  )$functional_diversity_indices %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Watershed")
}

fd_indices <- bind_rows(
  get_fd(watersheds_ok$fdis_fori, c("fdis", "fori")),
  get_fd(watersheds_ok$feve, "feve"),
  get_fd(watersheds_ok$hull, c("fric", "fdiv"))
) %>%
  group_by(Watershed) %>%
  summarise(across(everything(), ~ first(na.omit(.))), .groups = "drop") %>%
  right_join(
    tibble(Watershed = rownames(watershed_taxa)),
    by = "Watershed"
  ) %>%
  mutate(sp_richn = watershed_richness[Watershed])   # always fill richness

# ----------------------------
# Originality & Redundancy (fix NaN and patch richness=1)
# ----------------------------
asb_sp_w <- as.matrix(watershed_taxa > 0)
storage.mode(asb_sp_w) <- "numeric"
dist_matrix <- as.matrix(dist_comb)

stopifnot(identical(colnames(asb_sp_w), rownames(dist_matrix)))

species_orig <- funrar::distinctiveness(asb_sp_w, dist_matrix, relative = FALSE)

site_fd_df <- tibble(
  Watershed   = rownames(species_orig),
  Originality = rowMeans(species_orig, na.rm = TRUE),
  Redundancy  = rowMeans(1 - species_orig, na.rm = TRUE)
) %>%
  mutate(
    # Replace NaN with NA
    Originality = ifelse(is.nan(Originality), NA, Originality),
    Redundancy  = ifelse(is.nan(Redundancy),  NA, Redundancy),
    # Patch for richness == 1
    Originality = ifelse(watershed_richness[Watershed] == 1 & is.na(Originality), 1, Originality),
    Redundancy  = ifelse(watershed_richness[Watershed] == 1 & is.na(Redundancy), 0, Redundancy)
  )

# ----------------------------
# Compile all results
# ----------------------------
fd_all <- fd_indices %>%
  left_join(site_fd_df, by = "Watershed") %>%
  left_join(invasibility_df, by = "Watershed")

# ----------------------------
# Trait-space Vacancy
# ----------------------------
global_comm <- matrix(1, nrow = 1, ncol = nrow(sp_coords),
                      dimnames = list("Global", rownames(sp_coords)))

global_fric <- mFD::alpha.fd.multidim(
  sp_faxes_coord = sp_coords,
  asb_sp_w       = global_comm,
  ind_vect       = "fric",
  scaling        = FALSE
)$functional_diversity_indices$fric

local_fric <- fd_indices %>% select(Watershed, fric)

trait_vacancy <- local_fric %>%
  mutate(
    GlobalFRic   = global_fric,        # same for all watersheds
    LocalFRic    = fric,               # watershed-specific convex hull volume
    TraitVacancy = 1 - (fric / global_fric)
  ) %>%
  select(Watershed, LocalFRic, GlobalFRic, TraitVacancy)

fd_all <- fd_all %>%
  left_join(trait_vacancy, by = "Watershed")

# ----------------------------
# Taxonomic beta-diversity
# ----------------------------
comm_pa <- as.matrix(watershed_taxa > 0)
storage.mode(comm_pa) <- "numeric"

beta_tax_core <- betapart.core(comm_pa)
beta_tax_pair <- beta.pair(beta_tax_core, index.family = "sorensen")

tax_turnover <- rowMeans(as.matrix(beta_tax_pair$beta.sim), na.rm = TRUE)
tax_nested   <- rowMeans(as.matrix(beta_tax_pair$beta.sne), na.rm = TRUE)
tax_total    <- rowMeans(as.matrix(beta_tax_pair$beta.sor), na.rm = TRUE)

beta_tax_df <- data.frame(
  Watershed      = rownames(comm_pa),
  Tax_Turnover   = tax_turnover,
  Tax_Nestedness = tax_nested,
  Tax_Total      = tax_total
)

# Append to master dataframe
fd_all <- fd_all %>%
  left_join(beta_tax_df, by = "Watershed")

# ----------------------------
# Save final dataframe
# ----------------------------
write_csv(fd_all, "FD_indices_with_invasibility.csv")


