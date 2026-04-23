# ================================
# 0. Load packages
# ================================
library(tidyverse)
library(vegan)
library(ecotraj)
library(tibble)

# ================================
# 1. Create FAKE community data
# ================================

set.seed(123)

# Define structure
sites <- c("SiteA", "SiteB", "SiteC")
years <- 2018:2022
species <- paste0("Sp", 1:5)

# Create all combinations
meta <- expand.grid(site = sites, year = years) %>%
  arrange(site, year)

# Create fake species counts
fake_data <- meta %>%
  mutate(
    Sp1 = rpois(n(), lambda = 10),
    Sp2 = rpois(n(), lambda = 20),
    Sp3 = rpois(n(), lambda = 5),
    Sp4 = rpois(n(), lambda = 15),
    Sp5 = rpois(n(), lambda = 8)
  )

# Create rownames like "SiteA-2018"
fake_data <- fake_data %>%
  mutate(site_year = paste(site, year, sep = "-"))

# Build matrix (species only)
comm_matrix <- fake_data %>%
  select(site_year, all_of(species)) %>%
  column_to_rownames("site_year")

# ================================
# 2. Extract metadata
# ================================

meta_from_names <- data.frame(
  site_year = rownames(comm_matrix)
) %>%
  separate(site_year, into = c("site", "year"), sep = "-") %>%
  mutate(year = as.numeric(year)) %>%
  arrange(site, year)

# # Reorder matrix to match metadata
# comm_matrix <- comm_matrix[meta_from_names$site_year, ]

# ================================
# 3. Distance matrix
# ================================

dist_bray <- vegdist(comm_matrix, method = "bray")

# ================================
# 4. Define trajectories
# ================================

traj <- defineTrajectories(
  d = dist_bray,
  sites = meta_from_names$site,
  times = meta_from_names$year
)

# ================================
# 5. Plot CTA
# ================================

# Define colors for each site
site_colors <- c(
  "SiteA" = "blue",
  "SiteB" = "red",
  "SiteC" = "darkgreen"
)

png("figures/CTA/cta_example_plot.png", width = 6, height = 6, units = "in", res = 300)

trajectoryPCoA(
  x = traj,
  traj.colors = site_colors,
  survey.labels = TRUE,
  lwd = 2
)

legend("topleft",
       legend = names(site_colors),
       col = site_colors,
       lwd = 2,
       bty = "n")

dev.off()

# ================================
# 6. Compute metrics
# ================================

lengths <- trajectoryLengths(traj)
angles  <- trajectoryAngles(traj)
convergence <- trajectoryConvergence(traj)

# Convert rownames to column for export
lengths_df <- lengths %>% rownames_to_column("site")
angles_df  <- angles %>% rownames_to_column("site")

tau_df <- as.data.frame(as.table(convergence$tau))
colnames(tau_df) <- c("site_from", "site_to", "tau")

p_df <- as.data.frame(as.table(convergence$p.value))
colnames(p_df) <- c("site_from", "site_to", "p_value")

# Save outputs
write_csv(lengths_df, "tables/example_cta_lengths.csv")
write_csv(angles_df, "tables/example_cta_angles.csv")
write_csv(tau_df, "tables/example_cta_convergence_tau.csv")
write_csv(p_df, "tables/example_cta_convergence_pvalues.csv")
