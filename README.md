# Community Trajectory Analysis (CTA) Workshop

This repository contains materials, code, and examples for conducting **Community Trajectory Analysis (CTA)** in R using the `ecotraj` framework. The workshop focuses on analyzing multivariate community change through time using ordination-based trajectory methods.

---

## 📌 Overview

Community Trajectory Analysis (CTA) is a framework for quantifying and comparing how ecological communities change through time in multivariate space. It is commonly used to assess:

- Directionality of community change
- Magnitude of change (trajectory length)
- Convergence or divergence among sites
- Responses to disturbance or environmental gradients

This workshop walks through:
1. Building community datasets
2. Running ordinations (e.g., PCoA)
3. Generating trajectories
4. Extracting CTA metrics
5. Visualizing results

---

## 📦 Installation

Install required packages:

```r
install.packages(c("tidyverse", "vegan"))

# ecotraj (GitHub version recommended)
install.packages("devtools")
devtools::install_github("EMF-creaf/ecotraj")
