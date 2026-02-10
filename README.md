## Project Title
Optimizing Water Access in the DRC (Two-Stage MIP)

**TL;DR:** Built a two-stage mixed-integer optimization pipeline to plan water well repairs and new construction in the Democratic Republic of Congo, maximizing population coverage under budget constraints and then minimizing average walking distance using real population and infrastructure data.

## Goal
Design a data-driven decision framework to improve access to basic water services in the Democratic Republic of Congo by optimally allocating a limited budget across:
- repairing non-functional wells, and
- constructing new wells,
while prioritizing both **population coverage** and **accessibility**.

## Approach
- Formulated a **two-stage mixed-integer program (MIP)**:
  - **Stage 1:** maximize total population served under a fixed budget by selecting wells to repair or build
  - **Stage 2:** minimize population-weighted average walking distance while preserving the coverage achieved in Stage 1
- Explicitly modeled:
  - repair vs. new construction costs
  - budget constraints
  - service radius constraints
  - assignment of demand cells to wells

### Data
- **WorldPop DRC population raster**  
  Used to represent spatial population demand.
- **Water Point Data Exchange (WPDx)**  
  Used to identify functional and non-functional water points.

### Demand Modeling
- Constructed two demand representations:
  - a uniform square grid, and
  - a population-weighted **K-means demand grid**
- Selected the K-means grid to better reflect population density while keeping the optimization tractable.

### Candidate Well Generation
- Generated **3,000 candidate well locations** by:
  - ranking unserved demand cells by population,
  - placing candidates at demand centroids,
  - enforcing a minimum spacing constraint to avoid clustering.

### Evaluation / Analysis
- Budget sensitivity analysis across multiple budget levels.
- Comparison of Stage-1 vs. Stage-2 decisions.
- Examined diminishing returns in coverage and walking distance improvements.

## Results
- Early budget increases yield the largest gains in both population coverage and walking distance reduction.
- Later spending exhibits diminishing returns as interventions expand into sparsely populated regions.
- The second optimization stage slightly adjusts site selection to improve accessibility while preserving maximum coverage.
- Population-adaptive demand grids improve realism compared to a uniform grid.

## Repo structure
- `src/` — optimization models and preprocessing scripts  
- `notebooks/` — exploratory analysis and visualization  
- `data/` — processed inputs (raw data excluded)  
- `report/` — final LaTeX / Overleaf report  

## How to run
1. Download raw data:
   - WorldPop population raster (DRC)
   - Water Point Data Exchange (WPDx)
2. Preprocess demand and candidate wells:
   ```bash
   julia src/build_demand_grid.jl
   julia src/build_candidate_wells.jl
