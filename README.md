# kohonen 

The `kohonen` R package trains supervised and unsupervised self-organizing
maps (SOMs), maps new observations to trained maps, generates predictions,
and visualizes and interrogates map structure. It supports online and batch
training, parallel batch training, multiple data layers with separate distance
measures, and both rectangular and hexagonal grids.

## Installation

Install this fork from GitHub:

``` r
remotes::install_github("AngusWright/kohonen/kohonen")
```

## Using the package

The standard API is centred on `supersom()`, with `som()` and `xyf()` as
single-layer and supervised convenience wrappers. Train a map, map observations
with `map()` or `predict()`, and inspect it with `summary()`, `getCodes()`, and
`plot()`. The built-in plot types include codebook vectors, training changes,
mapping locations, counts, per-unit properties, representation quality, and
neighbour-distance (U-matrix) plots.

This fork also provides a high-level workflow for tabular catalogues:

```r
library(kohonen)

som_map <- kohtrain(
  data = training_data,
  train.expr = c("feature_1", "feature_2", "feature_1 - feature_2"),
  som.dim = c(20, 20)
)

som_map <- kohparse(som_map, new_data)
som_map <- generate.kohgroups(som_map, n.cluster.bins = 20)

group_properties <- generate.kohgroup.property(
  som_map,
  new_data,
  expression = "mean(data$feature_1, na.rm = TRUE)",
  expr.label = "mean_feature_1"
)

plot(som_map, type = "counts")
plot(som_map, type = "property",
     property = group_properties$property$mean_feature_1)
```

`train.expr` accepts column names and R expressions evaluated against `data`.
`kohtrain()` evaluates and robustly scales these features, handles configured
missing values and thresholds, can train from a sparse random or weighted
sample, and stores the scaling parameters on the returned map. Use
`kohparse()` to apply those same parameters to new data and assign each usable
row to its best-matching unit. Set `n.cores` when parallel parsing or training
is appropriate.

## Additions to base (CRAN) version

This repository is a fork of the CRAN package and retains its core SOM API.
Its fork-specific additions focus on catalogue-oriented preprocessing,
clustering, and visualisation:

- `kohtrain()` is a high-level training wrapper that prepares expression-based
  input features, optionally whitens them, supports sparse training samples,
  and records the preprocessing metadata needed for later parsing.
- `kohparse()` applies a trained map's saved whitening parameters to new data,
  maps valid rows to best-matching units, and refreshes existing group
  classifications when applicable.
- `generate.kohgroups()` hierarchically clusters codebook vectors into a
  requested number of groups and stores both cell-level and observation-level
  group assignments on the map. It can also parse and group a supplied new
  dataset in one call.
- `generate.kohgroup.property()` evaluates one or more expressions for every
  generated group and returns a group-level property table together with the
  updated map. This makes group statistics directly usable in downstream
  analysis and property plots.
- Plotting is extended for these group assignments: count plots can use
  clustered classifications and subsets, while property plots accept either
  cell-level or group-level values. The plotting helpers also add logarithmic
  colour scaling, percentile-based colour limits, configurable out-of-range
  colours, and heat-key controls.


## References

-  Wright, A.H., Hildebrandt, H., van den Busch, J.L., and Heymans, C. (2020), Photometric redshift calibration with self-organising maps, Astronomy and Astrophysics, 637, A100. https://doi.org/10.1051/0004-6361/201936782)
- Wehrens, R., & Kruisselbrink, J. (2018). Flexible self-organizing maps in kohonen 3.0. Journal of Statistical Software, 87, 1-18. https://doi.org/10.18637/jss.v087.i07
- Wehrens, R., & Buydens, L. M. C. (2007). Self- and Super-organizing Maps in R: The kohonen Package. Journal of Statistical Software, 21(5), 1–19. https://doi.org/10.18637/jss.v021.i05
