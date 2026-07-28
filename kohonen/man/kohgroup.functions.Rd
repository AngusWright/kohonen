\name{kohparse}
\alias{kohparse}
\alias{generate.kohgroups}
\alias{generate.kohgroup.property}
\alias{kohgroup.loop}
\title{Parse observations and summarise self-organising-map groups}
\description{
  These helpers apply a trained expression-based map to tabular data and
  aggregate the mapped observations. \code{kohparse} maps data using stored
  whitening parameters. \code{generate.kohgroups} hierarchically clusters map
  cells. \code{generate.kohgroup.property} evaluates summary expressions by
  group, and \code{kohgroup.loop} evaluates a single observation-level
  expression within each group.
}
\usage{
kohparse(som, data, train.expr, data.missing = NA,
         data.threshold = c(-Inf, Inf), quiet = FALSE, n.cores = 1,
         max.na.frac)

generate.kohgroups(som, n.cluster.bins = Inf, n.cores = 1, new.data,
                   subset, quiet = FALSE, method = "complete", hclust_in,
                   \dots)

generate.kohgroup.property(som, data, expression, expr.label = NULL,
                           n.cores = 1, n.cluster.bins, subset,
                           quiet = FALSE, returnMatrix = FALSE, \dots)

kohgroup.loop(som, data, expression, expr.label = NULL, n.cores = 1,
              n.cluster.bins, quiet = FALSE, \dots)
}
\arguments{
  \item{som}{A \code{kohonen} object returned by \code{\link{kohtrain}}.}
  \item{data}{Data frame or data table to map or summarise.}
  \item{train.expr}{Optional character vector of feature expressions. By
    default, the feature names stored in \code{som} are used.}
  \item{data.missing, data.threshold}{Missing-value sentinel and accepted
    input range, passed to \code{\link{kohwhiten}}.}
  \item{quiet}{Whether to suppress progress output.}
  \item{n.cores}{Number of cores used while mapping new data.}
  \item{max.na.frac}{Optional replacement for the map's maximum allowed
    missing-data fraction.}
  \item{n.cluster.bins}{Requested number of hierarchical groups. Values at
    least as large as the number of map cells assign each cell to its own
    group. Defaults to the grouping stored on \code{som}.}
  \item{new.data}{Optional data to parse before assigning groups.}
  \item{subset}{Optional logical vector or row indices identifying
    observations to include in the grouping or summary.}
  \item{method}{Linkage method passed to \code{\link{hclust}}.}
  \item{hclust_in}{Optional precomputed \code{hclust} object for the map
    codebook vectors.}
  \item{expression}{Character expression, evaluated with \code{data} bound
    to observations in one group. A vector of expressions is combined into a
    single group summary by \code{generate.kohgroup.property}.}
  \item{expr.label}{Optional names for values returned by \code{expression}.}
  \item{returnMatrix}{If \code{TRUE}, return group properties as a matrix;
    otherwise return a data frame with a \code{group.id} column.}
  \item{\dots}{For \code{generate.kohgroups}, arguments passed to
    \code{kohparse} when \code{new.data} is supplied. For the property and
    loop helpers, arguments passed to \code{generate.kohgroups} if grouping
    must be regenerated.}
}
\value{
  \code{kohparse} returns \code{som} with refreshed \code{data},
  \code{unit.classif}, and \code{distances}. It preserves the training
  assignments in \code{training.classif}.

  \code{generate.kohgroups} returns \code{som} with observation group
  assignments in \code{clust.classif}, cell group assignments in
  \code{cell.clust}, and the requested group count in \code{n.cluster.bins}.

  \code{generate.kohgroup.property} returns a list containing \code{som} and
  a group-level \code{property} data frame or matrix. \code{kohgroup.loop}
  returns a list containing \code{som} and a \code{values} vector with one
  result for each input row.
}
\details{
  Grouping clusters the first codebook layer using
  \code{hclust(dist(som$codes[[1]]))}. Property and loop expressions are
  evaluated dynamically; supply only trusted expressions. Each expression can
  refer to the current group as \code{data}, for example
  \code{"mean(data[[\"alcohol\"]], na.rm = TRUE)"}.
}
\seealso{
  \code{\link{kohtrain}}, \code{\link{kohwhiten}}, \code{\link{plot.kohonen}}.
}
\examples{
data(wines)

wine_features <- as.data.frame(
  wines[, c("alcohol", "magnesium", "flavonoids")]
)
wine_map <- kohtrain(
  wine_features,
  train.expr = colnames(wine_features),
  som.dim = c(3, 3),
  som.iter = 5,
  quiet = TRUE,
  seed = 1
)

wine_map <- kohparse(wine_map, wine_features, quiet = TRUE)
wine_map <- generate.kohgroups(wine_map, n.cluster.bins = 3, quiet = TRUE)

group_means <- generate.kohgroup.property(
  wine_map,
  wine_features,
  expression = "mean(data[[\"alcohol\"]], na.rm = TRUE)",
  expr.label = "mean_alcohol",
  quiet = TRUE
)
group_means$property

alcohol_residuals <- kohgroup.loop(
  group_means$som,
  wine_features,
  expression = "data[[\"alcohol\"]] - mean(data[[\"alcohol\"]])",
  quiet = TRUE
)
head(alcohol_residuals$values)
}
\keyword{classif}
