\name{kohtrain}
\alias{kohtrain}
\alias{kohwhiten}
\title{Train and whiten expression-based self-organising maps}
\description{
  \code{kohtrain} is a convenience wrapper around \code{\link{som}} for
  tabular data. It evaluates selected columns or expressions, optionally
  whitens the resulting features, creates a map grid, and trains the map.
  \code{kohwhiten} performs the feature evaluation and whitening step used by
  \code{kohtrain} and \code{\link{kohparse}}.
}
\usage{
kohtrain(data, train.expr, som.dim = c(10, 10),
         som.topo = "hexagonal", som.toroidal = TRUE, som.iter = 100,
         som.rate = c(0.05, 0.01), n.cores = 1, som.method = "pbatch",
         max.na.frac = 1, train.sparse = FALSE, sparse.frac = 0.1,
         sparse.min.density = 3, sparse.var = NULL, data.missing = NA,
         data.threshold = c(-Inf, Inf), whiten.input = TRUE, quiet = FALSE,
         seed, keep.data = TRUE, \dots)

kohwhiten(data, train.expr, whiten.param, data.missing, data.threshold,
          factor.weight, warn = TRUE, quiet=FALSE)
}
\arguments{
  \item{data}{A data frame or data table containing the training data.}
  \item{train.expr}{Character vector of at least two column names or R
    expressions evaluated against \code{data}. The expressions become the
    feature names in the trained map.}
  \item{som.dim}{Length-two vector giving the horizontal and vertical map
    dimensions.}
  \item{som.topo}{Grid topology passed to \code{\link{somgrid}}:
    \code{"hexagonal"} or \code{"rectangular"}.}
  \item{som.toroidal}{Whether the map grid wraps at its edges.}
  \item{som.iter, som.rate, som.method}{Training iterations, learning-rate
    range, and training mode passed to \code{\link{som}}.}
  \item{n.cores}{Number of cores supplied to \code{\link{som}}.}
  \item{max.na.frac}{Maximum allowed fraction of missing feature values for a
    training row.}
  \item{train.sparse}{If \code{TRUE}, train on a sampled subset of the rows.}
  \item{sparse.frac}{Fraction of rows sampled when \code{train.sparse} is
    \code{TRUE}.}
  \item{sparse.min.density}{Minimum requested sampled observations per map
    unit. The sampling fraction is increased when necessary.}
  \item{sparse.var}{Optional column whose values are used as sampling weights.}
  \item{data.missing}{Value that denotes a missing input measurement. The
    default, \code{NA}, treats \code{NA} values as missing.}
  \item{data.threshold}{Length-two range of accepted input values. Components
    outside this range are treated as missing before feature evaluation.}
  \item{whiten.input}{If \code{TRUE}, centre features by their median and
    divide by their MAD before training.}
  \item{quiet}{Whether to suppress progress output.}
  \item{seed}{Optional random seed, used for reproducible sampling and map
    initialisation.}
  \item{keep.data}{Whether to retain training data in the returned map.}
  \item{\dots}{Additional arguments passed to \code{\link{som}}.}
  \item{whiten.param}{Optional two-row matrix of previously calculated
    medians and MADs. Its columns must correspond to \code{train.expr}.}
  \item{factor.weight}{Optional positive feature weights applied to newly
    calculated MADs. Ignored when \code{whiten.param} is supplied.}
  \item{warn}{Whether to warn if supplied whitening parameters do not have
    matching column names.}
}
\value{
  \code{kohtrain} returns a \code{kohonen} object. In addition to the usual
  components returned by \code{\link{som}}, it records \code{whiten.param},
  \code{bad.index}, \code{training.classif}, and \code{n.cluster.bins}.

  \code{kohwhiten} returns a list with \code{data.white}, the whitened feature
  matrix, and \code{whiten.param}, the two-row matrix of medians and MADs used
  to produce it.
}
\details{
  Missing components and values outside \code{data.threshold} are converted
  to \code{NA} before whitening. Rows with no usable features are removed
  before training and their original row numbers are stored in
  \code{bad.index}. The whitening parameters can subsequently be reused by
  \code{\link{kohparse}} so that new observations are scaled consistently
  with the training data.
}
\seealso{
  \code{\link{som}}, \code{\link{somgrid}}, \code{\link{kohparse}}.
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

whitened <- kohwhiten(
  wine_features,
  train.expr = colnames(wine_features),
  data.missing = NA,
  data.threshold = c(-Inf, Inf)
)
dim(whitened$data.white)
}
\keyword{classif}
