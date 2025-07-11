% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coat-methods.R
\name{print.coat}
\alias{print.coat}
\alias{coef.coat}
\alias{plot.coat}
\alias{node_baplot}
\title{Methods for Conditional Method Agreement Trees (COAT)}
\usage{
\method{print}{coat}(
  x,
  digits = 2L,
  header = TRUE,
  footer = TRUE,
  title = "Conditional method agreement tree (COAT)",
  ...
)

\method{coef}{coat}(object, node = NULL, drop = TRUE, ...)

\method{plot}{coat}(x, terminal_panel = node_baplot, tnex = 2, drop_terminal = TRUE, ...)

node_baplot(
  obj,
  level = 0.95,
  digits = 2,
  pch = 1,
  cex = 0.5,
  col = 1,
  linecol = 4,
  lty = c(1, 2),
  bg = "white",
  confint = FALSE,
  B = 500,
  cilevel = 0.95,
  cicol = "lightgray",
  xscale = NULL,
  yscale = NULL,
  ylines = 3,
  id = TRUE,
  mainlab = NULL,
  gp = gpar()
)
}
\arguments{
\item{x, object, obj}{a \code{coat} object as returned by \code{\link[coat]{coat}}.}

\item{digits}{numeric. Number of digits used for rounding the displayed coefficients
or limits of agreement.}

\item{header, footer}{logical. Should a header/footer be printed for the tree?}

\item{title}{character with the title for the tree.}

\item{...}{further arguments passed to methods.}

\item{node}{integer. ID of the node for which the Bland-Altman parameters
(coefficients) should be extracted.}

\item{drop}{logical. Should the matrix attribute be dropped if the parameters
from only a single node are extracted?}

\item{terminal_panel}{a panel function or panel-generating function passed to
\code{\link[partykit]{plot.party}}. By default, \code{node_baplot} is used to
generate a suitable panel function for drawing Bland-Altman plots based on the
the provided \code{coat} object. It can be customized using the \code{tp_args} argument
(passed through \code{...}).}

\item{tnex}{numeric specification of the terminal node extension
relative to the inner nodes (default is twice the size).}

\item{drop_terminal}{logical. Should all terminal nodes be "dropped" to
the bottom row?}

\item{level}{numeric level for the limits of agreement.}

\item{pch, cex, col, linecol, lty, bg}{graphical parameters for the scatter plot and limits
of agreement in the Bland-Altman plot (scatter plot character, character extension, plot color,
line color, line types, and background color).}

\item{confint}{logical. Should nonparametric bootstrap percentile confidence intervals be plotted?}

\item{B}{numeric. Number of bootstrap samples to be used if \code{confint = TRUE}.}

\item{cilevel}{numeric. Level of the confidence intervals if \code{confint = TRUE}.}

\item{cicol}{color specification for the confidence intervals if \code{confint = TRUE}.}

\item{xscale, yscale}{numeric specification of scale of x-axis and y-axis, respectively.
By default the range of all scatter plots and limits of agreement across all nodes
are used.}

\item{ylines}{numeric. Number of lines for spaces in y-direction.}

\item{id}{logical. Should node IDs be plotted?}

\item{mainlab}{character or function. An optional title for the plots. Either
a character or a \code{function(id, nobs)}.}

\item{gp}{grid graphical parameters.}
}
\value{
The \code{print()} method returns the printed object invisibly.
The \code{coef()} method returns the vector (for a single node) or matrix (for multiple nodes) of estimated parameters (bias and standard deviation).
The \code{plot()} method returns \code{NULL}.
The \code{node_baplot()} panel-generating function returns a function that can be plugged into the \code{plot()} method.
}
\description{
Extracting information from or visualization of conditional method agreement trees.
Visualizations use trees with Bland-Altman plots in terminal nodes, drawn
via grid graphics.
}
\details{
Various methods are provided for trees fitted by \code{\link[coat]{coat}},
in particular \code{print}, \code{plot} (via \pkg{grid}/\pkg{partykit}) and
\code{coef}. The \code{plot} method draws Bland-Altman plots in the terminal panels by default,
using the function \code{node_baplot}.

In addition to these dedicated \code{coat} methods, further methods are inherited
from \code{\link[partykit]{ctree}} or \code{\link[partykit]{mob}}, respectively,
depending on which \code{type} of \code{coat} was fitted.
}
\examples{
\dontshow{ if(!requireNamespace("MethComp")) {
  if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
    stop("the MethComp package is required for this example but is not installed")
  } else q() }
}
## package and data (reshaped to wide format)
library("coat")
data("scint", package = "MethComp")
scint_wide <- reshape(scint, v.names = "y", timevar = "meth", idvar = "item", direction = "wide")

## conditional method agreement tree
tr <- coat(y.DTPA + y.DMSA ~ age + sex, data = scint_wide)

## illustration of methods (including some customization)

## printing
print(tr)
print(tr, header = FALSE, footer = FALSE)

## extracting Bland-Altman parameters
coef(tr)
coef(tr, node = 1)

## visualization (via grid with node_baplot)
plot(tr)
plot(tr, ip_args = list(id = FALSE),
  tp_args = list(col = "slategray", id = FALSE, digits = 3, pch = 19))

}
