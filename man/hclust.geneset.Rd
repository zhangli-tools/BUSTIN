\name{hclust.geneset}
\alias{hclust.geneset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
hclust.geneset
}
\description{
Build a hierarchical clustering tree for resistant genes.
}
\usage{
hclust.geneset(seurat.obj,geneset,maxK,p,percent)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seurat.obj}{
The Seurat object: single-cell RNA-seq data.
}
  \item{geneset}{
Vector: The gene symbols you intend to conduct hierarchical clustering analysis.
}
  \item{maxK}{
Integer: the maximal cluster numbers for hierarchical clustering, greater than 2, default = 10.
}
  \item{p}{
Numeric, the top percentage of expression for discretization, default = 0.05.
}
  \item{percent}{
Numeric, the percentage of non-expressed cells used for excluding genes with lowly expressed cells, default = 0.01.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list contains silhouette index, binary expression matrix, and distance matrix.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
hclust.geneset(seurat.obj,geneset,maxK=10,p=0.05,percent=0.01)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
