\name{test.contribution}
\alias{test.contribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
test.contribution
}
\description{
Test whether the cell types contribute to the overall gene expression significantly. 
}
\usage{
test.contribution(seurat.obj,geneset,bg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seurat.obj}{
Seurat object: The cell-type-defined single-cell object by Seurat package. The cell types can be accessed by Idents(seurat.obj)
}
  \item{geneset}{
A vector: The gene set of your interest, such as drug resistant genes.
}
  \item{bg}{
A vector: The genes from the whole genome.
}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with a vector of P-values and a vector of genes and corresponding cell type
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
test.contribution(seurat.obj,geneset,bg)
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
