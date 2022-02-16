\name{calculate.contribution}
\alias{calculate.contribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate.contribution
}
\description{
calculate the contribution of each cell type to the overall gene expression of the tissue
}
\usage{
calculate.contribution(seurat.obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seurat.obj}{
Seurat object: The cell-type-defined single-cell object by Seurat package. The cell types can be accessed by Idents(seurat.obj)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix with cell proportion of each gene contributed by the cell types. The rows and columns indicate the genes and cell types. 
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
calculate.contribution(seurat.obj)
}

