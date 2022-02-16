\name{identify.subpop}
\alias{identify.subpop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
identify.subpop
}
\description{
Identify subpopulations from a specified cell type
}
\usage{
identify.subpop(seurat.obj,geneset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seurat.obj}{
Seurat object: it has only one cell type.
}
  \item{geneset}{
A vector, The gene set of your interest, such as drug resistant genes.
}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Seurat object: The subpopulation associated with the geneset. 
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
identify.subpop(seurat.obj,geneset)

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
