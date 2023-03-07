\name{get.gene.modules}
\alias{get.gene.modules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
get.gene.modules
}
\description{
Extract the gene modules from the hierarchical clustering tree.
}
\usage{
get.gene.modules=function(k.out,minModuleSize=30,median.sihouette.cutoff=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k.out}{
List: a list returned by hclust.geneset().
}
  \item{minModuleSize}{
Integer: minimal module size, defualt = 30.
}
  \item{median.sihouette.cutoff}{
Numeric: the cutoff for the median Sihouette index, default = 0.01.
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector contains the class labels for the genes.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
get.gene.modules(k.out,minModuleSize=30,median.sihouette.cutoff=0.01)
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
