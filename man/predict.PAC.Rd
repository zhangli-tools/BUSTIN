\name{predict.PAC}
\alias{predict.PAC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
predict.PAC
}
\description{
Predict the phenotype-associated cells.
}
\usage{
predict.PAC(seurat.obj,k.out.list,prob=0.95,minModuleSize=30,weight=NULL,median.sihouette.cutoff=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seurat.obj}{
The Seurat object: single-cell RNA-seq data.
}
  \item{k.out.list}{
List: a list of the output returned by hclust.geneset(). 
}
  \item{prob}{
Numeric, the cutoff for the probability that the cell specifically expresses the gene module, default = 0.95. 
}
  \item{minModuleSize}{
Integer: minimal module size, defualt = 30.
}
  \item{weight}{
Numeric, the gene weight used for weighted analysis. Users can use the test statistics by gene differential expression analysis. 
}
  \item{median.sihouette.cutoff}{
Numeric: the cutoff for the median Sihouette index, default = 0.01.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
function (x) 
{
  }
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
