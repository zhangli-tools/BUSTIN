\name{run.limma}
\alias{run.limma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
run.limma
}
\description{
Identify the upregulated genes in resistant samples using bulk data by limma method
}
\usage{
run.limma(bulk.data, pdata, resistant=T, padj, log2fc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bulk.data}{
Matrix: The log2 normalized bulk gene expression data. The rows and columns indicate the genes and samples, respectively.
}
  \item{pdata}{
Vector: The variable indicates resistant or sensitive to drug treatment. 1 and 0 represent resistant and sensitive, respectively. 
}
  \item{resistant}{
Logic: TRUE and FALSE indicate return the resistant and sensitivie genes, respectively. 
}
  \item{padj}{
Numeric: The cutoff for the adjusted P-value, default=0.05. 
}
  \item{log2fc}{
Numeric: The cutoff for the log2 fold change, default=0.5. 
}
}


\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector of genes resistant/sensitive to drug treatment.
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
