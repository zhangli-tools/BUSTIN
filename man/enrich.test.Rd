\name{enrich.test}
\alias{enrich.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
enrich.test
}
\description{
Test the overrepresentation of one geneset in the other one. Fisher's exact test was applied. 
}
\usage{
enrich.test(geneset1,geneset2,bg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneset1}{
A vector: Gene set 1.
}
  \item{geneset2}{
A vector: Gene set 2.
}
  \item{bg}{
A vector: Gene set used as the background, usually the genes from the whole genome. 
}

}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The P-value measures the statistical significance. 
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
enrich.test(geneset1, geneset2, bg)
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
