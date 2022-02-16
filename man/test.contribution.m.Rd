\name{test.contribution.m}
\alias{test.contribution.m}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
test.contribution.m
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
test.contribution.m(seurat.obj,geneset,bg,sample)
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
A vector: Gene set used as the background, usually the genes from the whole genome. 
}
  \item{sample}{
character: The variable represents the sample label in Seurat-meta.data.
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

