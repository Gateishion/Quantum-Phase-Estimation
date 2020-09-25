\name{arb_sampling}
\alias{arb_sampling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample from an arbitrary continous probability density function
}
\description{
Given a probability density function, the lower bound and upper bound for the domain of the pdf.
The function returns a sample of specified size.
}
\usage{
arb_sampling(f_d, sz, l_s, u_s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f_d}{
A probability density function. f_d should be a function object.
}
\item{sz}{
Size of the sample.
}
\item{l_s}{
Lower bound of the domain.
}
\item{u_s}{
Upper bound of the domain.
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
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
