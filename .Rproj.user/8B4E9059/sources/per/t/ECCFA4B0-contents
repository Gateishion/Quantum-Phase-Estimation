\name{Hvar_Scheme}
\alias{Hvar_Scheme}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Bootstrap estimation for different quantum estimation strategies
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Hvar_Scheme(theta_real, par_space, n_boost, num_prob, n, a, strategy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta_real}{
%%     ~~Describe \code{theta_real} here~~
A numeric value secifying the phase to estimate
}
  \item{par_space}{
%%     ~~Describe \code{par_space} here
A numeric vector specifying the parametric space
}
  \item{n_boost}{
%%     ~~Describe \code{n_boost} here~~
A numeric value specifying the number of repetitions in the bootstrap
}
  \item{num_prob}{
%%     ~~Describe \code{num_prob} here~~
A numeric value specifying the number of probes
}
  \item{n}{
%%     ~~Describe \code{n} here~~
A numeric vector specifying the direction for the probe qubit
}
  \item{a}{
%%     ~~Describe \code{a} here~~
A numeric vector specifying the axis for the rotation
}
  \item{strategy}{
%%     ~~Describe \code{strategy} here~~
A numeric vector specifying the estimation strategy
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
strategy = 1, returns the Holevo variance for Covariant estimation strategy.
strategy = 2, returns the Holevo variance for AQSE strategy.
}
\value{
A number.
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
Marco A. Rodríguez
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

Hvar_Scheme(2,c(0,2*pi),1000,24,c(1/sqrt(3),1/sqrt(3),1/sqrt(3)), c(1/sqrt(2), 0, -1/sqrt(2)),1)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
