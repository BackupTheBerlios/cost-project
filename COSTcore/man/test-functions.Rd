\name{Tests}
\alias{checkTrue}
\alias{checkFalse}
\alias{checkRun}
\alias{checkFail}
\alias{checkIdentical}
\alias{checkEqual}
\alias{startTest}
\alias{finishTest}
\alias{tagTest}
\alias{setCon}
\title{Functions to be used in tests.}
\description{
These functions are used to set up unit tests simpler than those implemented by package \pkg{RUnit}. It's an attempt to implement methods that hopefully will allow the implementation of some sQuote{Extreme Programming} best pratices. The main idea with these functions is to set up tests for all funcionalities the system must have and all bugs found along the development process. This pratice will contribute to provide better quality software and will avoid old bugs to surface. For more about sQuote{Extreme Programming} check wikipedia's page \url{http://en.wikipedia.org/wiki/Extreme_Programming}.
}
\usage{
setCon(con="zz")
startTest(file="testReport.txt")
tagTest(tag="My tag !", con=getOption("con"))
checkTrue(x, con = getOption("con"))
checkFalse(x, con = getOption("con"))
checkIdentical(x, y, con = getOption("con"))
checkRun(x, con = getOption("con"))
checkFail(x, con = getOption("con"))
checkEqual(x, y, con = getOption("con"))
finishTest(con=getOption("con"))
}
\arguments{
  \item{x}{R object to be used in test. It maybe an expression for \code{checkTrue}, \code{checkFalse}, \code{checkRun} or \code{checkFail} or an object to compared to \code{y} with \code{checkIdentical} or \code{checkEqual}.}
  \item{y}{R object to be compared to \code{x} in tests \code{checkIdentical} or \code{checkEqual}.}
  \item{con}{Name of communication.}
  \item{file}{Report file name.}
  \item{tag}{A tag.}
}
\details{
\describe{
  \item{setCon}{Set name of connection to file.}
  \item{startTest}{Starts test by creating the report file and opening the connection to it defined by sQuote{con}.}
  \item{tagTest}{Tags the test for future reference. Automatically includes date and time.}
  \item{checkTrue}{Check \sQuote{x} expression is \code{TRUE} and report to file opened by connection defined in \sQuote{con}.}
  \item{checkFalse}{Check \sQuote{x} expression is \code{FALSE} and report to file opened by connection defined in \sQuote{con}.}
  \item{checkIdentical}{Check \sQuote{x} and \sQuote{y} are \code{identical} and report to file opened by connection defined in \sQuote{con}.}
  \item{checkRun}{Check \sQuote{x} expression runs succefully and report to file opened by connection defined in \sQuote{con}.}
  \item{checkFail}{Check \sQuote{x} expression fails to run and report to file opened by connection defined in \sQuote{con}.}
  \item{checkEqual}{Check \sQuote{x} and \sQuote{y} are \code{equal} and report to file opened by connection defined in \sQuote{con}.}
  \item{finishTest}{Closes connection defined by \sQuote{con}.}
}}
\value{A report with the results of the tests.}
\references{These functions were originally developed for FLR (\url{http://flr-project.org}).}
\author{FLCore Team and Ernesto Jardim \email{ernesto@ipimar.pt}}
\note{A simple way to deal with the report is to set up your tests so that the results are allways \code{TRUE}, that way it will be easy to look for results different from \code{TRUE}. That's where your problems are.

The test must be set up in a single file and to be used during \code{R CMD check} it should be placed on directory \sQuote{tests}. The procedure to set up a test is: 
  \enumerate{
     \item{1}{Set up the connection with \code{setCon()}, you may change the name of the connection but in that case you'll have to also change the connection in all functions.}
     \item{2}{Start the test with something like \code{zz <- startTest("clDataTest.txt")}}
     \item{3}{Tag your test for later reference with \code{tagTest("1, 2, testing ...")}}
     \item{4}{Test expressions and comparisons with the \code{check} functions, \emph{e.g.} \code{checkRun(mm <- rnorm(100))}}
     \item{5}{Finish test close the connection to the report file with \code{finishTest()}}
  }
}
\keyword{programming}
\keyword{error}
\keyword{utilities}
\keyword{debugging}

