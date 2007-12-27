\name{accessor-methods}
\docType{methods}
\alias{tr}
\alias{tr-methods}
\alias{tr,csData-method}
\alias{tr,csDataCons-method}
\alias{hh}
\alias{hh-methods}
\alias{hh,csData-method}
\alias{hh,csDataCons-method}
\alias{sl}
\alias{sl-methods}
\alias{sl,csData-method}
\alias{sl,csDataCons-method}
\alias{hl}
\alias{hl-methods}
\alias{hl,csData-method}
\alias{hl,csDataCons-method}
\alias{ca}
\alias{ca-methods}
\alias{ca,csData-method}
\alias{ca,csDataCons-method}
\alias{ce}
\alias{ce-methods}
\alias{ce,ceData-method}
\alias{ce,ceDataCons-method}
\alias{cl}
\alias{cl-methods}
\alias{cl,clData-method}
\alias{cl,clDataCons-method}
\alias{desc}
\alias{desc-methods}
\alias{desc,csData-method}
\alias{desc,clData-method}
\alias{desc,ceData-method}
\alias{desc,csDataCons-method}
\alias{desc,clDataCons-method}
\alias{desc,ceDataCons-method}
\title{Accessor methods in COSTcore}
\description{
	Accessor methods are used to get a slot from an object. It is supposed to implement "encapsulation" in S4, although not completely. Users are strongly advised to use these fucntions on their own code so that future changes on the class definition don't break it.
}
\section{Methods}{
\item{tr}{csData, csDataVal, csDataCons \code{tr} slot accessor.}
\item{hh}{csData, csDataVal, csDataCons \code{hh} slot accessor.}
\item{sl}{csData, csDataVal, csDataCons \code{sl} slot accessor.}
\item{hl}{csData, csDataVal, csDataCons \code{hl} slot accessor.}
\item{ca}{csData, csDataVal, csDataCons \code{ca} slot accessor.}
\item{cl}{clData, clDataVal, clDataCons \code{cl} slot accessor.}
\item{ce}{ceData, ceDataVal, ceDataCons \code{ce} slot accessor.}
\item{desc}{description slot accessor for all classes.}
}
\examples{
	obj <- csData()
	tr(obj)
	hh(obj)
	sl(obj)
	hl(obj)
	ca(obj)
	obj <- ceData()
	ce(obj)
	obj <- clData()
	cl(obj)
}
\keyword{methods}

