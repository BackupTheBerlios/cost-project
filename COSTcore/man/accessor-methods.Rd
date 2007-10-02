\name{accessor-methods}
\docType{methods}
\alias{tr}
\alias{tr-methods}
\alias{tr,csData-method}
\alias{hh}
\alias{hh-methods}
\alias{hh,csData-method}
\alias{sl}
\alias{sl-methods}
\alias{sl,csData-method}
\alias{hl}
\alias{hl-methods}
\alias{hl,csData-method}
\alias{ca}
\alias{ca-methods}
\alias{ca,csData-method}
\alias{ce}
\alias{ce-methods}
\alias{ce,ceData-method}
\alias{cl}
\alias{cl-methods}
\alias{cl,clData-method}
\alias{desc}
\alias{desc-methods}
\alias{desc,csData-method}
\alias{desc,clData-method}
\alias{desc,ceData-method}
\title{Accessor methods in COSTcore}
\description{
	Accessor methods are used to get a slot from an object. It is supposed to implement "encapsulation" in S4, although not completely. Users are strongly advised to use these fucntions on their own code so that future changes on the class definition don't break it.
}
\section{Methods}{
\item{tr}{csData tr slot accessor.}
\item{hh}{csData hh slot accessor.}
\item{sl}{csData sl slot accessor.}
\item{hl}{csData hl slot accessor.}
\item{ca}{csData ca slot accessor.}
\item{cl}{clData cl slot accessor.}
\item{ce}{ceData ce slot accessor.}
\item{desc}{csData, clData and ceData desc slot accessor.}
}
\keyword{methods}
