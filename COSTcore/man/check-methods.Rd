\name{check-methods}
\docType{methods}
\alias{checkNms}
\alias{checkNms-methods}
\alias{checkNms,data.frame,character-method}
\alias{checkTys}
\alias{checkTys-methods}
\alias{checkTys,data.frame,list-method}
\alias{checkTRnms}
\alias{checkTRnms-methods}
\alias{checkTRnms,data.frame-method}
\alias{checkHHnms}
\alias{checkHHnms-methods}
\alias{checkHHnms,data.frame-method}
\alias{checkSLnms}
\alias{checkSLnms-methods}
\alias{checkSLnms,data.frame-method}
\alias{checkHLnms}
\alias{checkHLnms-methods}
\alias{checkHLnms,data.frame-method}
\alias{checkCAnms}
\alias{checkCAnms-methods}
\alias{checkCAnms,data.frame-method}
\alias{checkCLnms}
\alias{checkCLnms-methods}
\alias{checkCLnms,data.frame-method}
\alias{checkCEnms}
\alias{checkCEnms-methods}
\alias{checkCEnms,data.frame-method}
\alias{checkTRpk}
\alias{checkTRpk-methods}
\alias{checkTRpk,data.frame-method}
\alias{checkHHpk}
\alias{checkHHpk-methods}
\alias{checkHHpk,data.frame-method}
\alias{checkSLpk}
\alias{checkSLpk-methods}
\alias{checkSLpk,data.frame-method}
\alias{checkHLpk}
\alias{checkHLpk-methods}
\alias{checkHLpk,data.frame-method}
\alias{checkCApk}
\alias{checkCApk-methods}
\alias{checkCApk,data.frame-method}
\alias{checkCEpk}
\alias{checkCEpk-methods}
\alias{checkCEpk,data.frame-method}
\alias{checkCLpk}
\alias{checkCLpk-methods}
\alias{checkCLpk,data.frame-method}
\title{Check methods in COSTcore}
\description{
	Check methods used to validate S4 classes. These are meant for internal usage only, not to be directly called by users.
}
\section{Methods}{
\describe{
	\item{checkNms}{\code{signature(object, names)}: Check names in a data.frame against a character vector.}
	\item{checkTys}{\code{signature(object, tys)}: Check column types in a data.frame against a list of class names.}
	\item{checkTRnms}{\code{signature(object)}: Wrapper for TR checkNms.}
	\item{checkHHnms}{\code{signature(object)}: Wrapper for HH checkNms.}
	\item{checkSLnms}{\code{signature(object)}: Wrapper for SL checkNms.}
	\item{checkHLnms}{\code{signature(object)}: Wrapper for HL checkNms.}
	\item{checkCAnms}{\code{signature(object)}: Wrapper for CA checkNms.}
	\item{checkCLnms}{\code{signature(object)}: Wrapper for CL checkNms.}
	\item{checkCEnms}{\code{signature(object)}: Wrapper for CE checkNms.}
	\item{checkTRpk}{\code{signature(object)}: Check TR primary key.}
	\item{checkHHpk}{\code{signature(object)}: Check HH primary key.}
	\item{checkSLpk}{\code{signature(object)}: Check SL primary key.}
	\item{checkHLpk}{\code{signature(object)}: Check HL primary key.}
	\item{checkCLpk}{\code{signature(object)}: Check CL primary key.}
	\item{checkCEpk}{\code{signature(object)}: Check CE primary key.}
}}
\keyword{methods}

