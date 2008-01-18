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
\item{\code{checkNms(object, names)}}{Check names in a data.frame against a character vector.}
\item{\code{checkTys(object, tys)}}{Check column types in a data.frame against a list of class names.}
\item{\code{checkTRnms(object)}}{Wrapper for TR checkNms.}
\item{\code{checkHHnms(object)}}{Wrapper for HH checkNms.}
\item{\code{checkSLnms(object)}}{Wrapper for SL checkNms.}
\item{\code{checkHLnms(object)}}{Wrapper for HL checkNms.}
\item{\code{checkCAnms(object)}}{Wrapper for CA checkNms.}
\item{\code{checkCLnms(object)}}{Wrapper for CL checkNms.}
\item{\code{checkCEnms(object)}}{Wrapper for CE checkNms.}
\item{\code{checkTRpk(object)}}{Check TR primary key.}
\item{\code{checkHHpk(object)}}{Check HH primary key.}
\item{\code{checkSLpk(object)}}{Check SL primary key.}
\item{\code{checkHLpk(object)}}{Check HL primary key.}
\item{\code{checkCLpk(object)}}{Check CL primary key.}
\item{\code{checkCEpk(object)}}{Check CE primary key.}
}
\keyword{methods}
