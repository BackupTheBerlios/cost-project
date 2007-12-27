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
\item{checkNms}{Check names in a data.frame against a character vector.}
\item{checkTys}{Check column types in a data.frame against a list of class names.}
\item{checkTRnms}{Wrapper for TR checkNms.}
\item{checkHHnms}{Wrapper for HH checkNms.}
\item{checkSLnms}{Wrapper for SL checkNms.}
\item{checkHLnms}{Wrapper for HL checkNms.}
\item{checkCAnms}{Wrapper for CA checkNms.}
\item{checkCLnms}{Wrapper for CL checkNms.}
\item{checkCEnms}{Wrapper for CE checkNms.}
\item{checkTRpk}{Check TR primary key.}
\item{checkHHpk}{Check HH primary key.}
\item{checkSLpk}{Check SL primary key.}
\item{checkHLpk}{Check HL primary key.}
\item{checkCLpk}{Check CL primary key.}
\item{checkCEpk}{Check CE primary key.}
}
\keyword{methods}
