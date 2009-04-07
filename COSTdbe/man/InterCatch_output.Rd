\name{InterCatch_output}
\alias{makeICdf}
\alias{makeICfile}

\title{ToDo}
\description{
ToDo}

\usage{
makeICdf(dbeOutput, filename=NA, output.df = TRUE, append = FALSE,CANUMtype = "length", Country, Species = NA, ReportingCategory,
  Usage = NA, DataToFrom = NA, SamplesOrigin,InfoFleet = NA, InfoStockCoordinator = NA, InfoGeneral = NA, Sex = "N", PlusGroup = -9, unitMeanWeight = "g",
  unitCANUM = "n", Maturity = NA)
  
makeICfile(ICdfs, filename, append=FALSE)
}

\arguments{
  \item{dbeOutput}{ToDo}
  \item{ICdfs}{ToDo}
  \item{filename}{ToDo}
  \item{output.df}{ToDo}
  \item{append}{ToDo}
  \item{CANUMtype}{ToDo}
  \item{Country}{ToDo}
  \item{Species}{ToDo}
  \item{ReportingCategory}{ToDo}
  \item{Usage}{ToDo}
  \item{DataToFrom}{ToDo}
  \item{SamplesOrigin}{ToDo}
  \item{InfoFleet}{ToDo}
  \item{InfoStockCoordinator}{ToDo}
  \item{InfoGeneral}{ToDo}
  \item{Sex}{ToDo}
  \item{PlusGroup}{ToDo}
  \item{unitMeanWeight}{ToDo}
  \item{unitCANUM}{ToDo}
  \item{Maturity}{ToDo}
  
}
%\details{ 
  
%}
\value{

ToDo
}

\author{David Maxwell}


\seealso{\code{\link{dbeOutput}} 
}
\examples{

# load example data set
data("LEMexample")
# Setup object for output
LEM.dbeOut.an = dbeObject(desc="Example",species="Microstomus kitt",param="landings",strataDesc=LEM.strat, catchCat="LAN",methodDesc="analytical")

## Run vesselRaise function
LEM.dbeOut.an = vesselRaise.an (csObject = LEM.CScon, clObject = LEM.CLcon, dbeOutp = LEM.dbeOut.an)

# In two steps
# Create list with HI, SI, SD dataframes
LEM.ICdf <- makeICdf (LEM.dbeOut.an, Country = "UKE", Species = "LEM", ReportingCategory = "R", SamplesOrigin = "M")
# edit dataframes in LEM.ICdf if required then write InterCatch file
# Not run
makeICfile (LEM.ICdf, filename ="IC_LEM.csv")

# In one step, straight to file. output.df=F stops dataframes being saved within R
# Not run
makeICdf (LEM.dbeOut.an, filename ="IC_LEM.csv", output.df = FALSE, Country = "UKE", Species = "LEM", ReportingCategory = "R", SamplesOrigin = "M")

}

\keyword{attribute}

