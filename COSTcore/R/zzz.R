.onLoad <- function(lib,pkg) {
	cat("--------------------------------------------------------------\n")
	pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION", package="COSTcore"), fields=c("Version", "Built")))
	cat("COSTcore - \"Zero\"\n")
	cat(paste("(Version: ", pkg.info["Version"], ". Built on: ", pkg.info["Built"], ")\n", sep=""))
	cat("-----------------------------------------------------------------\n")
}