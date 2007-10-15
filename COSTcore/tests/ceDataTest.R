#=====================================================================
#
# Date: 13/04/2007
# Version: 0.1-0
# Authors: Ernesto Jardim
#
# Short description: tests for FLlst
#
# ToDo:
#
# References (bibtex):
#
#!Notes:
#
#=====================================================================

library(COSTcore)

# start test
setCon()
zz <- startTest("ceDataTest.txt")
tagTest("ceData testing ...")

data(sole)
data(soleData)
checkTrue(is(sole.ce, "ceData"))

#! Checks
checkTrue(checkCEnms(ce))

#! Constructors
checkRun(ceobj <- new("ceData"))
checkTrue(is(ceobj, "ceData"))
checkRun(ceobj <- ceData())
checkTrue(is(ceobj, "ceData"))
checkEqual(ceobj@ce, ce(ceobj))
checkRun(ceobj <- ceData(ce[,-1]))
checkTrue(is(ceobj, "ceData"))

# Accessors
checkEqual(ceobj@ce, ce(ceobj))

#! check object is well formed
checkRun(o0 <- ce(ceobj))
o0 <- c(o0)
names(o0) <- NULL
o <- c(ce)
names(o) <- NULL
checkEqual(o[-1], o0)

#! Replacement
checkFail(sole.ce[ce(sole.ce)$quarter==1,"quarter"]<-5)
checkRun(sole.ce[,"daysAtSea"]<-10)
checkTrue(is(sole.ce,"ceData"))
checkFail(sole.ce[,"quarter"]<-10)

#! Selection
checkRun(sole.ce[ce(sole.ce)$quarter==1,])
checkTrue(is(sole.ce[ce(sole.ce)$quarter==1,], "ceData"))
checkRun(sole.ce[ce(sole.ce)$quarter==1,2])
checkTrue(is(sole.ce[ce(sole.ce)$quarter==1,2], "factor"))
checkRun(sole.ce[,2])
checkTrue(is(sole.ce[,2], "factor"))
checkRun(sole.ce[ce(sole.ce)$quarter==1,"quarter"])
checkTrue(is(sole.ce[ce(sole.ce)$quarter==1,"quarter"],"factor"))
checkRun(sole.ce[,"quarter"])
checkTrue(is(sole.ce[,"quarter"],"factor"))
checkRun(sole.ce[ce(sole.ce)$quarter==1,c(1:4)])
checkTrue(is(sole.ce[ce(sole.ce)$quarter==1,c(1:4)],"data.frame"))
checkRun(sole.ce[ce(sole.ce)$quarter==1,c("year","quarter")])
checkTrue(is(sole.ce[ce(sole.ce)$quarter==1,c("year","quarter")],"data.frame"))
checkRun(subset(sole.ce, quarter==1))
checkTrue(is(subset(sole.ce, quarter==1), "ceData"))
checkRun(subset(sole.ce, quarter==1 & area=="27.7.e"))
checkTrue(is(subset(sole.ce, quarter==1 & area=="27.7.e"), "ceData"))

#! utils methods
checkRun(head(sole.ce))
checkRun(tail(sole.ce))
checkRun(summary(sole.ce))
checkRun(dim(sole.ce))
checkTrue(is.ceData(sole.ce))

finishTest()
