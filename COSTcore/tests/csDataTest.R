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
zz <- startTest("csDataTest.txt")
tagTest("csData testing ...")

data(sole.cs)
data(soleData)
checkTrue(is(sole.cs, "csData"))

#! Checks

checkTrue(checkTRnms(tr))
checkTrue(checkHHnms(hh))
checkTrue(checkSLnms(sl))
checkTrue(checkHLnms(hl))

#! Constructors & accessors
checkRun(csobj <- new("csData"))
checkTrue(is(csobj, "csData"))
checkRun(csobj <- csData())
checkTrue(is(csobj, "csData"))
checkRun(csobj <- csData(tr[,-1], hh[,-1], sl[,-1], hl[,-1]))
checkTrue(is(csobj, "csData"))
checkEqual(csobj@tr, tr(csobj))
checkEqual(csobj@hh, hh(csobj))
checkEqual(csobj@sl, sl(csobj))
checkEqual(csobj@hl, hl(csobj))
checkEqual(csobj@ca, ca(csobj))

#! check object is well formed
# tr
checkRun(o0 <- tr(csobj))
o0 <- c(o0)
names(o0) <- NULL
o <- c(tr)
names(o) <- NULL
checkEqual(o[-1], o0)

# hh
checkRun(o0 <- hh(csobj))
o0 <- c(o0)
names(o0) <- NULL
o <- c(hh)
names(o) <- NULL
checkEqual(o[-1], o0)

# sl
checkRun(o0 <- sl(csobj))
o0 <- c(o0)
names(o0) <- NULL
o <- c(sl)
names(o) <- NULL
checkEqual(o[-1], o0)

# hl
checkRun(o0 <- hl(csobj))
o0 <- c(o0)
names(o0) <- NULL
o <- c(hl)
names(o) <- NULL
checkEqual(o[-1], o0)

#! Replacement

#! Selection

finishTest()
