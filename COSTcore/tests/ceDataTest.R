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

#! Checks

#! Constructors & accessors
checkRun(ceobj <- new("ceData"))
checkTrue(is(ceobj, "ceData"))
checkRun(ceobj <- ceData())
checkTrue(is(ceobj, "ceData"))
checkEqual(ceobj@ce, ce(ceobj))

#! Replacement

#! Selection

finishTest()
