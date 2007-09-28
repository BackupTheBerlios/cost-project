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
data(soleData)

# start test
setCon()
zz <- startTest("csDataTest.txt")
tagTest("csData testing ...")

#! Checks

checkTrue(checkTRnms(tr))
checkTrue(checkHHnms(hh))
checkTrue(checkSLnms(sl))
checkTrue(checkHLnms(hl))

#! Constructors

#! Replacement

#! Selection

finishTest()
