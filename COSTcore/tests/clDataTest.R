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
zz <- startTest("clDataTest.txt")
tagTest("clData testing ...")

#! Checks

#! Constructors & accessors
checkRun(clobj <- new("clData"))
checkTrue(is(clobj, "clData"))
checkRun(clobj <- clData())
checkTrue(is(clobj, "clData"))
checkEqual(clobj@cl, cl(clobj))

#! Replacement

#! Selection

finishTest()
