

setClass("PerformStats",
	representation(
    desc="character",                        # descriptor
    species="character",                     # recall of SL$spp (+ SL$sex)
	catchCat="character",                    # recall of the catch category (discards/landings)
	param="character",                       # recall of the parameter estimated (N, W, maturity, sex-ratio,...)
	strataDesc="strIni",                     # time, space and technical stratification considered
	methodDesc="character",                  # recall of the method (analytical, bootstrap, bayesian)
	nSamples = "numeric",                    # number of samples
    ageTrue = "data.frame",                  # True age structure (param-at-length)
    ageEst = "data.frame",                   # mean estimates of the age structure (param-at-length)
	ageAcc = "list",                         # a list with Accuracy statistics at age.
	ageBias = "list",                        # a list with Bias statistics at age.
	agePrec = "list",                        # a list with Precision statistics at age.
    lengthTrue = "data.frame",               # True age structure (param-at-length)
    lengthEst = "data.frame",                # mean estimates of the length structure (param-at-length)
	lengthAcc = "list",                      # a list with Accuracy statistics at length.
	lengthBias = "list",                     # a list with Bias statistics at length.
	lengthPrec = "list"                     # a list with Precision statistics at length.
	),
	prototype(
        desc        = "PerformStats Object",
		species     = as.character(NA),
		catchCat    = as.character(NA),
		strataDesc  = strIni(),
		methodDesc  = as.character(NA),
		nSamples    = as.integer(NA),
# AGE -----------------------------------------------------------------------
        ageTrue     = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA)),
        ageEst      = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA)),
        ageAcc    = list(mse  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         rmse = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         mae = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         smse = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         srmse = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         smae = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA))),
        ageBias   = list(me  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         sme = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         par = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         poe = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA))),
        agePrec   = list(var  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         cv = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA))),
# LENGTH -----------------------------------------------------------------------
        lengthTrue     = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA)),
        lengthEst      = data.frame(time     = as.character(NA),
                                space     = as.character(NA),
                                technical = as.character(NA),
                                value     = as.numeric(NA)),
        lengthAcc    = list(mse  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         rmse = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         mae = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         smse = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         srmse = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         smae = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA))),
        lengthBias   = list(me  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         sme = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         par = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         poe = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA))),
        lengthPrec   = list(var  = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)),
                         cv = data.frame(time      = as.character(NA),
                                           space     = as.character(NA),
                                           technical = as.character(NA),
                                           value     = as.numeric(NA)))
    )
)





