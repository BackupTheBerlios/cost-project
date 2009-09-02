#subset <- getMethod("subset","csData")
predict.true.bayes <- function (csdata, cldata, fit = NULL,  species,
    season.definition = "quarter",
    ageMin = 0, ageMax = 20, wgl.eqn = NULL, age.covariates = NULL,
    weight.covariates = NULL, prediction.cells = NULL, arealist = NULL,
    length.list)
{
    obsdata <- mldata <- NULL
    if (sum(csdata@tr$sampType == "M") > 0)
        mldata <- COSTcore:::subset(csdata, subset = (sampType == "M"))
    if (sum(csdata@tr$sampType == "S") > 0)
        obsdata <- COSTcore:::subset(csdata, subset = (sampType == "S"))
    if (is.null(mldata) & (is.null(obsdata))) {
        print("no data of type M or S")
        return
    }
    models <- make.reg.models(age.covariates, weight.covariates)
    l.int <- log(seq(length.list$minl, length.list$maxl, length.list$int))
    input.data <- read.cost.data(obsdata, mldata, species, wgl.eqn)
    nseas <- 0
    if (season.definition == "quarter")
        nseas <- 4
    if (season.definition == "month")
        nseas <- 12
    if (nseas == 0) {
        print("invalid time strata")
        return
    }
  #  if (nseas == 4) {
  #      input.data$season_obs <- 1 + floor((input.data$season_obs -
 #           1)/3)
 #       input.data$season_mland <- 1 + floor((input.data$season_mland -
 #           1)/3)
 #       input.data$seas <- 1 + floor((input.data$seas - 1)/3)
 #   }
    if (!is.null(input.data$sampsize_disc))
        input.data$sampsize_disc[is.na(input.data$sampsize_disc)] <- 0
    if (!is.null(input.data$haulsize_disc))
        input.data$haulsize_disc[is.na(input.data$haulsize_disc)] <- 0
    input.data$cens.mu = c(1, 30, 3.4)
    input.data$cens.tau = c(1, 1, 1)
    input.data$cens.pri = c(3.4, 0.001, 1e-04, 0.001)
    seas.code <- sim.codecov(input.data$seas)
    gear.code <- sim.codecov(input.data$gear)
    area.code <- sim.codecov(input.data$area, arealist)


    predict <- dbeObject.list <- NULL

        seas.pred <- sim.codecov(prediction.cells$seas, seas.code$covlist)
        gear.pred <- sim.codecov(prediction.cells$gear, gear.code$covlist)
        area.pred <- sim.codecov(prediction.cells$area, area.code$covlist)

        if (prediction.cells$seas == "ALL")
            seas.use <- seas.code$covlist
        else seas.use <- seas.pred$cov
        if (!models$agemodel$Int$seas)
            seas.use <- 1
            
        if (prediction.cells$gear == "ALL")
            gear.use <- gear.code$covlist
        else gear.use <- gear.pred$cov
        if (!models$agemodel$Int$gear)
            gear.use <- 1
            
        if (prediction.cells$area == "ALL")
            area.use <- area.code$covlist
        else area.use <- area.pred$cov
        if (!models$agemodel$Int$area)
            area.use <- 1
            
        year.pred <- 1
        newpred <- expand.grid(year.pred, seas.use, gear.use,
            area.use)
        landings <- NULL
        if (season.definition == "quarter")
            land.seas <- cldata@cl$quarter
        if (season.definition == "month")
            land.seas <- cldata@cl$month
        if (season.definition != "quarter" & season.definition !=
            "month") {
            print("invalid season.definition")
            return
        }
        for (i in 1:length(newpred[, 1])) {
            ind <- cldata@cl$taxon == species
            if (models$agemodel$Int$seas)
                ind <- ind & (land.seas == newpred[i, 2])
            if (models$agemodel$Int$gear)
                ind <- ind & (cldata@cl$foCatEu5 == newpred[i,
                  3])
            if (models$agemodel$Int$area)
                ind <- ind & (cldata@cl$area == newpred[i, 4])
            landings <- c(landings, sum(cldata@cl$landWt[ind]))
        }
        print(landings)
        print(cbind(year.pred, seas.pred$cov, gear.pred$cov,
            area.pred$cov))
        landings <- 1000 * landings
        predict <- predict.fit.COST(fit, fit$COST.list, newpred[,
            1], newpred[, 2], newpred[, 3], newpred[, 4], landings,
            newpred[, 2], t2.year = NULL, t2.seas = NULL, t2.gear = NULL,
            t2.area = NULL, burnin = 0, nMC = 100, l.int = l.int,
            par.haulsize = NULL)
        dbeObject.list <- mbe2dbe.dor(predict, species)

    list(predict = predict, dbeObject.list = dbeObject.list, seas = seas.pred$cov, gear = gear.pred$cov, area = area.pred$cov)
}
