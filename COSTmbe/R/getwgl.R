getwgl<-function(species)
{
       if(species == "Cod") {
               slope <- 2.8571
               a <- log(c(0.0176, 0.0169, 0.0168, 0.0163, 0.0174, 0.0172,
                       0.017, 0.0185, 0.018, 0.0181, 0.0182, 0.0177))
               const <- mean(a)
               seas <- a - const
       }
       if(species == "Haddock"|species=='Melanogrammus aeglefinus') {
               slope <- 2.8268
               a <- log(c(0.0155, 0.0153, 0.0145, 0.0148, 0.015, 0.0151,
                       0.016, 0.0163, 0.0164, 0.017, 0.0164, 0.0159))
               const <- mean(a)
               seas <- a - const
       }
list(const=const, slope=slope,seas=seas)
}
