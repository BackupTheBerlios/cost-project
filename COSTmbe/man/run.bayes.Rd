\name{run.bayes}
\alias{run.bayes}
\title{Run Bayesian analysis to catch at age model}
\description{
Using age, length and weight data, caa fits catch at age model to
data using a Bayesian approach and MCMC simulations.

The model contain three sub-models.

Age model:
Multinomial with probabilities given by a "multilogistic"
relation to covariates which can include intercept (always included),
year-effect, season-effect, gear-effect, area-effect and haul effect
(always included). 

Length given age model:
  \deqn{log(l_{t,f}) =
  \beta_{0,t}+\beta_{1,t}g(a_{t,f})+\varepsilon_{t,f}}
where \eqn{g(a)} is log-linear (scaled to be between zero and one).
Both \eqn{\beta_{0,t}} and \eqn{\beta_{1,t}} are linear in covariates.
\eqn{\beta_{0,t}} can include intercept (always included), year-effect,
season-effect, gear-effect, area-effect and haul effect, while
\eqn{\beta_{1,t}} includes intercept only.

Weight given length model:
  \deqn{log(w_{t,f}) =
  \delta_{0,t}+\delta_{1,t}log(l_{t,f})+\varepsilon_{t,f}}
Both \eqn{\delta_{0,t}} and \eqn{\delta_{1,t}} are linear in covariates.
\eqn{\delta_{0,t}} can include intercept (always included), year-effect,
season-effect, gear-effect, area-effect and haul effect, while
\eqn{\delta_{1,t}} includes intercept only.
}
\usage{run.bayes<-function(COST.data,species,timeStrata="quarter",
                    burnin=2000,thin=1,nmcmc=1000,ageMin=0,ageMax=20,
                    usewglrel=T,ageseas=F,agegear=F,agearea=F,
                    lgaseas=F,lgagear=F,lgaarea=F,
                    wglseas=F,wglgear=F,wglarea=F)
  }
\arguments{
  \item{COST.data}{An object of type ...}
  \item{species}{}
  \item{timeStrata}{Time resolution of sampling. Can be "quarter"
    (default) if sampled by season, or "month" if sampled by months.}
  \item{burnin}{Number of MCMC iterations before storing samples.}
  \item{thin}{Number of iterations between each time samples are
    stored.}
  \item{nmcmc}{Number of samples to be stored. The total number of
    iterations will be \env{burnin} + \env{thin} * env{nmcmc}.}
  \item{ageMin}{Minimum age in the model.}
  \item{ageMax}{Maximum age in the model.}
  \item{usewglrel}{}
  \item{ageseas}{logical; if 'TRUE', the component is included in the
    age model.}
  \item{agegear}{logical; if 'TRUE', the component is included in the
    age model.}
  \item{agearea}{logical; if 'TRUE', the component is included in the
    age model.}
  \item{lgaseas}{logical; if 'TRUE', the component is included in the
    lga model.}
  \item{lgagear}{logical; if 'TRUE', the component is included in the
    lga model.}
  \item{lgaarea}{logical; if 'TRUE', the component is included in the
    lga model.}
  \item{wglseas}{logical; if 'TRUE', the component is included in the
    wgl model.}
  \item{wglgear}{logical; if 'TRUE', the component is included in the
    wgl model.}
  \item{wglarea}{logical; if 'TRUE', the component is included in the
    wgl model.} 
  }
\value{An object \env{fit} from the estimation and an object \env{predict}
  }
