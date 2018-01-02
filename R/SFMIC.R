#' SFMIC
#' @export
SFMiC <- function(M = 0, CF1x = NULL, CF2x = NULL, t = NULL, N = NULL, class = NULL, time = 0, dpf=0){

  SFMiC <- SFMiC$new(M = M, CF1x = CF1x, CF2x = CF2x, t = t, N = N, class_id = class, time = time, dpf=dpf)
  desc <- "SFMiC structure, based on TCF (Zhou, 2008)"

  structure(list(description = desc, RObj = SFMiC), class = c("SFMiC"))

}
SFMiC <-
  setRefClass("SFMiC",
              fields = list(
                ## args
                M = "numeric",
                CF1x = "numeric",
                CF2x = "numeric",
                t = "numeric",
                class_id = "character", # cluster prototype
                N = "numeric", # number of associated examples
                ## Store
                c = "numeric",
                meant = "numeric",
                lastTimestamp = "numeric",
                dpf = "numeric"
              ),

              methods = list(
                #x = Data points
                #membership  = Data points's memberhips to the SFMIC
                initialize = function( x = NULL, membership = NULL , class=NULL,time = NULL, dpf=0) {

                  if(!is.null(x)){

                    ## Examples are weighted by their membership value!

                    # Membership can't be NULL
                    if(is.null(membership)) stop("Membership values can not be NULL")

                    M <<- sum(membership)
                    CF1x <<- membership*x
                    CF2x <<- membership*(x^2)
                    t <<- time
                    class_id <<- class
                    N <<- 1
                    c <<- CF1x/M
                    meant <<- t/N
                    lastTimestamp <<- time
                    dpf <<- dpf
                  }
                  else {
                    M <<- numeric()
                    CF1x <<- numeric()
                    CF2x <<- numeric()
                    t <<- numeric()
                    class_id <<- character()
                    c <<- numeric()
                    N <<- 0
                    meant <<- numeric()
                    lastTimestamp <<- time
                    dpf <<- dpf
                  }

                  .self
                }
              )
  )

SFMiC$methods(
  addPoints = function(x, membership = NULL, time,m,...){
    ## Examples are weighted by their membership value!

    # Membership can't be NULL
    if(is.null(membership)) stop("Membership values can not be NULL")

    M <<- M + sum(membership)
    CF1x <<- CF1x + membership*x
    CF2x <<- CF2x + membership*(x^2)
    t <<- t + time
    c <<- CF1x/M
    N <<- N + 1
    meant <<- t/N
    lastTimestamp <<- time

    dpf <<- dpf + (membership^m)*sum((x-c)^2)
    #r <<- sqrt((abs(CF2x)/w) - ((abs(CF1x)/w)^2))
  },

  merge = function(fmic, ...){

    M <<- M + fmic$get_M()
    CF1x <<- CF1x + fmic$get_CF1x()
    CF2x <<- CF2x + fmic$get_CF2x()
    t <<- t + fmic$get_t()
    c <<- CF1x/M
    N <<- N + fmic$get_N()
    dpf <<- dpf + fmic$get_dpf()
    meant <<- t/N
    if (lastTimestamp < fmic$get_lastTimestamp())
      lastTimestamp <<- fmic$get_lastTimestamp()
    #r <<- sqrt((abs(CF2x)/w) - ((abs(CF1x)/w)^2))

  },

  get_M = function(...) { M },
  get_dpf = function(...) { sqrt((1/N)*dpf) },
  get_CF1x = function(...) { CF1x },
  get_CF2x = function(...) { CF2x },
  get_t = function(...) { t },
  get_meant = function(...) {meant},
  get_center = function(...) { c },
  get_N = function(...) {N},
  get_weight = function(...) { M }, # M is the weight of the FMiC
  get_class = function(...) {class_id},
  get_lastTimestamp = function(...) {lastTimestamp},
  get_maxboundary = function(...){
    rmsd <- sqrt(sum(CF2x)/N - sum(CF1x^2)/N^2)
    return(rmsd)
  }
  #  get_radius = function(...) { r }
)

