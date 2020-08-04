\name{calcBlockW}
\title{Find disjoint spatial blocks of sites}
\alias{calcBlockW}
\description{
Partition the space into disjoint spatial blocks of sites. Call the C++ function \code{calcOrderSite_hpp}. For internal use.
}
\usage{
calcBlockW(W, site_order = NULL, dist.center = NULL, dist.max = Inf, 
         verbose = optionsMRIaggr("verbose"))
}
\arguments{
  \item{W}{the neighbourhood matrix. \emph{dgCMatrix}. REQUIRED.}
  \item{site_order}{a specific order to go all over the sites. \emph{integer vector}.}  
  \item{dist.center}{the distance between each point and a reference point. \emph{numeric vector}.}  
  \item{dist.max}{the neighbourhood range. \emph{numeric vector}.} 
  \item{verbose}{Should the process be verbose over iterations ? \emph{logical}.}
}
\details{
This function requires to have installed the \emph{Matrix} and the \emph{spam} package to work.
If no specific order is set, sites are visitating from the first to the last, according to the neighbourhood matrix.
}
\value{
  An \emph{list} containing :
  \itemize{
  \item\code{[[ls_groups]]} : a \emph{list} containing the index of the sites for each independant group.
  \item\code{[[size_groups]]} : a \emph{vector} containing the size of each independant group.
  \item\code{[[n_groups]]} : an \emph{integer} giving the number of independant groups.
  }
}
\examples{
#### spatial field
\dontrun{
n <- 100
}
\dontshow{
n <- 10
}
coords <- data.frame(which(matrix(0, nrow = n, ncol = n) == 0,arr.ind = TRUE), 1)
optionsMRIaggr(quantiles.legend = FALSE, axes = FALSE, num.main = FALSE, bg = "white")

#### 1- neighbourhood matrix (king) ####
W_king <- calcW(coords, range = 1.001, row.norm = TRUE)$W

#### find independant groups
Block_king <- calcBlockW(W_king)

## check groups
# diagonal : percent of neighborhing sites whithin group
# extra-diagonal : percent of neighborhing sites between groups
sapply(1:Block_king$n_groups, function(x){
  sapply(1:Block_king$n_groups, function(y){
    sum(spam::rowSums(W_king[Block_king$ls_groups[[x]], Block_king$ls_groups[[y]]] > 0) > 0)
  }) / length(Block_king$ls_groups[[x]])
}
)

## diplay sparse matrix
spam::image(W_king)
spam::image(W_king[unlist(Block_king$ls_groups), unlist(Block_king$ls_groups)])

## display site blocks
col_sites <- unlist(lapply(1:Block_king$n_groups, function(x){
	rep(rainbow(Block_king$n_groups)[x], Block_king$size_groups[x])
}))

multiplot(coords[unlist(Block_king$ls_groups),],
          xlim = c(0,30),ylim = c(0,30),
          col = col_sites, legend = FALSE)


#### 2- neighbourhood matrix (Queen) ####
W_queen <- calcW(coords, range = sqrt(2) + 0.001, row.norm = TRUE)$W

#### find independant groups
Block_queen <- calcBlockW(W_queen)

## check groups
# diagonal : percent of neighborhing sites whithin group
# extra-diagonal : percent of neighborhing sites between groups
sapply(1:Block_queen$n_groups, function(x){
  sapply(1:Block_queen$n_groups, function(y){
    sum(spam::rowSums(W_queen[Block_queen$ls_groups[[x]], Block_queen$ls_groups[[y]]] > 0) > 0)
  }) / length(Block_queen$ls_groups[[x]])
}
)

## diplay sparse matrix
spam::image(W_queen)
spam::image(W_queen[unlist(Block_queen$ls_groups), unlist(Block_queen$ls_groups)])

## display site blocks
col_sites <- unlist(lapply(1:Block_queen$n_groups, function(x){
	rep(rainbow(Block_queen$n_groups)[x], Block_queen$size_groups[x])
}))

multiplot(coords[unlist(Block_queen$ls_groups),],
          xlim = c(0,30), ylim = c(0,30),
          col = col_sites, legend = FALSE)

#### 3- neighbourhood matrix (Regional) ####
W_Regional <- calcW(coords, range = 3, row.norm = TRUE)$W

#### find independant groups
system.time(
  Block_Regional <- calcBlockW(W_Regional)
)

system.time(
Block_Regional_test1 <- calcBlockW(W_Regional, 
     dist.center = sqrt(spam::rowSums(sweep(coords, MARGIN = 2, 
	                                  STATS = apply(coords, 2, median), FUN = "-")^2))
     )
)
system.time(
  Block_Regional_test2 <- calcBlockW(W_Regional, 
     dist.center = sqrt(spam::rowSums(sweep(coords, MARGIN = 2,
                            	      STATS = apply(coords, 2, median), FUN = "-")^2)),
     dist.max = 3
  )
)
# all(unlist(Block_Regional_test1$ls_groups) == unlist(Block_Regional_test2$ls_groups))


## check groups
# diagonal : percent of neighborhing sites whithin group
# extra-diagonal : percent of neighborhing sites between groups
sapply(1:Block_Regional$n_groups,function(x){
   sapply(1:Block_Regional$n_groups,function(y){
    if(length(Block_Regional$ls_groups[[x]]) > 1){
      sum(spam::rowSums(as.matrix(W_Regional[Block_Regional$ls_groups[[x]],
	                                   Block_Regional$ls_groups[[y]]]) > 0) > 0)
    }else{
      sum(W_Regional[Block_Regional$ls_groups[[x]],
	                 Block_Regional$ls_groups[[y]]] > 0) > 0
    }
  }) / length(Block_Regional$ls_groups[[x]])
}
)
# clustering could be improved

## diplay sparse matrix
spam::image(W_Regional)
spam::image(W_Regional[unlist(Block_Regional$ls_groups), unlist(Block_Regional$ls_groups)])

## display site blocks
col_sites <- unlist(lapply(1:Block_Regional$n_groups, function(x){
rep(rainbow(Block_Regional$n_groups)[x], Block_Regional$size_groups[x])
}))

multiplot(coords[unlist(Block_Regional$ls_groups),],
          xlim = c(0,30), ylim = c(0,30),
          col = col_sites, legend = FALSE)
}
\concept{calc.}
\keyword{functions}

