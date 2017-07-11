#######################################################################################################
# Code provided by                                                                                    #
#                                                                                                     #
#    John Fox, Professor                                                                              #
#    Department of Sociology                                                                          #
#    McMaster University                                                                              #
#    Hamilton, Ontario, Canada                                                                        #
#    web: socserv.mcmaster.ca/jfox                                                                    #
#######################################################################################################

cochrane.orcutt <- function(mod, ...){
     UseMethod("cochrane.orcott")
     }

cochrane.orcutt.lm <- function(mod){
     X <- model.matrix(mod)
     y <- model.response(model.frame(mod))
     e <- residuals(mod)
     n <- length(e)
     names <- colnames(X)
     rho <- sum(e[1:(n-1)]*e[2:n])/sum(e^2)
     y <- y[2:n] - rho * y[1:(n-1)]
     X <- X[2:n,] - rho * X[1:(n-1),]
     mod <- lm(y ~ X - 1)
     result <- list()
     result$coefficients <- coef(mod)
     names(result$coefficients) <- names
     summary <- summary(mod, corr = F)
     result$cov <- (summary$sigma^2) * summary$cov.unscaled
     dimnames(result$cov) <- list(names, names)
     result$sigma <- summary$sigma
     result$rho <- rho
     class(result) <- 'cochrane.orcutt'
     result
     }

prais.winsten <- function(mod, ...){
     UseMethod("prais.winsten")
     }

prais.winsten.lm <- function(mod){
     X <- model.matrix(mod)
     y <- model.response(model.frame(mod))
     e <- residuals(mod)
     n <- length(e)
     names <- colnames(X)
     rho <- sum(e[1:(n-1)]*e[2:n])/sum(e^2)
     y <- c(y[1] * (1 - rho^2)^0.5, y[2:n] - rho * y[1:(n-1)])
     X <- rbind(X[1,] * (1 - rho^2)^0.5, X[2:n,] - rho * X[1:(n-1),])
     mod <- lm(y ~ X - 1)
     result <- list()
     result$coefficients <- coef(mod)
     names(result$coefficients) <- names
     summary <- summary(mod, corr = F)
     result$cov <- (summary$sigma^2) * summary$cov.unscaled
     dimnames(result$cov) <- list(names, names)
     result$sigma <- summary$sigma
     result$rho <- rho
     class(result) <- 'prais.winsten'
     result
     }
