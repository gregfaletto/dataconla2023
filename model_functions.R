## @knitr models

sim_data <- function(n, p, K, intercepts, beta_mat, dev_num=NA, dev_size=NA){

    stopifnot(K >= 3)
    stopifnot(length(intercepts) == K - 1)
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))

    stopifnot(is.numeric(beta_mat) | is.integer(beta_mat))
    stopifnot(nrow(beta_mat) == p)
    stopifnot(ncol(beta_mat) == K - 1)

    # If we don't observe at least one observation from each class, toss
    # out simulation and try again
    all_classes <- FALSE
    counter <- 0
    while(!all_classes){

        if(counter >= 10){
            return(NA)
        }
        
        # Generate uniformly distributed X
        if(p > 1){
            X <- matrix(runif(n*p, min=-1, max=1), nrow=n, ncol=p)
        } else{
            X <- runif(n, min=-1, max=1)
        }

        # Probabilities of each class
        probs <- matrix(as.numeric(NA), nrow=n, ncol=K)
        if(p > 1){
            probs[, 1] <- plogis(intercepts[1] + X %*% beta_mat[, 1])
        } else{
            probs[, 1] <- plogis(intercepts[1] + beta_mat[, 1]*X)
        }

        if(p > 1){
            for(k in 2:(K-1)){
                if(k > 2){
                    probs[, k] <- plogis(intercepts[k] + X %*% beta_mat[, k]) -
                        rowSums(probs[, 1:(k - 1)])
                } else{
                    probs[, k] <- plogis(intercepts[k] + X %*% beta_mat[, k]) -
                        probs[, 1]
                }
                # stopifnot(all(probs[, k] > 0))
            }
        } else{
            for(k in 2:(K-1)){
                if(k > 2){
                    probs[, k] <- plogis(intercepts[k] + beta_mat[, k]*X) -
                        rowSums(probs[, 1:(k - 1)])
                } else{
                    probs[, k] <- plogis(intercepts[k] + beta_mat[, k]*X) -
                        probs[, 1]
                }
                # stopifnot(all(probs[, k] > 0))
            } 
        }

        stopifnot(all(rowSums(probs[, 1:(K-1)]) <= 1))

        probs[, K] <- 1 - rowSums(probs[, 1:(K-1)])

        stopifnot(all(!is.na(probs)))
        stopifnot(all(abs(rowSums(probs) - 1) < 10^(-7)))

        if(any(probs < 0)){
            return(NA)
        }
        if(any(probs > 1)){
            return(NA)
        }

        stopifnot(all(probs >= 0))
        stopifnot(all(probs <= 1))

        # Generate y
        y <- rep(as.numeric(NA), n)
        for(j in 1:n){
            y[j] <- which(rmultinom(n=1, size=1, prob=probs[j,])[, 1] != 0)
        }
        stopifnot(all(!is.na(y)))
        stopifnot(all(y %in% 1:K))

        # Make sure at least one observation from each class; otherwise, toss 
        # out this simulation
        all_classes <- (length(unique(y)) == K)

        stopifnot(length(unique(y)) <= K)

        counter <- counter + 1
    }

    # y <- as.factor(y)
    y <- factor(y, levels=1:K, ordered=TRUE)

    
    return(list(X=X, y=y, probs=probs, beta_final=beta_mat))
}

sim_prop_odds <- function(n, p, K, intercepts, beta, nsim){

    stopifnot(K >= 3)
    stopifnot(length(intercepts) == K - 1)
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))

    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) == p)

    # Generate beta_mat

    beta_mat <- matrix(0, nrow=p, ncol=K - 1)
    for(k in 1:(K-1)){
        beta_mat[, k] <- beta
    }

    # List we'll return: will be length nsim, and every element will be a
    # named list with elements X and y.
    ret_list <- list()

    for(i in 1:nsim){
        ret_list[[i]] <- sim_data(n, p, K, intercepts, beta_mat)
    }        
    
    return(ret_list)
}

sim_relaxed_prop_odds <- function(n, p, K, intercepts, beta, nsim, dev_num,
    dev_size=0){

    stopifnot(K >= 3)
    stopifnot(length(intercepts) == K - 1)
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))

    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) == p)

    # List we'll return: will be length nsim, and every element will be a
    # named list with elements X and y.
    ret_list <- list()

    for(i in 1:nsim){

        # If beta_mat that is generated results in negative probabilities,
        # try again with another random draw, up to 10 times before giving up

        ret_i <- NA
        n_tries <- 0

        while(any(is.na(ret_i))){

            if(n_tries >= 10){
                print(paste("p =", p, ", K =", K, "dev_num =", dev_num,
                    "dev_size =", dev_size))
                print("intercepts:")
                print(intercepts)
                print("colMeans(probs):")
                print(colMeans(probs))
                stop("At least one probability was negative")
            }

            # Generate beta_mat

            beta_mat <- matrix(0, nrow=p, ncol=K - 1)
            beta_mat[, 1] <- beta

            # Generate sparse random deviations
            dev_mat <- matrix(0, nrow=p, ncol=K - 2)

            for(k in 1:(K - 2)){
                # Select random indices for deviations
                inds <- sample(p, dev_num)
                dev_mat[inds, k] <- dev_size
            }

            # Create random sign flips to multiply deviations by
            sign_mat_content <- sample(c(-1, 1), size=p*(K - 2), replace=TRUE)
            sign_mat <- matrix(sign_mat_content, nrow=p, ncol=K - 2)

            # Multiply deviations by random sign flips
            dev_mat <- dev_mat * sign_mat

            # Add sparse random deviations to beta_mat
            for(k in 2:(K-1)){
                beta_mat[, k] <- beta_mat[, k - 1] + dev_mat[, k - 1]
            }

            # Generate X, y, probabilities
            ret_i <- sim_data(n, p, K, intercepts, beta_mat, dev_num,
                dev_size)

            n_tries <- n_tries + 1

        }

        ret_list[[i]] <- ret_i
        
    }
    
    return(ret_list)
}

sim_relaxed_prop_odds_rand <- function(n, p, K, intercepts, beta, nsim,
    dev_prob, dev_size=0){

    stopifnot(K >= 3)
    stopifnot(length(intercepts) == K - 1)
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))

    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) == p)

    stopifnot(is.numeric(dev_prob) | is.integer(dev_prob))
    stopifnot(!is.na(dev_prob))
    stopifnot(length(dev_prob) == 1)
    stopifnot(dev_prob >= 0)
    stopifnot(dev_prob <= 1)


    # List we'll return: will be length nsim, and every element will be a
    # named list with elements X and y.
    ret_list <- list()

    for(i in 1:nsim){

        # If beta_mat that is generated results in negative probabilities,
        # try again with another random draw, up to 10 times before giving up

        ret_i <- NA
        n_tries <- 0

        while(any(is.na(ret_i))){

            if(n_tries >= 100){
                print(paste("p =", p, ", K =", K, "dev_num =", dev_num,
                    "dev_size =", dev_size))
                print("intercepts:")
                print(intercepts)
                print("colMeans(probs):")
                print(colMeans(probs))
                stop("At least one probability was negative")
            }

            # Generate beta_mat

            beta_mat <- matrix(0, nrow=p, ncol=K - 1)
            beta_mat[, 1] <- beta

            # Generate sparse random deviations
            dev_mat <- matrix(0, nrow=p, ncol=K - 2)

            for(k in 1:(K - 2)){
                # Select random indices for deviations
                inds <- which(as.logical(rbinom(p, size=1, prob=dev_prob)))
                dev_mat[inds, k] <- dev_size
            }

            # Create random sign flips to multiply deviations by
            sign_mat_content <- sample(c(-1, 1), size=p*(K - 2), replace=TRUE)
            sign_mat <- matrix(sign_mat_content, nrow=p, ncol=K - 2)

            # Multiply deviations by random sign flips
            dev_mat <- dev_mat * sign_mat

            # Add sparse random deviations to beta_mat
            for(k in 2:(K-1)){
                beta_mat[, k] <- beta_mat[, k - 1] + dev_mat[, k - 1]
            }

            # Generate X, y, probabilities
            ret_i <- sim_data(n, p, K, intercepts, beta_mat, dev_num,
                dev_size)

            n_tries <- n_tries + 1

        }

        ret_list[[i]] <- ret_i
        
    }
    
    return(ret_list)
}

sim_relaxed_gaussian_prop_odds <- function(n, p, K, intercepts, beta, nsim,
    dev_var=1/3){

    stopifnot(K >= 3)
    stopifnot(length(intercepts) == K - 1)
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))

    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) == p)

    # stopifnot(is.numeric(dev_prob) | is.integer(dev_prob))
    # stopifnot(!is.na(dev_prob))
    # stopifnot(length(dev_prob) == 1)
    # stopifnot(dev_prob >= 0)
    # stopifnot(dev_prob <= 1)


    # List we'll return: will be length nsim, and every element will be a
    # named list with elements X and y.
    ret_list <- list()

    # dev_num <- round(sqrt(p))

    for(i in 1:nsim){

        # If beta_mat that is generated results in negative probabilities,
        # try again with another random draw, up to 10 times before giving up

        ret_i <- NA
        n_tries <- 0

        while(any(is.na(ret_i))){

            if(n_tries >= 100){
                print(paste("p =", p, ", K =", K))
                print("intercepts:")
                print(intercepts)
                stop("At least one probability was negative")
            }

            # Generate beta_mat

            beta_mat <- matrix(0, nrow=p, ncol=K - 1)
            beta_mat[, 1] <- beta

            # Generate Gaussian deviations
            dev_mat <- matrix(rnorm(p*(K - 2), mean=0, sd=sqrt(dev_var)),
                nrow=p, ncol=K - 2)

            # Add Gaussian deviations to beta_mat
            for(k in 2:(K-1)){
                beta_mat[, k] <- beta_mat[, k - 1] + dev_mat[, k - 1]
            }

            # Generate X, y, probabilities
            ret_i <- sim_data(n, p, K, intercepts, beta_mat)

            n_tries <- n_tries + 1

            ret_list[[i]] <- ret_i

        }

    }
    
    return(ret_list)
}

sim_relaxed_unif_prop_odds <- function(n, p, K, intercepts, beta, nsim,
    dev_size){

    stopifnot(K >= 3)
    stopifnot(length(intercepts) == K - 1)
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))

    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) == p)

    # stopifnot(is.numeric(dev_prob) | is.integer(dev_prob))
    # stopifnot(!is.na(dev_prob))
    # stopifnot(length(dev_prob) == 1)
    # stopifnot(dev_prob >= 0)
    # stopifnot(dev_prob <= 1)


    # List we'll return: will be length nsim, and every element will be a
    # named list with elements X and y.
    ret_list <- list()

    # dev_num <- round(sqrt(p))

    for(i in 1:nsim){

        # If beta_mat that is generated results in negative probabilities,
        # try again with another random draw, up to 10 times before giving up

        ret_i <- NA
        n_tries <- 0

        while(any(is.na(ret_i))){

            if(n_tries >= 100){
                print(paste("p =", p, ", K =", K))
                print("intercepts:")
                print(intercepts)
                stop("At least one probability was negative")
            }

            # Generate beta_mat

            beta_mat <- matrix(0, nrow=p, ncol=K - 1)
            beta_mat[, 1] <- beta

            # Generate uniform deviations
            dev_mat <- matrix(runif(p*(K - 2), min=-dev_size, max=dev_size),
                nrow=p, ncol=K - 2)

            # Add uniform deviations to beta_mat
            for(k in 2:(K-1)){
                beta_mat[, k] <- beta_mat[, k - 1] + dev_mat[, k - 1]
            }

            # Generate X, y, probabilities
            ret_i <- sim_data(n, p, K, intercepts, beta_mat)

            n_tries <- n_tries + 1

            ret_list[[i]] <- ret_i

        }

    }
    
    return(ret_list)
}

prop_odds_model <- function(n, p, K, intercepts, beta){
    # Number of classes
    stopifnot(is.numeric(K) | is.integer(K))
    stopifnot(K == round(K))
    stopifnot(K >= 3)
    stopifnot(is.numeric(intercepts) | is.integer(intercepts))
    stopifnot(length(intercepts) == K - 1)
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))
    stopifnot(is.numeric(p) | is.integer(p))
    stopifnot(p == round(p))
    stopifnot(p >= 1)
    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) >= p)
    if(length(beta) > p){
        beta <- beta[1:p]
    }
    
    my_model <- new_model(name = "prop_odds_mod", 
                label = sprintf("Proportional odds model (n = %s,
                    p = %s, K = %s, rare intercept = %s)",
                    n, p, K, intercepts[K-1]),
                params = list(n = n, p = p, K = K, intercepts=intercepts,
                    beta = beta, dev_size=0),
                simulate = sim_prop_odds
    )
    return(my_model)
}

relax_prop_odds_model <- function(n, p, K, intercepts, beta, dev_size,
    dev_num){
    # Number of classes
    stopifnot(is.numeric(K) | is.integer(K))
    stopifnot(K == round(K))
    stopifnot(K >= 3)
    stopifnot(is.numeric(intercepts) | is.integer(intercepts))
    stopifnot(length(intercepts) >= K - 1)
    if(length(intercepts) > K - 1){
        intercepts <- intercepts[1:(K - 1)]
    }
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))
    stopifnot(is.numeric(p) | is.integer(p))
    stopifnot(p == round(p))
    stopifnot(p >= 1)
    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) >= p)
    if(length(beta) > p){
        beta <- beta[1:p]
    }

    stopifnot(is.numeric(dev_num) | is.integer(dev_num))
    stopifnot(!is.na(dev_num))
    stopifnot(length(dev_num) == 1)
    stopifnot(dev_num == round(dev_num))
    stopifnot(dev_num <= p)
    
    my_model <- new_model(name = "relax_prop_odds_mod", 
        label = sprintf("Relaxed proportional odds model (n = %s,
            p = %s, K = %s, rare intercept = %s, deviation size = %s)",
            n, p, K, intercepts[K-1], dev_size),
        params = list(n = n, p = p, K = K, intercepts=intercepts,
            beta = beta, dev_num=dev_num, dev_size = dev_size),
        simulate = sim_relaxed_prop_odds
    )
    return(my_model)
}

relax_prop_odds_model_rand <- function(n, p, K, intercepts, beta, dev_size,
    dev_prob){
    # Number of classes
    stopifnot(is.numeric(K) | is.integer(K))
    stopifnot(K == round(K))
    stopifnot(K >= 3)
    stopifnot(is.numeric(intercepts) | is.integer(intercepts))
    stopifnot(length(intercepts) >= K - 1)
    if(length(intercepts) > K - 1){
        intercepts <- intercepts[1:(K - 1)]
    }
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))
    stopifnot(is.numeric(p) | is.integer(p))
    stopifnot(p == round(p))
    stopifnot(p >= 1)
    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) >= p)
    if(length(beta) > p){
        beta <- beta[1:p]
    }

    stopifnot(is.numeric(dev_prob) | is.integer(dev_prob))
    stopifnot(!is.na(dev_prob))
    stopifnot(length(dev_prob) == 1)
    stopifnot(dev_prob >= 0)
    stopifnot(dev_prob <= 1)
    
    my_model <- new_model(name = "relax_prop_odds_mod", 
        label = sprintf("Relaxed proportional odds model (n = %s,
            p = %s, K = %s, rare intercept = %s, deviation size = %s)",
            n, p, K, intercepts[K-1], dev_size),
        params = list(n = n, p = p, K = K, intercepts=intercepts,
            beta = beta, dev_prob=dev_prob, dev_size = dev_size),
        simulate = sim_relaxed_prop_odds_rand
    )
    return(my_model)
}

relax_prop_odds_gaussian_model <- function(n, p, K, intercepts, beta, dev_var){
    # Number of classes
    stopifnot(is.numeric(K) | is.integer(K))
    stopifnot(K == round(K))
    stopifnot(K >= 3)
    stopifnot(is.numeric(intercepts) | is.integer(intercepts))
    stopifnot(length(intercepts) >= K - 1)
    if(length(intercepts) > K - 1){
        intercepts <- intercepts[1:(K - 1)]
    }
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))
    stopifnot(is.numeric(p) | is.integer(p))
    stopifnot(p == round(p))
    stopifnot(p >= 1)
    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) >= p)
    if(length(beta) > p){
        beta <- beta[1:p]
    }

    # stopifnot(is.numeric(dev_prob) | is.integer(dev_prob))
    # stopifnot(!is.na(dev_prob))
    # stopifnot(length(dev_prob) == 1)
    # stopifnot(dev_prob >= 0)
    # stopifnot(dev_prob <= 1)

    stopifnot(is.numeric(dev_var) | is.integer(dev_var))
    stopifnot(!is.na(dev_var))
    stopifnot(length(dev_var) == 1)
    stopifnot(dev_var > 0)
    
    my_model <- new_model(name = "relax_prop_odds_mod_gauss", 
        label = sprintf("Relaxed proportional odds model with Gaussian noise (n = %s,
            p = %s, K = %s, rare intercept = %s, deviation variance = %s)",
            n, p, K, intercepts[K-1], dev_var),
        params = list(n = n, p = p, K = K, intercepts=intercepts,
            beta = beta, dev_var=dev_var),
        simulate = sim_relaxed_gaussian_prop_odds
    )
    return(my_model)
}

relax_prop_odds_unif_model <- function(n, p, K, intercepts, beta, dev_size){
    # Number of classes
    stopifnot(is.numeric(K) | is.integer(K))
    stopifnot(K == round(K))
    stopifnot(K >= 3)
    stopifnot(is.numeric(intercepts) | is.integer(intercepts))
    stopifnot(length(intercepts) >= K - 1)
    if(length(intercepts) > K - 1){
        intercepts <- intercepts[1:(K - 1)]
    }
    stopifnot(all(intercepts == sort(intercepts, decreasing=FALSE)))
    stopifnot(is.numeric(p) | is.integer(p))
    stopifnot(p == round(p))
    stopifnot(p >= 1)
    stopifnot(is.numeric(beta) | is.integer(beta))
    stopifnot(length(beta) >= p)
    if(length(beta) > p){
        beta <- beta[1:p]
    }

    # stopifnot(is.numeric(dev_prob) | is.integer(dev_prob))
    # stopifnot(!is.na(dev_prob))
    # stopifnot(length(dev_prob) == 1)
    # stopifnot(dev_prob >= 0)
    # stopifnot(dev_prob <= 1)

    stopifnot(is.numeric(dev_size) | is.integer(dev_size))
    stopifnot(!is.na(dev_size))
    stopifnot(length(dev_size) == 1)
    stopifnot(dev_size > 0)
    
    my_model <- new_model(name = "relax_prop_odds_mod_unif", 
        label = sprintf("Relaxed proportional odds model with uniform noise (n = %s,
            p = %s, K = %s, rare intercept = %s, deviation size = %s)",
            n, p, K, intercepts[K-1], dev_size),
        params = list(n = n, p = p, K = K, intercepts=intercepts,
            beta = beta, dev_size=dev_size),
        simulate = sim_relaxed_unif_prop_odds
    )
    return(my_model)
}

gen_sim_data_app_freq <- function(df, test_prop, resp_names, final_resp_names, nsim){

    # Check inputs

    stopifnot(is.data.frame(df))

    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop > 0)
    stopifnot(test_prop < 1)

    stopifnot(is.character(resp_names))
    stopifnot(all(!is.na(resp_names)))
    stopifnot(all(resp_names %in% colnames(df)))

    K <- length(resp_names)

    stopifnot(K >= 3)

    stopifnot(is.character(final_resp_names))
    stopifnot(all(!is.na(final_resp_names)))
    stopifnot(length(final_resp_names) == K)


    # Create design matrix and response matrix for ordinalNet

    ret_list <- prepOrdNet(df=df, final_y_names=final_resp_names,
         cum_y_names=resp_names)

    X_ordnet <- ret_list$X_ordnet
    y_mat <- ret_list$y_mat

    rm(ret_list)

    p <- ncol(X_ordnet)

    n_obs <- sum(y_mat)

    # Split y_mat into training and test sets

    n_test <- round(test_prop*n_obs)

    stopifnot(n_test >= 1)
    stopifnot(n_test < n_obs)

    stopifnot(is.matrix(X_ordnet))
    stopifnot(ncol(X_ordnet) == p)
    stopifnot(all(!is.na(X_ordnet)))

    ret <- list()

    for(i in 1:nsim){

        # Want to ensure that every column has at least one observation (the
        # function subsampCountMat ensures this already for y_mat_test, so just
        # need to confirm for y_mat_train)
        all_cols_flag <- FALSE
        iter_count <- 0

        while(!all_cols_flag){
              
            if(iter_count >= 10){
                stop("Unable to create subsampled matrix with at least one observation in each column after 10 attempts; breaking loop here.")
            }

            y_mat_test <- subsampCountMat(count_mat=y_mat, sample_size=n_test)

            y_mat_train <- y_mat - y_mat_test

            stopifnot(all(y_mat_train <= y_mat))
            stopifnot(all(y_mat_train >= 0))
            stopifnot(all(dim(y_mat_train) == dim(y_mat)))
            stopifnot(all(colnames(y_mat_train) == colnames(y_mat)))
            stopifnot(sum(y_mat_train) == n_obs - n_test)

            all_cols_flag <- all(colSums(y_mat_train) > 0)
            iter_count <- iter_count + 1
        }

        # Create data.frame for MASS::polr function

        ret_list <- pre_df_polr(y_mat=y_mat_train, X_ordnet=X_ordnet)

        df_polr_train <- ret_list$df_polr
        weights_polr_train <- ret_list$weights_polr

        rm(ret_list)

        ret_list2 <- pre_df_polr(y_mat=y_mat_test, X_ordnet=X_ordnet)

        df_polr_test <- ret_list2$df_polr
        weights_polr_test <- ret_list2$weights_polr

        rm(ret_list2)

        # Check outputs

        stopifnot(ncol(df_polr_test) == ncol(df_polr_train))
        stopifnot(ncol(df_polr_train) == p + 1)

        ret[[i]] <- list(train_df_polr=df_polr_train,
            train_weights_polr=weights_polr_train, X_ordnet=X_ordnet,
            train_ymat=y_mat_train, test_df_polr=df_polr_test,
            test_weights_polr=weights_polr_test, test_ymat=y_mat_test)
    }

    return(ret)
}

sim_data_app <- function(df, test_prop, resp_names, final_resp_names, nsim){

    # Check inputs

    stopifnot(is.data.frame(df))

    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop > 0)
    stopifnot(test_prop < 1)

    stopifnot(is.character(resp_names))
    stopifnot(all(!is.na(resp_names)))
    stopifnot(all(resp_names %in% colnames(df)))

    K <- length(resp_names)

    stopifnot(K >= 3)

    stopifnot(is.character(final_resp_names))
    stopifnot(all(!is.na(final_resp_names)))
    stopifnot(length(final_resp_names) == K)


    # Create design matrix and response matrix for ordinalNet

    ret_list <- prepOrdNet(df=df, final_y_names=final_resp_names,
         cum_y_names=resp_names)

    X_ordnet <- ret_list$X_ordnet
    y_mat <- ret_list$y_mat

    rm(ret_list)

    p <- ncol(X_ordnet)

    n_obs <- sum(y_mat)

    # Split y_mat into training and test sets

    n_test <- round(test_prop*n_obs)

    stopifnot(n_test >= 1)
    stopifnot(n_test < n_obs)

    stopifnot(is.matrix(X_ordnet))
    stopifnot(ncol(X_ordnet) == p)
    stopifnot(all(!is.na(X_ordnet)))

    ret <- list()

    for(i in 1:nsim){

        # Want to ensure that every column has at least one observation (the
        # function subsampCountMat ensures this already for y_mat_test, so just
        # need to confirm for y_mat_train)
        all_cols_flag <- FALSE
        iter_count <- 0

        while(!all_cols_flag){
              
            if(iter_count >= 10){
                stop("Unable to create subsampled matrix with at least one observation in each column after 10 attempts; breaking loop here.")
            }

            y_mat_test <- subsampCountMat(count_mat=y_mat, sample_size=n_test)

            y_mat_train <- y_mat - y_mat_test

            stopifnot(all(y_mat_train <= y_mat))
            stopifnot(all(y_mat_train >= 0))
            stopifnot(all(dim(y_mat_train) == dim(y_mat)))
            stopifnot(all(colnames(y_mat_train) == colnames(y_mat)))
            stopifnot(sum(y_mat_train) == n_obs - n_test)

            all_cols_flag <- all(colSums(y_mat_train) > 0)
            iter_count <- iter_count + 1
        }

        # Create data.frame for MASS::polr function

        ret_list <- pre_df_polr(y_mat=y_mat_train, X_ordnet=X_ordnet)

        df_polr_train <- ret_list$df_polr
        weights_polr_train <- ret_list$weights_polr

        rm(ret_list)

        ret_list2 <- pre_df_polr(y_mat=y_mat_test, X_ordnet=X_ordnet)

        df_polr_test <- ret_list2$df_polr
        weights_polr_test <- ret_list2$weights_polr

        rm(ret_list2)

        # Check outputs

        stopifnot(ncol(df_polr_test) == ncol(df_polr_train))
        stopifnot(ncol(df_polr_train) == p + 1)

        ret[[i]] <- list(train_df_polr=df_polr_train,
            train_weights_polr=weights_polr_train, X_ordnet=X_ordnet,
            train_ymat=y_mat_train, test_df_polr=df_polr_test,
            test_weights_polr=weights_polr_test, test_ymat=y_mat_test)
    }

    return(ret)
}

data_app_model <- function(df, test_prop, resp_names, final_resp_names){

    stopifnot(is.data.frame(df))
    
    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop > 0)
    stopifnot(test_prop < 1)
    
    my_model <- new_model(name = "dat_app_mod", 
                label = "Data Application",
                params = list(df=df, test_prop=test_prop, resp_names=resp_names,
                    final_resp_names=final_resp_names),
                simulate = sim_data_app
    )

    return(my_model)
}

gen_data_app_model <- function(df, test_prop, resp_names, final_resp_names){

    stopifnot(is.data.frame(df))
    
    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop > 0)
    stopifnot(test_prop < 1)
    
    my_model <- new_model(name = "dat_app_mod", 
                label = "Data Application",
                params = list(df=df, test_prop=test_prop, resp_names=resp_names,
                    final_resp_names=final_resp_names),
                simulate = gen_sim_data_app_freq
    )

    return(my_model)
}

gen_sim_data_app <- function(df, test_prop, resp_name, resp_levels, x_formula,
    omitted_x, nsim){

    # Check inputs

    stopifnot(is.data.frame(df))

    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop > 0)
    stopifnot(test_prop < 1)

    stopifnot(is.character(resp_levels))
    stopifnot(all(!is.na(resp_levels)))

    stopifnot(is.character(resp_name))
    stopifnot(length(resp_name) == 1)
    stopifnot(resp_name %in% colnames(df))

    K <- length(resp_levels)

    stopifnot(K >= 3)

    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop >= 0)
    stopifnot(test_prop < 1)

    if(any(!is.na(omitted_x))){
        stopifnot(is.character(omitted_x))
        stopifnot(omitted_x %in% colnames(df))
    } else{
        stopifnot(all(is.na(omitted_x)))
    }
    

    # Create design matrix and response matrix for ordinalNet

    X_ordnet <- model.matrix(x_formula, data=df)

    stopifnot(nrow(X_ordnet) == nrow(df))

    X_ordnet <- X_ordnet[, colnames(X_ordnet) != "(Intercept)"]

    p <- ncol(X_ordnet)

    n_obs <- nrow(X_ordnet)

    stopifnot(is.matrix(X_ordnet))
    stopifnot(all(!is.na(X_ordnet)))
    stopifnot(n_obs == nrow(df))

    ret <- list()

    for(i in 1:nsim){

        # Split y_mat into training and test sets

        ret_list <- splitIntoTwoSetsVec(X=X_ordnet, y=df[, resp_name],
            prop=test_prop, n_attempts=100)

        X_ordnet_test <- ret_list$X_prop

        X_ordnet_train <- ret_list$X_comp

        y_test <- ret_list$y_prop

        y_train <- ret_list$y_comp

        train_inds <- ret_list$comp_inds

        rm(ret_list)

        # Create data.frame for MASS::polr function

        if(all(!is.na(omitted_x))){
            df <- df[, !(colnames(df) %in% omitted_x)]
        }
        
        df_polr_train <- df[train_inds, ]
        # No need for test df--will just extract coefficients

        # Check outputs

        stopifnot(ncol(X_ordnet_train) == ncol(X_ordnet_test))

        stopifnot(nrow(X_ordnet_train) == length(y_train))
        stopifnot(nrow(df_polr_train) == length(y_train))

        stopifnot(nrow(X_ordnet_test) == length(y_test))

        stopifnot(is.factor(y_train))
        stopifnot(is.factor(y_test))

        stopifnot(length(unique(y_train)) == K)
        stopifnot(length(unique(y_test)) == K)

        ret[[i]] <- list(train_df_polr=df_polr_train,
            X_ordnet_train=X_ordnet_train, train_y=y_train,
            test_y=y_test, X_ordnet_test=X_ordnet_test)
    }

    return(ret)
}

gen_data_app_model <- function(df, tune_prop, test_prop, resp_name,
    resp_levels, x_formula, omitted_x){

    stopifnot(is.data.frame(df))
    stopifnot(all(!is.na(df)))

    stopifnot(is.numeric(test_prop) | is.integer(test_prop))
    stopifnot(length(test_prop) == 1)
    stopifnot(test_prop >= 0)
    stopifnot(test_prop < 1)

    stopifnot(is.numeric(tune_prop) | is.integer(tune_prop))
    stopifnot(length(tune_prop) == 1)
    stopifnot(tune_prop >= 0)
    stopifnot(tune_prop < 1)

    stopifnot(is.character(resp_levels))
    stopifnot(all(!is.na(resp_levels)))

    K <- length(resp_levels)

    stopifnot(K >= 3)

    stopifnot(is.character(resp_name))
    stopifnot(length(resp_name) == 1)
    stopifnot(resp_name %in% colnames(df))

    if(any(!is.na(omitted_x))){
        stopifnot(is.character(omitted_x))
        stopifnot(omitted_x %in% colnames(df))
    } else{
        stopifnot(all(is.na(omitted_x)))
    }
    
    my_model <- new_model(name = "dat_app_mod", 
                label = "Data Application",
                params = list(df=df, test_prop=test_prop, resp_name=resp_name,
                    resp_levels=resp_levels, tune_prop=tune_prop,
                    x_formula=x_formula, omitted_x=omitted_x),
                simulate = gen_sim_data_app
    )

    return(my_model)
}

