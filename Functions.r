## Collection of functions for the scripts of Project Hydropower

# Check whether data has correct form
check_riv_data <-  function(riv){
  
  # riv can be output of getRiver function or a similar list
  
  if(is.data.frame(riv$d)==FALSE){
    stop("please make sure data is a dataframe")
  }
  if(ncol(riv$d)!=2){
    stop("data has incorrect shape. make sure data has 2 columns")
  }
  if(inherits(riv$d$time, "Date")==FALSE){
    stop("please make sure columns <date> has Date-format")
  }
  print("data is in correct form")
}

longest_coh_ts <- function(d.ts){
  stopifnot(is.data.frame(d.ts),
            length(d.ts$time) == nrow(d.ts))
  d.ts <- na.omit(d.ts)
  
  n <- nrow(d.ts)
  gaps <- as.numeric(d.ts$time[-1] - d.ts$time[-n], units="days")
  
  index <- which(gaps > 1)
  
  window <- c(0, diff(index))
  i.max <- which.max(window)
  end  <-  index[i.max]
  start <- end - window[i.max] + 1
  
  d.ts[start:end, ]
}

# Produce beautiful time series plots
plott <- function(df, type="l", reset_par=TRUE, log="y", ylab = expression("discharge" * (m^3/sec)), main = NULL, sub=NULL, 
                  mar = c(4,4, 1+3*hasTit, 5) + 0.1 + marP, marP = rep(0,4), mgp = c(2, 0.6, 0),
                  col.grid  = "lightgray", lty.grid = "dotted", ...)
{
  hasTit <- !is.null(main) && is.character(main) && nzchar(main)
  opar <- par(mar=mar, mgp=mgp)
  if (reset_par) on.exit(par(opar))
  plot(df, type=type, log=log, main=main, sub=sub, ylab=ylab, yaxt = "n", ...)
  ## y-axis :
  aty <- axTicks(2) # logarithmic
  eaxis(2, at=aty, sub10=2)
  ## x-axis = time + grid lines every year 'beg.yr'
  nt <- as.numeric(tt <- df$time)
  ltt <- as.POSIXlt(tt)
  Jan.1st <- ltt $ yday == 0 # yday == 0:  January 1st
  beg.yr <- nt[Jan.1st]
  axis(1, at = beg.yr, labels = 1900 + ltt$year[Jan.1st]) # or eaxis()  [in new version]
  ## grid lines for both axes (h , v) :
  abline(h = aty,
         v = beg.yr, col=col.grid, lty=lty.grid)
}

#Transform extraction data from array do dataframe
meteo2df <- function(met) {
  within(as.data.frame.table(met$arr, responseName = met$var), {
    lon <- as.numeric(as.character(lon))
    lat <- as.numeric(as.character(lat))
    time <- met$dates # == as.Date(as.character(time))
  })
}

#Read all of extracted meteo .rds files from data directory and store them in a list
read_meteo <- function(){
  files <- list.files(pattern = paste0("^", "Extracted_"),
                      full.names = TRUE)
  collected_files <- list()
  for(i in files){
    meteo <- readRDS(i)
    var <- meteo$var
    collected_files[[var]] <- meteo
  }
  return(collected_files)
}

# Extract individual ts
get_ts_ind <- function(data, lon, lat, start_time="1979-01-01", 
                       end_time="2018-12-31", range=TRUE, pool=FALSE, 
                       pool_grid=c(3,3), pool_fun=c(mean, max), make_df=TRUE){
  
  #check if lon/lat in range
  if(any(lon<data$lon_range[1])||any(lon>data$lon_range[2])||
     any(lat<data$lat_range[1])||any(lat>data$loat_range[2])){
    stop("lon or lat coordinates out of range. lon range:", 
         paste(data$lon_range, collapse = "-"), 
         ". lat range:", paste(data$lon_range, collapse = "-"), ".")
  }
  array <- aperm(data$arr, c(2,1,3)) #arrange array lat*lon*time to match earth grid
  if(range==TRUE){
    stopifnot(length(lon) == 2, is.numeric(lon), abs(diff(lon)) >= 0
              && length(lat) == 2, is.numeric(lat), abs(diff(lat)) >= 0)
    if(lon[1]>lon[2]||lat[1]<lat[2]){
      stop("lon-vector should stretch west-east and
           lat-vector should stretch north-south.")
    }
    ilon1 <- match(lon[1], colnames(array))
    ilon2 <- match(lon[2], colnames(array))
    ilat1 <- match(lat[1], rownames(array))
    ilat2 <- match(lat[2], rownames(array))
    itime1 <- match(start_time, dimnames(array)$time)
    itime2 <- match(end_time, dimnames(array)$time)
    sliced_data <- array[ilat1:ilat2,ilon1:ilon2,itime1:itime2]
  }else{
    stopifnot(length(lon)==length(lat))
    ilons <- match(lon,colnames(array))
    ilats <- match(lat, rownames(array))
    itime1 <- match(start_time, dimnames(array)$time)
    itime2 <- match(end_time, dimnames(array)$time)
    sliced_data <- array[ilats,ilons,itime1:itime2]
  }
  
  # start of pooling operations
  if(pool==TRUE){
    if(length(pool_grid)!=2){
      stop("Please enter dimensions of pool_grid as a 2-dimensional vector")
    }
    # initialize pooled array
    pooled <- array(NA_real_, dim=c(pool_grid,length(dimnames(sliced_data)$time)))
    
    # check if components of pool_grid are ==1
    # else split lat/lon
    if(pool_grid[1]==1){
      lat_parts <- list(lat=dimnames(sliced_data)$lat)
    }else{
      lat_parts <- split_vec(dimnames(sliced_data)$lat, pool_grid[1])
    }
    if(pool_grid[2]==1){
      lon_parts <- list(lat=dimnames(sliced_data)$lon)
    }else{
      lon_parts <- split_vec(dimnames(sliced_data)$lon, pool_grid[2])
    }
    # filling pooled array
    p.lat.ind <- 0
    for(i in lat_parts){
      p.lon.ind <- 0
      p.lat.ind <- p.lat.ind + 1
      for(j in lon_parts){
        p.lon.ind <- p.lon.ind + 1
        p.ilons <- match(j,colnames(sliced_data))
        p.ilats <- match(i, rownames(sliced_data))
        pre_pool <- sliced_data[p.ilats,p.ilons,]
        pooled[p.lat.ind,p.lon.ind,] <- apply(pre_pool, 3, pool_fun)
      }
    }
    
    lat_parts_name <- c()
    for(i in lat_parts){
      name <- paste(c(i[1], i[length(i)]), collapse = "-")
      lat_parts_name <- c(lat_parts_name, name)
    }
    
    lon_parts_name <- c()
    for(i in lon_parts){
      name <- paste(c(i[1], i[length(i)]), collapse = "-")
      lon_parts_name <- c(lon_parts_name, name)
    }
    
    dimnames(pooled) = list(lat = lat_parts_name,
                            lon = lon_parts_name,
                            time = as.character(dimnames(sliced_data)$time))
  }
  
  # turn array into df's
  if(make_df == TRUE){
    df_full <-  data.frame(index=as.Date(dimnames(sliced_data)$time))
    for(i in dimnames(sliced_data)$lat){
      for(j in dimnames(sliced_data)$lon){
        df_full[paste(i,j , collapse="-")]=unname(sliced_data[i,j,])
      }
    }
    if(pool == TRUE){
      df_pooled <-  data.frame(index=as.Date(dimnames(pooled)$time))
      for(i in dimnames(pooled)$lat){
        for(j in dimnames(pooled)$lon){
          df_pooled[paste(i,j , collapse="-")]=unname(pooled[i,j,])
        }
      }
    }
    
  }
  
  if(make_df==TRUE && pool == TRUE){
    return(list(var=data$var, array=sliced_data, df_full=df_full,
                pooled_array=pooled, df_pooled=df_pooled))
  }else if(make_df==FALSE && pool == FALSE){
    return(list(var=data$var, array=sliced_data))
  }else if(make_df==TRUE && pool == FALSE){
    return(list(var=data$var,array=sliced_data, df_full=df_full))
  }else if(make_df==FALSE && pool == TRUE){
    return(list(var=data$var,array=sliced_data,
                pooled_array=pooled))
  }
  
}

# Extract all necessary covariate time series
get_ts_all <- function(data, lon, lat, start_time="1979-01-01", 
                       end_time="2018-12-31", range=TRUE, pool=FALSE, 
                       pool_grid=c(3,3), pool_fun=c(mean, max), make_matrix=FALSE,
                       lag=FALSE, lags =c(1,2,3), lag_pool=FALSE, lag_pool_span=5, standardize=TRUE){
  
  # Takes as input whole list of read_meteo(). Outputs a list 
  
  
  # initiate final list which holds smaller lists, one per variable
  final_list <- list()
  
  #start loop over variables
  for(i in data){
    # initiate variable specific list to merge to final list in the end
    variable_list <- list()
    
    #check if lon/lat in range
    if(any(lon<i$lon_range[1])||any(lon>i$lon_range[2])||
       any(lat<i$lat_range[1])||any(lat>i$loat_range[2])){
      stop("lon or lat coordinates out of range. lon range:", 
           paste(i$lon_range, collapse = "-"), 
           ". lat range:", paste(i$lon_range, collapse = "-"), ".")
    } #end check
    
    print(paste(i$var, ": check done"))
    
    # get date-sliced array
    array <- aperm(i$arr, c(2,1,3)) #arrange array lat*lon*time to match earth grid
    if(range==TRUE){
      stopifnot(length(lon) == 2, is.numeric(lon), abs(diff(lon)) >= 0
                && length(lat) == 2, is.numeric(lat), abs(diff(lat)) >= 0)
      if(lon[1]>lon[2]||lat[1]<lat[2]){
        stop("lon-vector should stretch west-east and
           lat-vector should stretch north-south.")
      }
      ilon1 <- match(lon[1], colnames(array))
      ilon2 <- match(lon[2], colnames(array))
      ilat1 <- match(lat[1], rownames(array))
      ilat2 <- match(lat[2], rownames(array))
      itime1 <- match(start_time, dimnames(array)$time)
      itime2 <- match(end_time, dimnames(array)$time)
      sliced_data <- array[ilat1:ilat2,ilon1:ilon2,itime1:itime2]
    }else{
      stopifnot(length(lon)==length(lat))
      ilons <- match(lon,colnames(array))
      ilats <- match(lat, rownames(array))
      itime1 <- match(start_time, dimnames(array)$time)
      itime2 <- match(end_time, dimnames(array)$time)
      sliced_data <- array[ilats,ilons,itime1:itime2]
    } # end unpooled 
    
    if(i$var == "Precipitation_Flux"){
      sliced_data <- log(sliced_data+0.0001, base = 10)
    }
    
    if(standardize==TRUE){
      for(u in 1:dim(sliced_data)[1]){
        for(v in 1:dim(sliced_data)[2])
          sliced_data[u,v,] <- scale(sliced_data[u,v,])
      }
    }
    attr(sliced_data, "pooled") <- "unpooled"
    attr(sliced_data, "lag_pooled") <- FALSE
    attr(sliced_data, "name") <- paste0(i$var, "_unpooled_arr")
    variable_list[["var_unpooled_arr"]] <- sliced_data
    
    print(paste(i$var, ": sliced_data done"))
    
    #define some general metrics
    len <- dim(sliced_data)[3]
    time_dim <- dimnames(sliced_data)$time
    
    # start of pooling operations
    if(pool==TRUE){
      if(length(pool_grid)!=2){
        stop("Please enter dimensions of pool_grid as a 2-dimensional vector")
      }
      
      # initialize pooled array
      pooled <- array(NA_real_, dim=c(pool_grid, len))
      # check if components of pool_grid are ==1
      # else split lat/lon
      if(pool_grid[1]==1){
        lat_parts <- list(lat=dimnames(sliced_data)$lat)
      }else{
        lat_parts <- split_vec(dimnames(sliced_data)$lat, pool_grid[1])
      }
      if(pool_grid[2]==1){
        lon_parts <- list(lat=dimnames(sliced_data)$lon)
      }else{
        lon_parts <- split_vec(dimnames(sliced_data)$lon, pool_grid[2])
      }
      # filling pooled array
      p.lat.ind <- 0
      for(k in lat_parts){
        p.lon.ind <- 0
        p.lat.ind <- p.lat.ind + 1
        for(j in lon_parts){
          p.lon.ind <- p.lon.ind + 1
          p.ilons <- match(j,colnames(sliced_data))
          p.ilats <- match(k, rownames(sliced_data))
          pre_pool <- sliced_data[p.ilats,p.ilons,]
          pooled[p.lat.ind, p.lon.ind,] <- apply(pre_pool, 3, pool_fun)
        }
      }
      
      # assign correct dimnames to pooled array
      
      lat_parts_name <- c()
      for(lp in lat_parts){
        name <- paste(c(lp[1], lp[length(lp)]), collapse = "-")
        lat_parts_name <- c(lat_parts_name, name)
      }
      
      lon_parts_name <- c()
      for(lp in lon_parts){
        name <- paste(c(lp[1], lp[length(lp)]), collapse = "-")
        lon_parts_name <- c(lon_parts_name, name)
      }
      dimnames(pooled) = list(lat = lat_parts_name,
                              lon = lon_parts_name,
                              time = as.character(dimnames(sliced_data)$time))
      
      attr(pooled, "name") <- paste0(i$var, "_pooled_arr")
      attr(pooled, "pooled") <- "pooled"
      attr(pooled, "lag_pooled") <- FALSE
      variable_list[["var_pooled_arr"]] <- pooled
      print(paste(i$var, ": pooling done"))
    } #end pooling operation
    
    # start lagging operations
    if(lag==TRUE){
      if(is.unsorted(lag)){
        lags <- sort(lags)
      }
      
      max_lag <- lags[length(lags)]
      lagged_arrays <- list()
      counter =0
      names <- names(variable_list)
      for(k in variable_list){
        counter = counter +1
        for(l in lags){
          lag_array <- k[,,1:(len-l)]
          dimnames(lag_array) <- list(lat= rownames(k), 
                                      lon = colnames(k), 
                                      time=dimnames(sliced_data)$
                                        time[(1+l):len])
          attr(lag_array, "name") <- paste0(names[counter],"_lag_", l)
          attr(lag_array, "lag") <- l
          attr(lag_array, "pooled") <- attributes(k)$pooled
          attr(lag_array, "lag_pooled") <- FALSE
          lagged_arrays[[paste0(names[counter],"_lag_", l)]] <- lag_array
        }
      }
      
      if(lag_pool==TRUE){
        unpooled_arrays <- list()
        pooled_arrays <- list()
        
        
        for(l in lagged_arrays){
          if(attributes(l)$pooled=="pooled"){
            pooled_arrays[[attributes(l)$name]] <- l
          }else{
            unpooled_arrays[[attributes(l)$name]] <- l
          }
        }
        
        parts_unpooled <- split_vec(unpooled_arrays, round((max(lags)-min(lags)+1)/lag_pool_span))
        parts_pooled <- split_vec(pooled_arrays, round((max(lags)-min(lags)+1)/lag_pool_span))
        
        unpooled_lagpooled_arrays <- list()
        for(j in parts_unpooled){
          first_lag <- attributes(j[[1]])$lag
          last_lag <- attributes(j[[length(j)]])$lag
          summedarray_unpooled <- array(0, dim=c(length(dimnames(j[[1]])$lat), length(dimnames(j[[1]])$lon), (len-max(lags))))
          counter <- 0
          for(k in j){
            #print(str(k[,,(length(dimnames(k)$time)-(len-max(lags))+1):length(dimnames(k)$time)]))
            summedarray_unpooled <- summedarray_unpooled + k[,,(length(dimnames(k)$time)-(len-max(lags))+1):length(dimnames(k)$time)]
            counter <- counter + 1
          }
          avg_array <- summedarray_unpooled/counter
          attr(avg_array, "name") <- paste0("var_unpooled_arr_lag_pooled_", first_lag, "_", last_lag)
          attr(avg_array, "first_lag") <- first_lag
          attr(avg_array, "last_lag") <- last_lag
          attr(avg_array, "lag") <- paste0(first_lag, "_", last_lag)
          attr(avg_array, "pooled") <- "unpooled"
          attr(avg_array, "lag_pooled") <- TRUE
          lagged_arrays[[paste0("var_unpooled_arr_lag_pooled_", first_lag, "_", last_lag)]] <- avg_array
          
        }
        
        pooled_lagpooled_arrays <- list()
        for(j in parts_pooled){
          first_lag <- attributes(j[[1]])$lag
          last_lag <- attributes(j[[length(j)]])$lag
          summedarray_pooled <- array(0, dim=c(length(dimnames(j[[1]])$lat), length(dimnames(j[[1]])$lon), (len-max(lags))))
          counter <- 0
          for(k in j){
            #print(str(k[,,(length(dimnames(k)$time)-(len-max(lags))+1):length(dimnames(k)$time)]))
            summedarray_pooled <- summedarray_pooled + k[,,(length(dimnames(k)$time)-(len-max(lags))+1):length(dimnames(k)$time)]
            counter <- counter + 1
          }
          avg_array <- summedarray_pooled/counter
          attr(avg_array, "name") <- paste0("var_pooled_arr_lag_pooled_", first_lag, "_", last_lag)
          attr(avg_array, "first_lag") <- first_lag
          attr(avg_array, "last_lag") <- last_lag
          attr(avg_array, "lag") <- paste0(first_lag, "_", last_lag)
          attr(avg_array, "pooled") <- "pooled"
          attr(avg_array, "lag_pooled") <- TRUE
          lagged_arrays[[paste0("var_pooled_arr_lag_pooled_", first_lag, "_", last_lag)]] <- avg_array
          
        }
      }
      
      
      variable_list[["var_lagged_arrays"]] <- lagged_arrays
      print(paste(i$var, ": lagging done"))
    } #end lagging operations

    # turn array into matrices
    if(make_matrix == TRUE){
      matrices <- list()
      names <- names(variable_list)
      counter <- 0
      for(n in variable_list){
        counter <- counter + 1
        if(names[counter] != "var_lagged_arrays"){
          if(attributes(n)$pooled=="pooled"){
            matrix_var <- matrix(NA_real_, nrow=len, ncol= (length(dimnames(n)$lat)*length(dimnames(n)$lon)))
            col_counter <- 1
            col_dim <- c()
            for(a in dimnames(n)$lat){
              for(b in dimnames(n)$lon){
                matrix_var[,col_counter]=unname(n[a,b,])
                col_dim[col_counter] <- paste0(strtrim(i$var,3), ":", a,"/",b)
                col_counter <- col_counter + 1
              }
            }
            matrix_var <- matrix_var[(nrow(matrix_var)-(len-max_lag)+1):nrow(matrix_var),]
            dimnames(matrix_var) <- list(time=time_dim[(max_lag+1):len], variables=col_dim)
            attr(matrix_var, "pooled") <- attributes(n)$pooled
            attr(matrix_var, "name") <- paste0(paste0(names[counter],"_matrix"))
            attr(matrix_var, "lag") <- attributes(n)$lag
            attr(matrix_var, "lag_pooled") <- attributes(n)$lag_pooled
            matrices[[paste0(names[counter],"_matrix")]] <- matrix_var
          }
        }else{
          names_lag <- names(n)
          #counter_lag <- 0
          for(m in n){
            if(attributes(m)$pooled=="pooled"){
              #counter_lag <- counter_lag + 1
              matrix_var <- matrix(NA_real_, nrow=length(dimnames(m)$time), ncol= (length(dimnames(m)$lat)*length(dimnames(m)$lon)))
              col_counter <- 1
              col_dim <- c()
              for(a in dimnames(m)$lat){
                for(b in dimnames(m)$lon){
                  matrix_var[, col_counter]=unname(m[a,b,])
                  col_dim[col_counter] <-paste0(strtrim(i$var,3), "lag ", attributes(m)$lag,": ", a,"/",b)
                  col_counter <- col_counter + 1
                }
              }
              dimnames(matrix_var) <- list(time=dimnames(m)$time, variables=col_dim)
              attr(matrix_var, "name") <- paste0(attributes(m)$pooled,"_matrix")
              attr(matrix_var, "pooled") <- attributes(m)$pooled
              attr(matrix_var, "lag") <- attributes(m)$lag
              attr(matrix_var, "lag_pooled") <- attributes(m)$lag_pooled
              matrices[[paste0("var", attributes(m)$pooled,"_lag_",attributes(m)$lag ,"_matrix")]] <- matrix_var
            }
          }
        }
      }
      
      # make master matrices one pooled (unpooled not done anymore)
      variable_master_matrix_pooled <- matrix(nrow=(len-max_lag)) #NA_real_, nrow=(len-max_lag)
      col_dim_pooled <- c()
      for(mat in matrices){
        variable_master_matrix_pooled <- cbind2(variable_master_matrix_pooled, mat[(nrow(mat)-(len-max_lag)+1):nrow(mat),])
        col_dim_pooled <- c(col_dim_pooled, dimnames(mat)$variables)
      }
      variable_master_matrix_pooled <- variable_master_matrix_pooled[,-1]
      dimnames(variable_master_matrix_pooled) <- list(time=time_dim[(max_lag+1):len], variables=col_dim_pooled)
      attr(variable_master_matrix_pooled, "lag_pooled") <- "mixed"
      matrices[["var_master_matrix_pooled"]] <- variable_master_matrix_pooled
      variable_list[["var_matrices"]] <- matrices
      ###### variable_master_matrix_lagpooled_pooled
      variable_master_matrix_lagpooled_pooled <- matrix(nrow=(len-max_lag)) #NA_real_, nrow=(len-max_lag)
      col_dim_lagpooled <- c()
      for(mat in matrices){
        if(attributes(mat)$lag_pooled==TRUE){
          variable_master_matrix_lagpooled_pooled <- cbind2(variable_master_matrix_lagpooled_pooled, mat[(nrow(mat)-(len-max_lag)+1):nrow(mat),])
          col_dim_lagpooled <- c(col_dim_lagpooled, dimnames(mat)$variables)
        }
      }
      variable_master_matrix_lagpooled_pooled <- variable_master_matrix_lagpooled_pooled[,-1]
      dimnames(variable_master_matrix_lagpooled_pooled) <- list(time=time_dim[(max_lag+1):len], variables=col_dim_lagpooled)
      attr(variable_master_matrix_lagpooled_pooled, "lag_pooled") <- TRUE
      matrices[["var_master_matrix_lagpooled_pooled"]] <- variable_master_matrix_lagpooled_pooled
      variable_list[["var_matrices"]] <- matrices
      
      print(paste0(i$var, ": making matrices done"))
    }
    
    final_list[[i$var]] <- variable_list
  }

  if(make_matrix==TRUE){
    master_matrix_pooled <-  matrix(NA_real_, nrow=(len-max_lag))
    col_dim_master_pooled <- c()
    for(f in final_list){
      master_matrix_pooled <- cbind2(master_matrix_pooled, f$var_matrices$var_master_matrix_pooled)
      col_dim_master_pooled <- c(col_dim_master_pooled, dimnames(f$var_matrices$var_master_matrix_pooled)$variables)
    }
    master_matrix_pooled <- master_matrix_pooled[,-1]
    dimnames(master_matrix_pooled) <- list(time=time_dim[(max_lag+1):len], variables=col_dim_master_pooled)
    
    ####### master_matrix_lagpooled_pooled
    
    master_matrix_lagpooled_pooled <-  matrix(NA_real_, nrow=(len-max_lag))
    col_dim_master_lagpooled <- c()
    for(f in final_list){
      master_matrix_lagpooled_pooled <- cbind2(master_matrix_lagpooled_pooled, f$var_matrices$var_pooled_arr_matrix)
      col_dim_master_lagpooled <- c(col_dim_master_lagpooled, dimnames(f$var_matrices$var_pooled_arr_matrix)$variables)
      master_matrix_lagpooled_pooled <- cbind2(master_matrix_lagpooled_pooled, f$var_matrices$var_master_matrix_lagpooled_pooled)
      col_dim_master_lagpooled <- c(col_dim_master_lagpooled, dimnames(f$var_matrices$var_master_matrix_lagpooled_pooled)$variables)
    }
    master_matrix_lagpooled_pooled <-master_matrix_lagpooled_pooled[,-1]
    dimnames(master_matrix_lagpooled_pooled) <- list(time=time_dim[(max_lag+1):len], variables=col_dim_master_lagpooled)
    
    final_list[["master_matrix_pooled"]] <- master_matrix_pooled
    final_list[["master_matrix_lagpooled_pooled"]] <- master_matrix_lagpooled_pooled
    
  }
  
  #return the final list
  return(final_list)
}

# Adjust discharge time series to covariates
adjust_discharge_to_meteo <- function(riv) {
  
  # riv is a list of a dataframe consisting of two columns (time, discharge), a longtitude value and a latitude value
  
  start <- 1
  end <- nrow(riv$d)
  if (as.Date(riv$d$time[1]) < as.Date("1979-01-01")) {
    start <- which(riv$d$time=="1979-01-01")
  }
  if (as.Date(riv$d$time[length(riv$d$time)]) > as.Date("2018-12-31"))
    end <- which(riv$d$time=="2018-12-31")
  
  riv$d[start:end,]
}


# Produce stacked time series plots of discharge and covariates
compare_discharge <- function(river, met_var, name=met_var$var, var_unit, log_scale=FALSE, scale_down_by=0.05, transp=0.05) {
  
  # river: output of getRiver function
  # met_var: output of get_ts_ind function
  # log_scale: indicator whether the met_var variable should be modelled on the log scale
  # scale_down_by: scales down the meteo variable time series ( > 0 and might need to be high if log_scale = TRUE)
  # transp: regulates the transparency of the curve of the meteo variables (needs to be between 0-1)
  
  # 1. Subset the discharge such that it is in the range of the meteo variable
  riv_rel <- adjust_discharge_to_meteo(river)
  
  # 2. Find the longest coherent discharge time series in that interval
  riv_rel <- longest_coh_ts(riv_rel)
  
  # 3. Adjust the weather variable to that interval
  start_index <- match(as.character(riv_rel$time[1]), dimnames(met_var$pooled_array)$time)
  end_index <- match(as.character(riv_rel$time[length(riv_rel$time)]), dimnames(met_var$pooled_array)$time)
  var <- met_var$pooled_array[,,start_index:end_index]
  
  if (log_scale) var <- var + 0.5*min(var[var != 0])
  
  # 4. Plot
  
  lats <- dim(var)[1]
  lons <- dim(var)[2]
  ylim_min <- max(min(var)-mean(var)*scale_down_by, 0.00001)
  ylim_max <- max(var)+mean(var)*scale_down_by
  
  plott(riv_rel, reset_par=FALSE, main=paste("Discharge vs.", tolower(name), "across a", lats, "x", lons,  "grid"))
  par(new=T)
  
  if (log_scale){
    plot(riv_rel$time, var[1,1,], log="y", ylim=c(ylim_min, ylim_max), type = "n", axes = FALSE, bty = "n", xlab = "", ylab = "", col=rgb(0.9, 0, 0, alpha=0.2))
  }
  else {
    plot(riv_rel$time, var[1,1,], ylim=c(ylim_min, ylim_max), type = "n", axes = FALSE, bty = "n", xlab = "", ylab = "", col=rgb(0.9, 0, 0, alpha=0.2))
    
  }
  
  mtext(paste(tolower(name), "(", var_unit, ")", sep=""),side=4,col=rgb(0.9, 0, 0, alpha=0.5),line=3)
  eaxis(side=4, 
        col.ticks=rgb(0.9, 0, 0, alpha=0.5),
        col=rgb(0.9, 0, 0, alpha=0.9), 
        col.axis=rgb(0.9, 0, 0, alpha=0.5),
        at.small=F) #rgb(0.9, 0, 0, alpha=0.9)
  
  for (i in 1:lats) {
    for (j in 1:lons) {
      lines(riv_rel$time, var[i,j,], type = "l", col=rgb(0.9, 0, 0, alpha=transp)) #, axes = FALSE, bty = "n", xlab = "", ylab = "", 
    }
  }
}

# Run linear regression on with the given input
river_reg_all <- function(meteo_data, dis_data, 
                          pooled_lags = TRUE, max_lag_used){
  
  # Input: Output of get_ts_all() and getRiver() or longest_coh_ts() -> size has to match
  # Output_ plot of discharge values with fitted time series regression values
  
  
  # #get the design matrix of the regression
  # if(pooled_lags==TRUE){
  #     design_mat <- meteo_data$master_matrix_lagpooled_pooled
  # }else{
  #     design_mat <- meteo_data$master_matrix_pooled
  # }
  
  design_mat <- meteo_data
  # get the response of the regression
  response <- log10(dis_data$discharge)[(1+max_lag_used):length(dis_data$discharge)]
  all_data <- as.data.frame(cbind(response, design_mat))
  # fit the linear model
  fit <- lm(response~., data=all_data)
  
  # return the regression results and the summary of it
  return(list(reg_fit=fit, dis_data = dis_data, max_lag=max_lag_used))
}


# Plot the fitted values of the linear regression model
plot_reg <- function(reg){
  
  # Takes as input the output of river_reg_all and version (might be shortened) of get_river output
  
  plott(reg$dis_data, main = "Log-Discharge with fitted Time Series Regression",
        sub="All variables included")
  lines(reg$dis_data$time[(1+reg$max_lag):length(reg$dis_data$time)], 10^(fitted(reg$reg_fit)), col="#DF536B", cex=0.5)
  legend(reg$dis_data$discharge[10], 100, 
         c("log-discharge", "fitted"), 
         col = c("black", "#DF536B"),
         lty = c(1, 1), cex=0.8, bg = "gray90")
}


# Plot the residual plots and other relevant time series plots
plot_res_reg <- function(reg, add_ts_diagnostics=TRUE, only_ts_diagnostics=FALSE){
  par(mar=c(4,4,4,2))
  if (!only_ts_diagnostics) {
    TA_plot <- plot(reg$reg_fit,1, id.n=0)
    NormPlot <- plot(reg$reg_fit,2, id.n=0)
    ScaleLoc <- plot(reg$reg_fit,3, id.n=0)
    ResLevPlot <- plot(reg$reg_fit,5)
  }
  if (add_ts_diagnostics | only_ts_diagnostics) {
    ACF <- acf(ts(reg$reg_fit$residuals), main="ACF Residuals of fit", lag.max = reg$max_lag)
    PACF <- pacf(ts(reg$reg_fit$residuals), main="PACF Residuals of fit", lag.max = reg$max_lag)
    # ResTS <- plot(reg$dis_data$time[(reg$max_lag+1):length(reg$dis_data$time)], resid(reg$reg_fit), type="l",
    #              main = "Residuals vs time", xlab = "time", ylab = "residual")
    
    # Uncomment last lines, if residual vs time plot needed
  }
}


# Arima Regression
arima_reg <- function(dis_data, meteo_data, order, max_lag_used) {
  
  # Extract design matrix
  design_mat <- meteo_data
  
  ts <- ts(log10(dis_data$discharge)[(1+max_lag_used):length(dis_data$discharge)], start=dis_data$time[1])
  
  fit <- arima(ts,
               order = order,
               xreg=design_mat)
  
  return(list(reg_fit=fit, dis_data=dis_data, var_names=colnames(design_mat), max_lag=max_lag_used, design_mat=design_mat, model_order=order))
}

# Plot Arima fit
plot_arima <- function(arima,sub="All variables included"){
  plott(arima$dis_data, main = "Log-Discharge with fitted Time Series Regression",
        sub=sub)
  lines(arima$dis_data$time[(1+arima$max_lag):length(arima$dis_data$time)], 
        (arima$dis_data$discharge[(1+arima$max_lag):length(arima$dis_data$time)]-arima$reg_fit$residuals), col="#DF536B", cex=0.5)
  legend(arima$dis_data$discharge[10], 100, 
         c("log-discharge", "fitted"), 
         col = c("black", "#DF536B"),
         lty = c(1, 1), cex=0.8, bg = "gray90")
}

arima_pred <- function(model, n.ahead, new_meteo){
  
  design_mat <- model$design_mat
  order <- model$model_order
  predictions <- predict(model$reg_fit, n.ahead=n.ahead, newxreg=new_meteo, se.fit=TRUE)
  return(predictions)
}

# Plot prediction of Arima model
plot_pred <- function(pred, model, discharge, order){
  
  # construction of prediction intervals
  upper_limit <- unname(pred$pred)+1.96*pred$se
  lower_limit <- unname(pred$pred)-1.96*pred$se
  
  # define prediction range
  pred_range <- seq(from=as.Date(attributes(pred$pred)$name[1]), 
                    to=as.Date(attributes(pred$pred)$name[length(attributes(pred$pred)$name)]), 
                    by="days")
  
  plott(discharge, 
        col=alpha(rgb(0,0,0), 0.7), 
        main = paste0("Predictive performance. ARIMA(", paste(order, collapse = ","),")"))
  # 2. plot the fitted model
  lines(model$dis_data$time[(1+model$max_lag):length(model$dis_data$time)], 
        (model$dis_data$discharge[(1+model$max_lag):length(model$dis_data$time)]/10^(model$reg_fit$residuals)), col="#DF536B", cex=0.5)
  # 3. plot the predicted value (year forward)
  lines(pred_range, 10^(pred$pred), col="orange", lwd=2.5)
  # 4. plot the prediction interval
  lines(pred_range, 10^(lower_limit), col=alpha("#2297E6", 0.4), lwd=2)  
  lines(pred_range, 10^(upper_limit), col=alpha("#2297E6", 0.4), lwd=2) 
  polygon(c(pred_range, rev(pred_range)), c(10^(upper_limit), rev(10^(lower_limit))), col=rgb(0, 0, 1,0.2), border=NA)
  # 5. plot legend
  legend("bottomleft", 10, 
         c("log-discharge", "fitted", "prediction", "95% prediction int"), 
         col = c("black", "#DF536B", "orange", "#2297E6" ),
         lty = c(1, 1), cex=0.8, bg = "gray90")
}


# Forward selection with AIC for the arima function
step.arima.forw <- function(ts, order, xreg, verbose=FALSE, numCores=1, tol=0.01){
  
  # input: time series, data matrix (can be taken from get_ts_all)
  
  
  # define all variable names which are eligible to be added
  var.names <- dimnames(xreg)$variables
  # fit the null-model
  fit_null <-  arima(ts, order=order)
  # initiate the round counter
  round_counter <- 0
  # initiate the AIC to beat as the AIC of the null-model
  aic <- fit_null$aic
  imp <- 1
  # initiate the matrix where the, in turn, the best variables get stored
  best_var <- matrix(1, nrow=nrow(xreg))
  chosen_var <-  c()
  # initiate while loop which controls that the AIC get continously improved
  while(imp > 0 && round_counter<ncol(xreg)){
    round_counter = round_counter + 1
    best_var_round <-  c()
    imp <- 0
    best_round_aic <- aic
    var_counter <- 0
    registerDoParallel(numCores)
    aics <- foreach(i = var.names, .combine=rbind) %dopar% {
      variables <- best_var
      i.col <- match(i, dimnames(xreg)$variables)
      fit <- arima(ts, order=order, xreg=cbind2(variables, xreg[,i.col])[,-1])
      c(i, fit$aic)
    }
    stopImplicitCluster()
    min <- aics[which.min(aics[,2]),]
    if((as.numeric(min[2])+tol)<best_round_aic){
      imp <- imp + 1
      improvement <- best_round_aic-as.numeric(min[2])
      best_round_aic <- as.numeric(min[2])
      best_var_round <- min[1]
      chosen_var <- c(chosen_var, best_var_round)
      best_var <- cbind2(best_var, xreg[, match(best_var_round, dimnames(xreg)$variables)])
      var.names <- var.names[var.names != best_var_round]
      aic <- best_round_aic
      print(paste0("Variable added: ", best_var_round, ". Improved AIC to:", round(aic,2), " by ", round(improvement,2), "."))
      print(paste0("Round ", round_counter, " finished. Variables tested:", nrow(aics)))
      print(paste0("Current model size: ", length(chosen_var)))
      print("-------------------------------")
    }else{
      best_var <- best_var[,-1]
      dimnames(best_var) <- list(time=dimnames(xreg)$time , variables = chosen_var)
      print("AIC no longer improved or all variables added")
    }
  }
  return(best_var)
}

# Backward selection with AIC for the arima function
step.arima.backw <- function(ts, order, xreg, verbose=FALSE, numCores=1, tol=0.01){
  
  # define all variable names which are eligible to be added
  var.names <- dimnames(xreg)$variables
  # fit the full-model
  fit_full <-  arima(ts, order=order, xreg=xreg)
  #initiate null-model
  fit_null <-  arima(ts, order=order) 
  # initiate the round counter
  round_counter <- 0
  # initiate the AIC to beat as the AIC of the null-model
  aic <- fit_full$aic
  imp <- 1
  # initiate the matrix where the, in turn, the best variables get stored
  best_var <- xreg
  dropped_var <-  c()
  # initiate while loop which controls that the AIC get continously improved
  while(imp > 0 && round_counter<ncol(xreg)){
    round_counter = round_counter + 1
    best_var_round <-  c()
    imp <- 0
    best_round_aic <- aic
    if(round_counter == ncol(xreg)){
      break
    }
    registerDoParallel(numCores)
    aics <- foreach(i = var.names, .combine=rbind) %dopar% {
      variables <- best_var
      i.col <- match(i, var.names)
      fit <- arima(ts, order=order, xreg=variables[,-i.col])
      c(i, fit$aic)
    }
    stopImplicitCluster()
    min <- aics[which.min(aics[,2]),]
    if((as.numeric(min[2])+tol)<best_round_aic){
      imp <- imp + 1
      improvement <- best_round_aic-as.numeric(min[2])
      best_round_aic <- as.numeric(min[2])
      best_var_round <- min[1]
      dropped_var <- c(dropped_var, best_var_round)
      best_var <- best_var[, -match(best_var_round, var.names)]
      var.names <- var.names[var.names != best_var_round]
      aic <- best_round_aic
      print(paste0("Variable dropped: ", best_var_round, ". Improved AIC to:", round(aic,2), " by ", round(improvement,2), "."))
      print(paste0("Round ", round_counter, " finished. Variables tested:", nrow(aics)))
      print(paste0("Current model size: ", (ncol(xreg)-round_counter)))
      print("-------------------------------")
    }else{
      best_var <- best_var[,-1]
      dimnames(best_var) <- list(time=dimnames(xreg)$time , variables = var.names)
      print("AIC no longer improved")
    }
  }
  if(round_counter == ncol(xreg)){
    if(fit_null$aic<best_round_aic){
      print("All variables dropped")
    }else{
      print("All but one variable dropped")
      return(best_var)
    }
  }
  return(best_var)
}

convert_discharge_to_power <- function(discharge) {
  
  # converts a given discharge value into produced electricity (in MWh) for the provided transformation function
  
  cap <- numeric(length(discharge))
  for (i in 1:length(discharge)) {
    if (discharge[i] > 10) cap[i] <- 3.2
    if (discharge[i] < 0.1) cap[i] <- 0
    if (discharge[i] <= 10 & discharge[i] >= 0.1) cap[i] <- 0.3*discharge[i] +0.2
  }
  
  turb_eff <- exp(-(-4.570433 + 73.54497*(cap/3*100)^(-2.4)))
  gen_eff <- exp(-(-4.582434 + 12.92002*(cap/3*100)^(-1.8)))
  inv_eff <- exp(-(-4.582434 + 12.92002*(cap/3*100)^(-1.8)))
  func <- cap * turb_eff / 100 * gen_eff / 100 * inv_eff / 100 * 24
  return(func)
}

plott_dur <- function(pred, type="l", reset_par=TRUE, ylab = expression("discharge" * (m^3/sec)), main = NULL, sub=NULL, 
                      mar = c(4,4, 1+3*hasTit, 5) + 0.1 + marP, marP = rep(0,4), mgp = c(2, 0.6, 0), col="orange",
                      col.grid  = "lightgray", lty.grid = "dotted", ...)
{
  hasTit <- !is.null(main) && is.character(main) && nzchar(main)
  opar <- par(mar=mar, mgp=mgp)
  if (reset_par) on.exit(par(opar))
  df <- as.data.frame(cbind(1:length(pred$pred), sort(10^(as.vector(pred$pred)), decreasing = TRUE)))
  plot(df, type=type, main=main, sub=sub, col=col, ylab=ylab, xlab="days", ylim=c(0.1, 100), lwd=2, ...) # , yaxt = "n"
  
  abline(v=seq(from=1, to=length(pred$pred), by= 20), col=col.grid, lty=lty.grid) # vertical gridlines
  abline(h=10, col="purple", lty="dashed") #upper cutoff for electricity production
  abline(h=0.1, col="purple", lty="dashed") # lower cutoff for electricity production
  
  # PI interval
  upper_limit <- unname(pred$pred)+1.96*pred$se
  lower_limit <- unname(pred$pred)-1.96*pred$se
  lines(1:length(pred$pred), sort(10^(upper_limit), decreasing = TRUE), col="#2297E6")
  lines(1:length(pred$pred), sort(10^(lower_limit), decreasing = TRUE), col="#2297E6")
  polygon(c(1:length(pred$pred), rev(1:length(pred$pred))),
          c(sort(10^(upper_limit), decreasing = TRUE), rev(sort(10^(lower_limit), decreasing = TRUE))), 
          col=rgb(0, 0, 1,0.2), border=NA)
  legend("topright", 100, 
         c("prediction", "95% prediction int", "true", "cutoffs for elec. prod."), 
         col = c(col, "#2297E6", "black", "purple"),
         lty = c(1, 1, 1, 2), cex=0.8, bg = "gray90")
}

plott_power <- function(pred, type="l", reset_par=TRUE, ylab = expression("electricity(in MWh)"), main = NULL, sub=NULL, 
                        mar = c(4,4, 1+3*hasTit, 5) + 0.1 + marP, marP = rep(0,4), mgp = c(2, 0.6, 0), col="orange",
                        col.grid  = "lightgray", lty.grid = "dotted", ...)
{
  hasTit <- !is.null(main) && is.character(main) && nzchar(main)
  opar <- par(mar=mar, mgp=mgp)
  if (reset_par) on.exit(par(opar))
  df <- as.data.frame(cbind(1:length(pred$pred), convert_discharge_to_power(sort(10^(as.vector(pred$pred)), decreasing = TRUE))))
  plot(df, type=type, main=main, sub=sub, col=col, ylab=ylab, xlab="days", ylim=c(1, 100), lwd=2, ...)
  
  abline(v = seq(from=1, to=length(pred$pred), by= 20), col=col.grid, lty=lty.grid) # vertical gridlines
  
  # PI interval
  upper_limit <- unname(pred$pred)+1.96*pred$se
  lower_limit <- unname(pred$pred)-1.96*pred$se
  lines(1:length(pred$pred), convert_discharge_to_power(sort(10^(upper_limit), decreasing = TRUE)), col="green")
  lines(1:length(pred$pred), convert_discharge_to_power(sort(10^(lower_limit), decreasing = TRUE)), col="green")
  polygon(c(1:length(pred$pred), rev(1:length(pred$pred))),
          c(convert_discharge_to_power(sort(10^(upper_limit), decreasing = TRUE)), rev(convert_discharge_to_power(sort(10^(lower_limit), decreasing = TRUE)))), 
          col=rgb(0, 1, 0,0.2), border=NA)
  
  legend("bottomleft", 100, 
         c("prediction", "95% prediction int", "true"), 
         col = c(col, "green", "black"),
         lty = c(1, 1), cex=0.8, bg = "gray90")
}


#Helper functions

split_vec <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 


hitReturn <- function(message)
{
  invisible(readline(sprintf("%s -- hit return", message)));
}

plotAllColumns <- function(dframe, waitFunc = hitReturn)
{
  for (n in colnames(dframe))
  {
    plot(dframe[[n]]);
    waitFunc(n);
  }
}
