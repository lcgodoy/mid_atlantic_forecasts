functions {
  // Gaussian temperature dependence function
  real T_dep(real sbt, real Topt, real width, real exp_yn){
    if(exp_yn==1){
      return exp(-0.5 * ((sbt - Topt)/width)^2); // gaussian temperature-dependent function; used for T-dependent recruitment and mortality?
    } else{
      // return (-0.1 * ((sbt - Topt)/width)^2); // gaussian temperature-dependent function, not exponentiated, for temperature-dependent process that are exponentiated later (like movement)
      // I'm taking this out because I can't get it to be non-negative without the exp! 
      
      // return normal_lpdf(sbt | Topt, width);
      return log((1 / sqrt(2 * pi() * width)) * exp(-pow(sbt - Topt,2) / (2 * pow(width,2))));
    }
  }
  
  matrix age_at_length_key(real loo, real l0, real k, real cv, int n_lbins, int n_ages){
    
    vector[n_lbins] mean_age_at_length;
    vector[n_lbins] sigma_age_at_length;
    matrix[n_lbins, n_ages] prob_age_at_length;
    
    // make vector of mean ages for each length 
    for(i in 1:n_lbins){
      mean_age_at_length[i] = log((loo - i) / (loo - l0)) / -k; 
      sigma_age_at_length[i] = mean_age_at_length[i] * cv; 
    }
    
    for(i in 1:n_lbins){
      for(j in 1:n_ages){
        if(j < n_ages){
          prob_age_at_length[i,j] = normal_cdf(j+1, mean_age_at_length[i], sigma_age_at_length[i]) - normal_cdf(j, mean_age_at_length[i], sigma_age_at_length[i]);  // analog of pnorm in R
        }
        else{
          prob_age_at_length[i,j] = normal_cdf(j, mean_age_at_length[i], sigma_age_at_length[i]);
        } // close if/else
      } // close ages
    } // close lengths
    return prob_age_at_length;
  } // close function
  
  vector colSums(matrix M){
    int ncol; 
    vector[cols(M)] sums; 
    
    ncol = cols(M); 
    
    for(i in 1:ncol){
      sums[i] = sum(M[,i]); //sums[i] = sum(col(M,i)); 
    }
    return(sums);
  }
  
  
  // cumulative sum function
  // slightly different than the built-in Stan one -- this takes a 1D array of numbers and returns a single real number that sums elements 1:x of that array
  real csum(real[] to_sum, int x){
    return(sum(to_sum[1:x])); 
  }
  
  // functions to calculate range quantiles
  // https://discourse.mc-stan.org/t/how-to-find-the-location-of-a-value-in-a-vector/19768/2
  
  // first: a function that calculates how many matches there are in vector x to real y (need this for counting along the abundance vector to estimate quantiles)
  int num_matches(vector x, real y) {
    int n = 0;
    for (i in 1:rows(x))
    if (x[i] == y)
    n += 1;
    return n;
  }
  
  // next: a function that returns the positions of those matches along vector x 
  int[] which_equal(vector x, real y) {
    int match_positions[num_matches(x, y)];
    int pos = 1;
    for (i in 1:num_elements(x)) { // example code used size(), got error "no matches for available argument signatures for size"; but beware of applying this function to >1D arrays 
    if (x[i] == y) {
      match_positions[pos] = i;
      pos += 1;
    }
    }
    return match_positions;
  }
  
  // finally: function to calculate range quantiles
  
  real calculate_range_quantile(int np, vector patches, real[] dens_by_patch, real quantile_out){
    vector[np] csum_dens;
    vector[np] csum_to_edge; 
    real quant_position; 
    real cutoff; 
    int cutoff_id; 
    
    cutoff = sum(dens_by_patch[1:np]) * quantile_out; // where along the vector of counts does the quantile cutoff fall? (then need to translate that into a patch position, below)  
    
    for(i in 1:np){
      csum_dens[i] = csum(dens_by_patch[1:np], i); // calculate cumulative sum of density along each patch 
      if(csum_dens[i] <= cutoff){
        csum_to_edge[i] = csum_dens[i]; // keep only the patches below the edge 
      } else {
        csum_to_edge[i] = 0; 
      }
    }
    
    cutoff_id = min(which_equal(csum_to_edge, 0)); // get lowest patch with a 0 (that's the patch where the edge falls)
    
    quant_position = ((cutoff - max(csum_to_edge)) / dens_by_patch[cutoff_id]) + cutoff_id - 1; // calculate what proportion of the edge-containing patch is "filled in" by the actual fish up to the weighted quantile (that's the decimal) and then add in the other patches 
    
    return(quant_position); 
    
  } // close quantile function
  
  // FUNCTIONS TO RUN POPULATION DYNAMICS MODEL 
  
  // note that I replaced ny_train with just ny so this can be used for both the model fit and projection 
  
  matrix[] pop_dy(int np, int ny, int n_ages, int n_lbins, vector bin_mids, real[,] sbt, real m, real loo, vector maturity_at_age, vector wt_at_age, real[,] f, matrix adj_m, real h, int age_at_maturity, // data 
  int process_error_toggle, int spawner_recruit_relationship, int T_dep_mortality, int T_dep_movement, int T_dep_recruitment, int exp_yn, // toggles 
  vector rec_dev, real log_r0, real sigma_r_raw, real p_length_50_sel, real sel_delta, real log_mean_recruits, real Topt, real width, real d, real beta_t, vector init_dep, real beta_rec, real[] raw, real alpha // parameters 
  ) {
    //**********************
    matrix[np, n_ages] n_at_age_hat[ny];
    real mean_recruits; 
  //  vector[ny-1] rec_dev;
    //***************************
    real r0; 
    real sigma_r; 
    real ssb0; 
    matrix[np, ny] ssb;
    vector[n_ages] unfished;
    real length_50_sel;
    vector[n_lbins] selectivity_at_bin; 
    matrix[np, ny] T_adjust;
    real surv[np, n_ages, ny];
    matrix[np, np] diff_m; // dispersal matrix (setting up for temperature-dependent movement, but doesn't actually need to vary by year)
    matrix[np, np] mov_inst_m[ny]; // array containing instantaneous movement matrices (dispersal + temperature-dependent taxis)
    matrix[np, np] mov_m[ny]; // Dan thinks this is the wrong syntax -- need to fix  // we sorted this out right? 
    matrix[np, np] T_adjust_m[ny]; // array containing temperature preference matrices
    matrix[np, np] tax_m[ny]; // array containing taxis between patches every year 
    vector[np] v_in; // vector for matrix multiplication
    vector[np] v_out; // vector for matrix multiplication 
    
    r0 = exp(log_r0);
    sigma_r = sigma_r_raw * process_error_toggle;
    //    ssb0 = -999;
    unfished[1] = r0;
    
    if(spawner_recruit_relationship==1){
      for(a in 2:n_ages){
        
        unfished[a] = unfished[a-1] * exp(-m);
      }
      
      ssb0 = sum(unfished .* maturity_at_age .* wt_at_age);
    }
    length_50_sel = loo * p_length_50_sel; // Dan made a note to change this sometime
    selectivity_at_bin = 1.0 ./ (1 + exp(-log(19) * ((bin_mids - length_50_sel) / sel_delta))); // selectivity ogive at age
    mean_recruits = exp(log_mean_recruits);
    
    
    // calculate temperature-dependence correction factor for each patch and year depending on sbt (not used for all process models)
    for(p in 1:np){
      for(y in 1:ny){ 
        T_adjust[p,y] =   T_dep(sbt[p,y], Topt, width, exp_yn);  
      } // close years
    } // close patches
    
    
    // calculate total annual mortality from instantaneous natural + fishing mortality data 
    
    for(p in 1:np){
      for(a in 1:n_ages){
        for(y in 1:ny){
          
          if(T_dep_mortality==1){
            surv[p,a,y] = exp(-((f[a,y] + m) * T_adjust[p,y]));
          }
          
          if(T_dep_mortality==0){
            surv[p,a,y] = exp(-(f[a,y] + m)) ;
          }
          
        }
      }
    }
    
    
    //  calculate annual movement matrix if using spatially varying dispersal
    if(T_dep_movement==1){
      
      diff_m = adj_m * d; // set up matrix as adjacency * dispersal
      diff_m = add_diag(diff_m, -1 * colSums(diff_m)); // rescale the columns so they will sum to 1 at the end
      
      // set up array of matrices containing 
      for(y in 1:ny){
        for(i in 1:np){
          for(j in 1:np){
            // in R this is just outer(np, np, "-") 
            
            if(exp_yn==1){
              T_adjust_m[y,i,j] = exp(beta_t * (log(T_adjust[i,y]) - log(T_adjust[j,y]))); 
              
            } else {
              
              T_adjust_m[y,i,j] = fmin(500,exp(beta_t * (T_adjust[i,y] - T_adjust[j,y]))); 
            }
          }
        }
        tax_m[y] = adj_m .* T_adjust_m[y]; 
        
        tax_m[y] = add_diag(tax_m[y], -1 * colSums(tax_m[y])); // fill in the diagonal with within-patch "taxis" so everything sums to 1 
        
        
        mov_inst_m[y] =  diff_m + tax_m[y]; // movement as a sum of diffusion and taxis (can cancel each other out)
        mov_m[y] = matrix_exp(mov_inst_m[y]); // matrix exponentiate, although see https://discourse.mc-stan.org/t/matrix-exponential-function/9595
        
        
        if ((sum(colSums(mov_m[y])) / np - 1) > .001 ){
          print("Something has gone very wrong, movement matrix columns do not sum to 1");
          print(colSums(mov_m[y]));
          print("width is", width);
          print("Topt is", Topt);
          print(diagonal(mov_inst_m[y]));
          
        }
      }
      //   } else {
        //     for(z in 1:np){
          //       for(x in 1:np){
            //         diff_m[z,x] = 999; 
            //       }
            //     }
            
    } // close else (end of movement model)
    
    
    // fill in year 1 of n_at_age_hat, initialized with mean_recruits 
    for(p in 1:np){
      for(a in 1:n_ages){
        if(a==1){
          
          if(T_dep_recruitment==1 && spawner_recruit_relationship==0){
            n_at_age_hat[1,p,a] = init_dep[p] * mean_recruits * beta_rec * T_adjust[p,1] * exp(sigma_r * raw[1] - pow(sigma_r,2) / 2); // initialize age 0 with mean recruitment in every patch
          }
          if(T_dep_recruitment==0 && spawner_recruit_relationship==0){
            n_at_age_hat[1,p,a] = init_dep[p] * mean_recruits * exp(sigma_r * raw[1] - pow(sigma_r,2) / 2); // initialize age 0 with mean recruitment in every patch
          }
          if(T_dep_recruitment==0 && spawner_recruit_relationship==1){
            n_at_age_hat[1,p,a] = init_dep[p] * r0 *  exp(sigma_r * raw[1] - pow(sigma_r,2) / 2); // scale it down a bit -- historical fishing was still occurring
          }
          if(T_dep_recruitment==1 && spawner_recruit_relationship==1){
            n_at_age_hat[1,p,a] = init_dep[p] *r0 *  exp(sigma_r *raw[1] - pow(sigma_r,2) / 2) * T_adjust[p,1] * beta_rec;
          }
        } // close age==1 case
        else{
          n_at_age_hat[1,p,a] = n_at_age_hat[1,p,a-1] * surv[p,a-1,1]; // initialize population with mean recruitment propogated through age classes with mortality
        }
        
      } // close ages
      
      ssb[p,1] = sum(to_vector(n_at_age_hat[1,p,1:n_ages]) .* maturity_at_age .* wt_at_age);
      
    } // close patches (end of year 1 code) 
    
    // main pop dy loop for years 2+ 
    
    for (y in 2:ny){
      
      // describe population dynamics
      for(p in 1:np){
        
        // density-independent, temperature-dependent recruitment of age 1
        
        if(T_dep_recruitment==1 && spawner_recruit_relationship==0){
          n_at_age_hat[y,p,1] = mean_recruits * exp(rec_dev[y-1] - pow(sigma_r,2)/2) * T_adjust[p,y-1] * beta_rec;
        }
        if(T_dep_recruitment==0 && spawner_recruit_relationship==0){
          n_at_age_hat[y,p,1] = mean_recruits * exp(rec_dev[y-1] - pow(sigma_r,2)/2) ;
        }
        
        if(T_dep_recruitment==0 && spawner_recruit_relationship==1){
          n_at_age_hat[y,p,1] = (0.8 * r0 * h * ssb[p, y-1]) / (0.2 * ssb0 * (1-h) + ssb[p, y-1] * (h - 0.2));
          
          n_at_age_hat[y,p,1] =  n_at_age_hat[y,p,1] *  exp(rec_dev[y-1] - pow(sigma_r,2)/2);
          
        }
        if(T_dep_recruitment==1 && spawner_recruit_relationship==1){
          n_at_age_hat[y,p,1] = ((0.8 * r0 * h * ssb[p, y-1]) / (0.2 * ssb0 * (1-h) +  ssb[p, y-1] * (h - 0.2))) * T_adjust[p,y-1];
          
          n_at_age_hat[y,p,1] =  n_at_age_hat[y,p,1] *  exp(rec_dev[y-1] - pow(sigma_r,2)/2);
        }
        // 
        // why estimate raw and sigma_r? we want to estimate process error
        // if we got rid of sigma_r, we would be saying that raw could be anything
        // that means raw could pick any value, and it would pick deviates to perfectly match abundance index
        // sigma_r scales the amount of process error that we say is possible
        // letting it go to infinity means the model will always fit the data perfectly 
        // raw is the realized process error 
        // exp(rec_dev[y-1] - pow(sigma_r,2)/2) = random variable with mean 0 and SD sigma_r
        // allows us to have a different recruitment deviation every year even though sigma_r is constant 
        
        // pop dy for non-reproductive ages 
        if(age_at_maturity > 1){ // confirm that there are non-reproductive age classes above 1
        
        n_at_age_hat[y,p,2:(age_at_maturity-1)] = (n_at_age_hat[y-1,p, 1:(age_at_maturity-2)]) .* to_row_vector(surv[p,1:(age_at_maturity-2),y-1]);
        
        } // close if 
      } // close patches 
      
      // pop dy for reproductive adults
      
      
      if(T_dep_movement==0){
        
        // mortality and dispersal are happening simultaneously here, between generations
        // because neither is patch-specific I don't think the order matters
        for(p in 1:np){
          for(a in age_at_maturity:n_ages){
            // edge cases -- edges are reflecting
            if(p==1){
              n_at_age_hat[y,p,a] = n_at_age_hat[y-1,p, a-1] * surv[p,a-1,y-1] * (1-d) + n_at_age_hat[y-1,+1, a-1] * surv[p+1,a-1,y-1] * d;
            } // close patch 1 case 
            
            else if(p==np){
              n_at_age_hat[y,p,a] = n_at_age_hat[y-1,p, a-1] * surv[p,a-1,y-1] * (1-d) + n_at_age_hat[y-1,p-1, a-1] * surv[p-1,a-1,y-1] * d;
            } // close highest patch
            
            else{
              n_at_age_hat[y,p,a] = n_at_age_hat[y-1,p, a-1] * surv[p,a-1,y-1] * (1-2*d) + n_at_age_hat[y-1,p-1, a-1] * surv[p-1,a-1,y-1] * d + n_at_age_hat[y-1,p+1, a-1] * surv[p+1,a-1,y-1] * d;
              
            } // close if/else for all other patches
            
          }// close ages
        } // close patches 
      } // close T-dep movement if 
      
      // this code block calculates adult population size based on survival and directional movement 
      if(T_dep_movement==1){
        
        
        for (p in 1:np){
          
          n_at_age_hat[y,p, age_at_maturity:n_ages] = n_at_age_hat[y-1,p,  (age_at_maturity - 1):(n_ages - 1)] .* to_row_vector(surv[p, (age_at_maturity - 1):(n_ages - 1),y-1]);
          
        }
        
        for(a in age_at_maturity:n_ages){
          
          // some acrobatics required here, because Stan won't do matrix multiplication with an array of reals like n_at_age_hat
          // instead we do the matrix multiplication with a placeholder vector and then populate n_at_age_hat 
          
          // fill in placeholder vector with reproductive ages across patches, and do mortality   
          for(p in 1:np){
            v_in[p] = n_at_age_hat[y,p, a]; 
          }
          v_out = mov_m[y] * v_in; // redistribute each age among patches according to the movement matrix 
          
          // fill in n_at_age_hat
          for(p in 1:np){
            n_at_age_hat[y,p,a] = v_out[p]; 
            
          }
          
        } // close ages
      }// close T-dep movement if 
      
      for (p in 1:np){
        
        ssb[p,y]  =  sum(to_vector(n_at_age_hat[y,p,1:n_ages]) .* maturity_at_age .* wt_at_age);
        
      }
      
      
    } // close year 2+ loop
    
    return n_at_age_hat; 
    
  } // close pop dy function 
  
} // close functions block

data {
  
  // survey data 
  
  int n_ages; // number of ages
  
  int np; // number of patches
  
  vector[np] patches; // easier than creating it in Stan
  
  int ny_train; // years for training
  
  int ny_proj; // number of years to forecast 
  
  int n_lbins; // number of length bins (here just the range of cm values)
  
  matrix[n_ages, n_lbins] l_at_a_key;
  
  vector[n_ages] wt_at_age;
  
  real dens[np, ny_train]; // MEAN density of individuals of any age in each haul; used for rescaling the abundance to fit to our data
  
  // environmental data 
  
  real sbt[np, ny_train]; // temperature data for training
  
  real sbt_proj[np, ny_proj];
  
  // fish data
  
  real m;  // total mortality 
  
  real k;
  
  real loo;
  
  real t0;
  
  real cv;
  
  real f[n_ages, ny_train]; 
  
  real f_proj[n_ages, (ny_proj+1)];
  
  real length_50_sel_guess;
  
  vector<lower=0>[n_lbins] bin_mids;
  
  int sel_100; // age at which selectivity is 1 
  
  int age_at_maturity;
  
  int<lower = 0, upper = 1> do_dirichlet;
  
  int<lower = 0, upper = 1> T_dep_recruitment;
  
  int<lower = 0, upper = 1> T_dep_mortality;
  
  int<lower = 0, upper = 1> T_dep_movement;
  
  int<lower = 0, upper = 1> eval_l_comps;
  
  int<lower = 0, upper = 1> spawner_recruit_relationship;
  
  int<lower = 0, upper = 1> run_forecast;
  
  int<lower = 0, upper = 1> exp_yn;
  
  int<lower = 0, upper = 1> process_error_toggle;
  
  int number_quantiles;
  
  real quantiles_calc[number_quantiles];
  
  int n_at_length[np, n_lbins, ny_train]; // SUM number of individuals in each length bin, patch, and year; used for age composition only, because the magnitude is determined by sampling effort
  
  real sigma_obs_cv; // cv of sigma_obs prior
  
  real h; // steepness
}

transformed data{
  
  vector[n_ages] maturity_at_age; // vector of probabilities of being mature at each age, currently binary (0/1) and taken as known
  
  matrix[np, np] adj_m; // adjacency matrix for patches 
  
  matrix[np, np] outer; // outer difference matrix 
  
  //int exp_yn; 
  
  // in R this is just outer(np, np, "-") 
  for(i in 1:np){
    for(j in 1:np){
      outer[i, j] = i - j; // should just create a matrix that counts from 1 to np along rows and columns
      if(fabs(outer[i,j])==1) {
        adj_m[i,j] = 1; // fill with 1s for adjacent patches and 0s otherwise 
      }else{
        adj_m[i,j] = 0; 
      }
    }
  }
  
  // print("the outer matrix is ",outer); 
  // print("the adjacency matrix is ",adj_m); 
  
  for(a in 1:n_ages){
    if(a < age_at_maturity){
      maturity_at_age[a]=0;
    }else{
      maturity_at_age[a]=1;
    }
  }
  
  // exp_yn=1;
  
  // matrix[n_lbins, n_ages] prob_age_at_length;
  
  // int abund[np, n_ages, ny_train]; 
  
  // vector[n_ages] age_dist[n_lbins]; // create array of vectors 
  
  // prob_age_at_length = age_at_length_key(
    //   loo=loo,
    //   l0=l0,
    //   k=k,
    //   cv=cv,
    //   n_lbins=n_lbins,
    //   n_ages=n_ages ); 
    
    // for(p in 1:np){
      //   for(l in 1:n_lbins){
        //     for(y in 1:ny_train){
          //       age_dist[l] = prob_age_at_length[l,] * n_at_length[p, l, y]; // get vector of probabilities for each length, times counts
          //     }
          //   }
          // }
          
          // for(p in 1:np){
            //   for(a in 1:n_ages){
              //     for(y in 1:ny_train){
                //       abund[p,a,y] = sum(age_dist[,a]); // not sure I'm indexing age_dist correctly, and still need to round to int somehow!
                //     }
                //   }
                // }
}

parameters{
  
  // real<lower = 1e-6> sigma_total;
  
  real<lower = 1e-6> sigma_obs;
  
  real<lower = 1e-6> sigma_r_raw;
  
  
  real<lower=0.5> width; // sensitivity to temperature variation
  
  real Topt; //  temp at which recruitment is maximized
  
  real<lower = 0, upper = 0.99> alpha; // autocorrelation term
  
  real  log_mean_recruits; // log mean recruits per patch, changed to one value for all space/time
  
  vector[ny_train] raw; // array of raw recruitment deviates, changed to one value per year
  
  real<upper = 0.8> p_length_50_sel; // length at 50% selectivity
  
  real<lower = 1e-6> sel_delta;
  
  real<lower=0> beta_obs; // controls how fast detection goes up with abundance
  
  // real<lower=0, upper=0.333> d; // dispersal fraction (0.333 = perfect admixture)
  
  real<lower=0, upper=1> d; // increasing bounds on this for the temperature-dependent movement model 
  
  vector<lower = 0, upper = 1>[np] init_dep;
  
  real<lower = 0, upper = 2> theta_d;
  
  real beta_t; // responsiveness of movement to temperature
  
  real beta_rec; // responsivenses of mean recruits to temperature
  
  real beta_obs_int; // intercept of detection probability
  
  real log_r0;
  
}

transformed parameters{
  
  matrix[np, ny_train] theta; // Bernoulli probability of encounter  
  
  matrix[np, n_ages] n_at_age_hat[ny_train];
  
  matrix[np, n_lbins] n_at_length_hat[ny_train]; // array number of years containing matrices with numbers at patch, length bin, and year 
  
  real density_hat [np, ny_train]; // for tracking sum density 
  
  vector[ny_train-1] rec_dev; // array of realized recruitment deviates, also now only 1/yr (it's a good or bad year everywhere)
  
  vector[n_lbins] selectivity_at_bin; // mean selectivity at length bin midpoint
  
  vector[ny_train] centroid; 
  
  vector[ny_train-1] abund_lr; 
  
  matrix[number_quantiles, ny_train] range_quantiles; 
  
        // calculate recruitment deviates every year (not patch-specific)

  for(y in 2:ny_train){
    
      if (y == 2){ 
        rec_dev[y-1]  =  sigma_r_raw * process_error_toggle * raw[y]; // initialize first year of rec_dev with raw (process error) -- now not patch-specific. sigma_r_raw * process_error_toggle = sigma_r 
      } // close y==2 case  
      else {
        
        rec_dev[y-1] =  alpha * rec_dev[y-2] +  sqrt(1 - pow(alpha,2)) *  sigma_r_raw * process_error_toggle * raw[y]; 
        
      } // close ifelse
      
  }
  
  n_at_age_hat = pop_dy(np=np, ny_train=ny, n_ages=n_ages, n_lbins= n_lbins, bin_mids= bin_mids, sbt= sbt, m= m, loo= loo, maturity_at_age= maturity_at_age, wt_at_age= wt_at_age, f= f, adj_m= adj_m, h= h, age_at_maturity= age_at_maturity,   process_error_toggle= process_error_toggle, spawner_recruit_relationship= spawner_recruit_relationship, T_dep_mortality= T_dep_mortality, T_dep_movement= T_dep_movement, T_dep_recruitment= T_dep_recruitment, exp_yn= exp_yn,   rec_dev= rec_dev, log_r0= log_r0, sigma_r_raw= sigma_r_raw, p_length_50_sel= p_length_50_sel, sel_delta= sel_delta, log_mean_recruits= log_mean_recruits, Topt= Topt, width= width, d= d, beta_t= beta_t, init_dep= init_dep, beta_rec= beta_rec, raw= raw, alpha=alpha); 
  
  for(y in 1:ny_train){
    for(p in 1:np){
      
      n_at_length_hat[y,p,1:n_lbins] = ((l_at_a_key' * to_vector(n_at_age_hat[y,p,1:n_ages])) .* selectivity_at_bin)'; // convert numbers at age to numbers at length. The assignment looks confusing here because this is an array of length y containing a bunch of matrices of dim p and n_lbins
      // see https://mc-stan.org/docs/2_18/reference-manual/array-data-types-section.html
      
      density_hat[p,y] = sum((to_vector(n_at_length_hat[y,p,1:n_lbins])));
      
      theta[p,y] = ((1/(1+exp(-(beta_obs_int + beta_obs*log(density_hat[p,y] + 1e-6))))));
      
    } // close patches
    
    for(q in 1:number_quantiles){
      // calculate every range quantile q for every year y
      range_quantiles[q, y] = calculate_range_quantile(np, patches, density_hat[,y], quantiles_calc[q]);
    }
    
    centroid[y] = sum(to_vector(density_hat[,y]) .* patches) / sum(to_vector(density_hat[,y])); // calculate center of gravity
    
    if(y>1) {
      abund_lr[y-1] = log(sum(to_vector(density_hat[,y])) / sum(to_vector(density_hat[,y-1])));
    }
  }
  
  
} // close transformed parameters block

model {
  
  real n;
  
  vector[n_lbins] prob_hat;
  
  vector[n_lbins] prob;
  
  real dml_tmp;
  
  real test;
  
  init_dep ~ beta(1.5,3);
  
  beta_obs ~ normal(0.001,0.1); 
  
  beta_obs_int ~ normal(-100,4);
  
  raw ~ normal(0, 1);
  
  sigma_r_raw ~ normal(.2,.1);
  
  sigma_obs ~ normal(0.2,sigma_obs_cv);
  
  log_mean_recruits ~ normal(7,5);
  
  log_r0 ~ normal(15,5);
  
  Topt ~ normal(18, 2);
  
  width ~ normal(4, 2); 
  
  beta_t ~ normal(0,2);
  
  beta_rec ~ normal(0,2);
  
  alpha ~  beta(12,20); // concentrated around 0.4
  
  d ~ normal(0.1, 0.1); // diffusion rate as a proportion of total population size within the patch
  
  p_length_50_sel ~ normal(length_50_sel_guess/loo, .2);
  
  sel_delta ~ normal(2,4);
  
  theta_d ~ normal(0.5,.5);
  
  
  for(y in 2:ny_train) {
    
    for(p in 1:np){
      
      if((dens[p,y]) > 0) {
        
        if(eval_l_comps==1){
          if (sum(n_at_length[p,1:n_lbins,y]) > 0) {
            
            if (do_dirichlet == 1){
              
              prob_hat = (to_vector(n_at_length_hat[y,p,1:n_lbins])  / sum(to_vector(n_at_length_hat[y,p,1:n_lbins])));
              
              prob = (to_vector(n_at_length[p,1:n_lbins,y])  / sum(to_vector(n_at_length[p,1:n_lbins,y])));
              
              n = sum(n_at_length[p,1:n_lbins,y]);
              
              dml_tmp = lgamma(n + 1) -  sum(lgamma(n * prob + 1)) + lgamma(theta_d * n) - lgamma(n + theta_d * n) + sum(lgamma(n * prob + theta_d * n * prob_hat) - lgamma(theta_d * n * prob_hat)); // see https://github.com/merrillrudd/LIME/blob/9dcfc7f7d5f56f280767c6900972de94dd1fea3b/src/LIME.cpp#L559 for log transformation of dirichlet-multinomial in Thorston et al. 2017
              
              target += (dml_tmp);
              
            } else {
              
              (n_at_length[p,1:n_lbins,y]) ~ multinomial((to_vector(n_at_length_hat[y,p,1:n_lbins])  / sum(to_vector(n_at_length_hat[y,p,1:n_lbins]))));
              
            } // close dirichlet statement
            
          } // close if any positive length comps
          
        } // close eval_length_comps
        
        log(dens[p,y]) ~ normal(log(density_hat[p,y] + 1e-6), sigma_obs); 
        
        1 ~ bernoulli(theta[p,y]);
        
        
      } else { // only evaluate density if there are length comps to evaluate
      
      0 ~ bernoulli(theta[p,y]);
      
      } // close else 
    } // close patch loop
    
  } // close year loop
  
  
}


generated quantities {
  
  // obs denotes the actual observed data as opposed to the true state
  // anything indexed over ny_proj+1 starts the year before the first forecast (i.e., the last year of the model fit) 
  
  matrix[np, ny_train] dens_pp;
//  matrix[np, n_ages] n_at_age_obs_proj[ny_proj+1]; // how should we add this back in? 
  matrix[np, n_lbins] n_at_length_obs_proj[ny_proj+1];
  matrix[np, n_ages] n_at_age_proj[ny_proj+1];
  matrix[np, n_lbins] n_at_length_proj[ny_proj+1];
  real density_obs_proj[np, ny_proj]; 
  real density_proj[np, ny_proj]; 
  real T_adjust_proj[np, ny_proj];
  vector[ny_proj-1] rec_dev_proj;
  vector[ny_proj] raw_proj;
  real surv_proj[np, n_ages, (ny_proj+1)];
  matrix[np, n_lbins] proj_n_at_length_hat[ny_proj];
  matrix[np, np] T_adjust_m_proj[ny_proj]; 
  matrix[np, np] tax_m_proj[ny_proj]; 
  matrix[np, np] mov_inst_m_proj[ny_proj]; 
  matrix[np, np] mov_m_proj[ny_proj]; 
  matrix[np, ny_proj] ssb_proj;
  vector[np] v_in_proj; // pretty sure we could reuse v_in here but just in case
  vector[np] v_out_proj; 
  vector[ny_proj] centroid_proj; 
  vector[ny_proj-1] abund_lr_proj; 
  matrix[number_quantiles, ny_proj] range_quantiles_proj; 
  matrix[np, ny_proj] theta_proj; // Bernoulli probability of encounter  
  
  
  
  // generate posterior predictive distributions for training data
  
  for (y in 1:ny_train){
    
    for (p in 1:np){
      // ignoring error around length sampling for now
      dens_pp[p,y] = bernoulli_rng(theta[p,y]) * exp(normal_rng(log(density_hat[p,y] + 1e-6), sigma_obs)); 
      
    }
    
  }
  
  if(run_forecast==1){
    
    // initialize with final year(s) of our model
    n_at_age_proj[1,,] = n_at_age_hat[ny_train,,];
    
    raw_proj[1] = raw[ny_train];
    
    // run pop dy 
      n_at_age_proj[2:ny_proj+1,,] = pop_dy(np=np, ny_proj=ny, n_ages=n_ages, n_lbins= n_lbins, bin_mids= bin_mids, sbt_proj= sbt, m= m, loo= loo, maturity_at_age= maturity_at_age, wt_at_age= wt_at_age, f= f, adj_m= adj_m, h= h, age_at_maturity= age_at_maturity,   process_error_toggle= process_error_toggle, spawner_recruit_relationship= spawner_recruit_relationship, T_dep_mortality= T_dep_mortality, T_dep_movement= T_dep_movement, T_dep_recruitment= T_dep_recruitment, exp_yn= exp_yn,   rec_dev= rec_dev, log_r0= log_r0, sigma_r_raw= sigma_r_raw, p_length_50_sel= p_length_50_sel, sel_delta= sel_delta, log_mean_recruits= log_mean_recruits, Topt= Topt, width= width, d= d, beta_t= beta_t, init_dep= init_dep, beta_rec= beta_rec, raw= raw, alpha= alpha ); 
    
    for(p in 1:np){
      for(y in 1:ny_proj){
        T_adjust_proj[p,y] = T_dep(sbt_proj[p,y], Topt, width, exp_yn);
      } // close years
    } // close patches
    
    for(p in 1:np){
      for(a in 1:n_ages){
        for(y in 1:ny_proj+1){
          
          if(T_dep_mortality==1){
            
            // fill in with last year of training data 
            if(y==1){
              surv_proj[p,a,y] = exp(-((f_proj[a,y] + m) * T_adjust[p,ny_train]));
            }
            
            if(y>1){
              surv_proj[p,a,y] = exp(-((f_proj[a,y] + m) * T_adjust_proj[p,y-1]));
            }
          }
          
          if(T_dep_mortality==0){
            surv_proj[p,a,y] = exp(-(f_proj[a,y] + m)) ;
          }
          
        } // close years
      } // close ages
      ssb_proj[p,1] = sum(to_vector(n_at_age_proj[1,p,1:n_ages]) .* maturity_at_age .* wt_at_age); // initialize ssb_proj 
      
    } // close patches
    
    
    
    if(T_dep_movement==1){
      
      for(y in 1:ny_proj){
        for(i in 1:np){
          for(j in 1:np){
            
            if(exp_yn==1){
              T_adjust_m_proj[y,i,j] = exp(beta_t * (log(T_adjust_proj[i,y]) - log(T_adjust_proj[j,y]))); 
              
            } else {
              
              T_adjust_m_proj[y,i,j] = fmin(500,exp(beta_t * (T_adjust_proj[i,y] - T_adjust_proj[j,y]))); 
            }
            
          }
        }
        tax_m_proj[y] = adj_m .* T_adjust_m_proj[y]; 
        
        tax_m_proj[y] = add_diag(tax_m_proj[y], -1 * colSums(tax_m_proj[y])); 
        
        
        mov_inst_m_proj[y] =  diff_m + tax_m_proj[y];
        
        mov_m_proj[y] = matrix_exp(mov_inst_m_proj[y]); 
        
        if ((sum(colSums(mov_m_proj[y])) / np - 1) > .001 ){
          print("Something has gone very wrong, movement matrix columns do not sum to 1");
          print(colSums(mov_m_proj[y]));
          print("width is", width);
          print("Topt is", Topt);
          print(diagonal(mov_inst_m_proj[y]));
        }
      }
    } 
    
    
    
  } // close patches
  
  
  for (y in 2:ny_proj){
    
    raw_proj[y] = normal_rng(0, 1); // draw a raw value 
    
    if(y==2){
      rec_dev_proj[y-1] = alpha * rec_dev[ny_train-1] +  sqrt(1 - pow(alpha,2)) *  sigma_r * raw_proj[y]; // initialize with last year of rec_dev
    }
    else{
      rec_dev_proj[y-1] =  alpha * rec_dev_proj[y-2] +  sqrt(1 - pow(alpha,2)) *  sigma_r * raw_proj[y]; 
    }
    
    // describe population dynamics
    for(p in 1:np){
      
      // density-independent, temperature-dependent recruitment of age 1
      
      if(T_dep_recruitment==1 && spawner_recruit_relationship==0){
        n_at_age_proj[y,p,1] = mean_recruits * exp(rec_dev_proj[y-1] - pow(sigma_r,2)/2) * T_adjust_proj[p,y-1] * beta_rec;
      }
      if(T_dep_recruitment==0 && spawner_recruit_relationship==0){
        n_at_age_proj[y,p,1] = mean_recruits * exp(rec_dev_proj[y-1] - pow(sigma_r,2)/2) ;
      }
      
      if(T_dep_recruitment==0 && spawner_recruit_relationship==1){
        n_at_age_proj[y,p,1] = (0.8 * r0 * h * ssb_proj[p, y-1]) / (0.2 * ssb0 * (1-h) + ssb_proj[p, y-1] * (h - 0.2));
        
        n_at_age_proj[y,p,1] =  n_at_age_proj[y,p,1] *  exp(rec_dev_proj[y-1] - pow(sigma_r,2)/2);
        
      }
      if(T_dep_recruitment==1 && spawner_recruit_relationship==1){
        n_at_age_proj[y,p,1] = ((0.8 * r0 * h * ssb_proj[p, y-1]) / (0.2 * ssb0 * (1-h) +  ssb_proj[p, y-1] * (h - 0.2))) * T_adjust_proj[p,y-1];
        
        n_at_age_proj[y,p,1] =  n_at_age_proj[y,p,1] *  exp(rec_dev_proj[y-1] - pow(sigma_r,2)/2);
      }
      // pop dy for non-reproductive ages 
      if(age_at_maturity > 1){ // confirm that there are non-reproductive age classes above 1
      
      n_at_age_proj[y,p,2:(age_at_maturity-1)] = (n_at_age_proj[y-1,p, 1:(age_at_maturity-2)]) .* to_row_vector(surv_proj[p,1:(age_at_maturity-2),y-1]);
      
      } // close if 
    } // close patches 
    
    // pop dy for reproductive adults
    
    if(T_dep_movement==0){
      
      for(p in 1:np){
        for(a in age_at_maturity:n_ages){
          // edge cases -- edges are reflecting
          if(p==1){
            n_at_age_proj[y,p,a] = n_at_age_proj[y-1,p, a-1] * surv_proj[p,a-1,y-1] * (1-d) + n_at_age_proj[y-1,p+1, a-1] * surv_proj[p+1,a-1,y-1] * d;
          } // close patch 1 case 
          
          else if(p==np){
            n_at_age_proj[y,p,a] = n_at_age_proj[y-1,p, a-1] * surv_proj[p,a-1,y-1] * (1-d) + n_at_age_proj[y-1,p-1, a-1] * surv_proj[p-1,a-1,y-1] * d;
          } // close highest patch
          
          else{
            n_at_age_proj[y,p,a] = n_at_age_proj[y-1,p, a-1] * surv_proj[p,a-1,y-1] * (1-2*d) + n_at_age_proj[y-1,p-1, a-1] * surv_proj[p-1,a-1,y-1] * d + n_at_age_proj[y-1,p+1, a-1] * surv_proj[p+1,a-1,y-1] * d;
            
          } // close if/else for all other patches
          
        }// close ages
      } // close patches 
    } // close T-dep movement if 
    
    // this code block calculates adult population size based on survival and directional movement 
    if(T_dep_movement==1){
      
      
      for (p in 1:np){
        
        n_at_age_proj[y,p, age_at_maturity:n_ages] = n_at_age_proj[y-1,p,  (age_at_maturity - 1):(n_ages - 1)] .* to_row_vector(surv_proj[p, (age_at_maturity - 1):(n_ages - 1),y-1]);
        
      }
      
      for(a in age_at_maturity:n_ages){
        
        // some acrobatics required here, because Stan won't do matrix multiplication with an array of reals like n_at_age_hat
        // instead we do the matrix multiplication with a placeholder vector and then populate n_at_age_hat 
        
        // fill in placeholder vector with reproductive ages across patches, and do mortality   
        for(p in 1:np){
          v_in_proj[p] = n_at_age_proj[y,p, a]; 
        }
        v_out_proj = mov_m_proj[y] * v_in_proj; // redistribute each age among patches according to the movement matrix 
        
        // fill in n_at_age_hat
        for(p in 1:np){
          
          n_at_age_proj[y,p,a] = v_out_proj[p]; 
          
        }
        
      } // close ages
    }// close T-dep movement if 
    
    for (p in 1:np){
      
      ssb_proj[p,y]  =  sum(to_vector(n_at_age_proj[y,p,1:n_ages]) .* maturity_at_age .* wt_at_age);
      
    }
    
    
  } // close year 2+ loop
  
  for(y in 1:ny_proj){
    
    
    for(p in 1:np){
      
      n_at_length_proj[y,p,1:n_lbins] = ((l_at_a_key' * to_vector(n_at_age_proj[y,p,1:n_ages])) .* selectivity_at_bin)'; // convert numbers at age to numbers at length. The assignment looks confusing here because this is an array of length y containing a bunch of matrices of dim p and n_lbins
      // see https://mc-stan.org/docs/2_18/reference-manual/array-data-types-section.html
      
      // NOTE, ignoring length sampling process at this point. In theory, need multinomial_rng or a custom dirichlet-multinomial rng here to generate simulated length comps
      
      n_at_length_obs_proj[y,p,1:n_lbins] =  n_at_length_proj[y,p,1:n_lbins]; // this is where length sampling process would go
      
      density_proj[p,y] = sum((to_vector(n_at_length_proj[y,p,1:n_lbins]))); // true population
      
      theta_proj[p,y] = ((1/(1+exp(-(beta_obs_int + beta_obs*log(density_proj[p,y] + 1e-6))))));
      
      density_obs_proj[p,y] = bernoulli_rng(theta_proj[p,y]) * exp(normal_rng(log(density_proj[p,y] + 1e-6), sigma_obs));
      // the observed densitieis as opposed to the true densities
      
    } // close patches 
    
    for(q in 1:number_quantiles){
      // calculate every range quantile q for every year y
      range_quantiles_proj[q, y] = calculate_range_quantile(np, patches, density_obs_proj[,y], quantiles_calc[q]);
    }
    
    
    centroid_proj[y] = sum(to_vector(density_obs_proj[,y]) .* patches) / sum(to_vector(density_obs_proj[,y])); // calculate center of gravity
    
    if(y>1) {
      abund_lr_proj[y-1] = log(sum(to_vector(density_obs_proj[,y])) / sum(to_vector(density_obs_proj[,y-1]))); 
    }
    
  } // close run forecast 
  
} // close generated quantities block


