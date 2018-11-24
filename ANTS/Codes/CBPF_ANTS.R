#### Applying the framework of "An adaptive color-based particle filter" to our pngFrames of Ants
#### Caleb Miller

#### Functions in (roughly) order of appearance in the paper 
#### Plus auxilliary function needed for processing

##############################################################################
### Pixel Weighting function 
### k(xc,yc,x,y,theta)
### Globals needed: rx, ry
### Aux functions: rotation
### Inputs: 
###   xc - x-coordinate of the rectangles center
###   yc - y-coordinate of the rectangles center
###   x - x coordinate of pixel being investigated
###   y - y coordinate of pixel being investigated
###   theta - angle off x-axis of rectangle
### Outputs:
###   val - weighted value
##############################################################################

# Aux: rotation function
rotate<-function(theta){
  R<-matrix(0,2,2)
  R[1,1] <- cos(theta)
  R[1,2] <- -sin(theta)
  R[2,1] <- -R[1,2]
  R[2,2] <- R[1,1]
  R  
}

# weight function should have rx and ry predefined as well as rx^2+ry^2
k<-function(xc,yc,theta,x,y){
  centpoint=rotate(-theta)%*%rbind(x-xc,y-yc)
    xstar=centpoint[1,]
    ystar=centpoint[2,]
  val = (rx^2+ry^2-xstar^2-ystar^2)/(rx^2+ry^2)
  return(val)
}


##############################################################################
### Form Tracking Rectangle Data Frame
### Globals needed: rx,ry
### Aux functions: k, frame_to_xy
### Inputs: 
###   xc - x-coordinate of the rectangles center
###   yc - y-coordinate of the rectangles center
###   theta - angle off x-axis of rectangle
###   framexy - converted xy
### Outputs:
###   tr - data frame
###     tr$color - contains the color values of the tracking rectangle
###     tr$x - contains the interior x-points of the tracking rectangle
###     tr$y - contains the interior y-points of the tracking rectangle
###     tr$bins - labeled bins for the data
##############################################################################

# Aux frame_to_xy
frame_to_xy<-function(frame){
  framexy<-t(frame1[nrow(frame):1,])
  return(framexy)
}

get_rect<-function(xc,yc,theta,framexy){
  
  # Parameters used to form rectangle with center at (0,0)
  trackboxcenter<-c(xc,yc)
  trackboxwidth<-rx
  trackboxheight<-ry
  trackboxangle<-theta
  
  ## Graphing the Tracking Rectangle
  
  # Corners centered at (0,0)
  trackboxx<-c(-0.5*trackboxwidth,0.5*trackboxwidth,0.5*trackboxwidth,-0.5*trackboxwidth)
  trackboxy<-c(-0.5*trackboxheight,-0.5*trackboxheight,0.5*trackboxheight,0.5*trackboxheight)
  
  # Rotate by given theta and recenter using xc and yc
  trackpoints<-rotate(trackboxangle)%*%as.matrix(rbind(trackboxx,trackboxy))+trackboxcenter
  
  # Connect corners for the red graphed rectangle 
  polygon(trackpoints[1,],trackpoints[2,],border="red")
  
  ## Accessing the actual points inside the tracking rectangle
  
  # Get points inside of the tracking rectangle
  interiorpointsx <-seq(-0.5*trackboxwidth,0.5*trackboxwidth-1,1)
  interiorpointsy <-seq(-0.5*trackboxheight,0.5*trackboxheight-1,1)
  interiorpoints  <-expand.grid(interiorpointsx,interiorpointsy)
  
  interiorpoints  <-rotate(trackboxangle)%*%as.matrix(t(interiorpoints))+trackboxcenter
  interiorpointsx <-matrix(interiorpoints[1,],trackboxheight,trackboxwidth,byrow = T)
  interiorpointsy <-matrix(interiorpoints[2,],trackboxheight,trackboxwidth,byrow = T)
  
  ## Averaging and Binning
  
  # Create the pixel-averaged rectangle of points
  track_rect<-matrix(0,trackboxheight,trackboxwidth)
  fracx<-interiorpointsx-floor(interiorpointsx)
  fracy<-interiorpointsy-floor(interiorpointsy)
  
  for(i in 1:trackboxwidth){
    for(j in 1:trackboxheight){
      trakc_rect[j,i]<-(1-fracx[j,i])*(1-fracy[j,i])*framexy[interiorpointsx[j,i],interiorpointsy[j,i]] +fracx[j,i]*(1-fracy[j,i])*framexy[interiorpointsx[j,i]+1,interiorpointsy[j,i]]+(1-fracx[j,i])*fracy[j,i]*framexy[interiorpointsx[j,i],interiorpointsy[j,i]+1]+fracx[j,i]*fracy[j,i]*framexy[interiorpointsx[j,i]+1,interiorpointsy[j,i]+1]
    }  
  }
  
  
  # Bin everything up (maybe be adjusted later for better representation but this just cuts up evenly)
  delta <- 1/10
  nlevels <- 1/delta
  cut_temp <- cut(track_rect, breaks = seq(0.0,1.0,delta))
  
  
  # Flatten data to put into the data frame
  dim_tot <- trackboxwidth*trackboxheight
  
  dim(template)        <- c(dim_tot,1)
  dim(interiorpointsx) <- c(dim_tot,1)
  dim(interiorpointsy) <- c(dim_tot,1)
  dim(cut_temp)        <- c(dim_tot,1)
  
  # Make that data frame (tr for tracking rectangle)
  tr<- data.frame(
    color  <- track_rect,
    x      <- interiorpointsx,
    y      <- interiorpointsy,
    bins   <- cut_temp
  )
  return(tr)
}

##############################################################################
### Get Color distribution of the interior of the tracking rectangle
### Globals needed: rx,ry
### Aux functions: k, get_rect
### Inputs: 
###   tr - x-coordinate of the rectangles center
### Outputs:
###   p - the color distribution
##############################################################################

# See about combining xc,yc,theta into the tr data frame

get_p<-function(xc,yc,theta,tr){
  
  
  # The color distribution to be filled
  p <- matrix(0,nlevels,1)
  
  # Get bin_lvls for logical check in the loop
  bins_lvl <- levels(tr$bins)
  
  # For each level sum weights of points assigned to that level
  for(i in 1:nlevels){
    cond <- tf$bins == bins_lvl[i]
    p[i] <- sum(k(xc,yc,theta,tr$x[cond],tr$y[cond]))
  }
  
  # Normalize
  p=p/(sum(p))
  
  return(p)
}

##############################################################################
### Bhattacharyya Distance
### Inputs: 
###   p - proposed distribution
###   q - template distribution
### Outputs:
###   d - Bhattacharyya Distance
##############################################################################

get_d<-function(p,q){
  # Get Bhattacharryya Coefficient
  rho <- sum(sqrt(p*q))
  # Get Bhattacharryya Distance
  d <- sqrt(1-rho)
  return(d)
}

##############################################################################
### Motion Model
### Sourced from "MCMC-Based Particle Filtering..." - Zia Khan 
### Inputs: 
###   xc - x- coordinate of center of rectangle
###   yc - y-coordinate of center of rectangle
###   theta - angle of head of ant rotated from postive x-axis, origin at xc,yc
### Outputs:
###   new - new[1] =xc, new[2]=yc, new[3]=theta
###     xc - updated/ pushed through motion model
###     yc - updated
###     theta - updated
##############################################################################

motion<-function(xc,yc,theta){
  # Make center
  center = c(xc,yc)
  
  # Draw changes from Normal distributions 
  delta_x      = rnorm(1,0,4)
  delta_y      = rnorm(1,0,8)
  delta_theta  = rnorm(1,0,0.4)
  
  # Get updates
  updates = rotate(theta+delta_theta)%*%c(delta_x,delta_y)
  new = center+updates
  
  xc=new[1]
  yc=new[2]
  theta = theta+delta_theta
  new<-c(new,theta)
  return(new)
}


##############################################################################
### Sample Weighting Function
###  
### Inputs:
###   q - template distribution
###   p - proposed distribution
### Outputs:
###   spi - sample weight
##############################################################################

get_spi<-function(p,d,q){
  # dont know what to put for this sigma at the moment
  simga_sq <- 1
  spi<-1/sqrt(2*pi*sigma_sq)*exp(-d^2/(2*sigma_sq))
}


##############################################################################
### Selection of Samples  
###  
### Inputs:
###   S - array of samples
###   Pi - array of weights
### Outputs:
###   S_prime - selected samples
##############################################################################

select <-function(S,Pi){
  # Length
  N=length(S)
  # Get cumulative probabilities and normalize
  cp <- matrix(0,N+1,1)
  for(i in 1:N){
    cp[i+1]=sum(cp)+Pi[i]
    print(cp)
  }
  cp<-cp/cp[N+1]
  
  # Get uniforms
  U <- runif(N,0,1)
  
  # Bin up data based and get counts of data in each bin
  hist_u<-hist(U, breaks = cp,plot = FALSE)
  bin_counts<-hist_u$counts
  
  # Form the new samples
  s_prime<- matrix(0,N,1)
  
  # For each number in the bin count total a sample of that is to appear in s_prime
  j<-1
  for(i in 1:N){
    while(bin_counts[i]!=0 && j<N+1){
      bin_counts[i]=bin_counts[i]-1
      s_prime[j] <- s[i]
      j<-j+1
    }
  }
  return(s_prime)
}

##############################################################################
### Iteration Step
### Functions needed: CBPF_ANTS has all the functions
### Have all the frames loaded...probably
### Inputs: 
###   q - the template distribution
###   S - samples (N-by-3) xc,yc,theta columns 
###   Pi - weights for proposed distributions
###   step_num - what step of the iteration are we on?
### Outputs:
###   
##############################################################################

iterate<- function(q,S,Pi,framexy){
  # Get length
  N = length(S)
  # SELECT samples/proposals based on weights
  S_prime <- select(S,Pi)
  # PROPOGATE the samples forward 
  for(i in 1:N){
    S_prime[i,]<-motion(S_prime[i,1],S_prime[i,2],S_prime[i,3])
  }
  # OBSERVE color distributions
  # Get tracking rectangles get_rect to get the rectangles (*** mbc to do some plots here)
  rects<-matrix(0,N,1)
  for(i in 1:N){
    rects[i]<-get_rect(S_prime[i,1],S_prime[i,2],S_prime[i,3],framexy)
  }
  # Get color distributions (*** seconds dimension (10) could be changed to vairable because of the bin#s)
  dists<-matrix(0,N,10)
  for(i in 1:N){
    dists[i]<-get_p(S_prime[i,1],S_prime[i,2],S_prime[i,3],rects[i])
  }
  # Get Bhattacharyya distance
  ds<-matrix(0,N,1)
  for(i in 1:N){
    ds[i]<-get_d(dists[i],q)
  }
  # Weight each sample
  for(i in 1:N){
    Pi[i]<-get_spi(dists[i],ds[i],q)
  }
  # ESTIMATE
  Expec<-sum(Pi*S_prime)
  
  return(S_prime,Pi,Expec)
}

