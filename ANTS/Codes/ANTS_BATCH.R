######################################################################################################
### ANTS BATCH
###
### We're going to string together all the functions to try to track some ants movement
### All the functions come from CBPF_ANTS.R
###
### Much of this was written in a Starbucks in Old Town Alexandria as I waited for Em to get off work.
### Why am I getting nostalgic...
###
### Caleb Miller
######################################################################################################
set.seed(1)
library(png)

### Global Parameters Needed
rx<-30
ry<-10
nlevels<-10
nsamples<-40
steps <-100

########################
### Get The Original ###
########################

## Original Paramters
xc<-201
yc<-317
theta<- 4*pi/5


## Read First Image
frame<-image_read(1)
framexy<-frame_to_xy(frame)

## Template Rectangle (this has polygon but not graphing, I should put graphing in here)
tq <- get_rect(xc,yc,theta,framexy)

## Template Distribution
q <- get_p(xc,yc,theta,tq,nlevels)

#################################
### Push the Original Forward ###
#################################

## Get S and Pi

# This function will help form our S and Pi (*** put this in function file later)
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

S<-rep.row(c(xc,yc,theta),nsamples)

Pi<-matrix(1/nsamples,nsamples,1)

## Push it 


for(j in 1:steps){
# get new image
frame<-image_read(j)
framexy<-frame_to_xy(frame)
# use iterate step
# holder
receivelist<-list()
receivelist<-iterate(q,S,Pi,framexy,nsamples)
# results
S<-receivelist[[1]]
Pi<-receivelist[[2]]
Expec<-receivelist[[3]]

frame<-image_read(j)
get_rect(Expec[1],Expec[2],Expec[3],framexy)
}
