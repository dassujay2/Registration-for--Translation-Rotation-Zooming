remove(list = ls())
library(matlab) #Version 1.0.2
library(DRIP) #Version 1.7
library(jpeg) #Version 0.1.9
library(OpenImageR) #Version 1.2.1

packageVersion("matlab")
packageVersion("DRIP")
packageVersion("jpeg")
packageVersion("OpenImageR")


############################
##Reading the image into R##
############################
im1<-readJPEG("aral1.jpg")
im2<-readJPEG("aral_rz.jpg") ##Rotated and zoomed

###Adding translation componenet
im2<-translation(im2,shift_rows = 2,shift_cols = 2)


####################
##Edge detection####
####################
edge_ref<-which(stepEdgeLL2K(image=im1,bandwidth=2,thresh=0.15,plot=FALSE)==1,arr.ind = TRUE)
edge_zoomed<-which(stepEdgeLL2K(image=im2,bandwidth=2,thresh=0.15,plot=FALSE)==1,arr.ind = TRUE)

#################################################################################
#################################Nearest point to the center#####################
#################################################################################
edge_mod<-edge_zoomed

edge_dist<-NULL
for(i in 1:nrow(edge_mod))
{
  edge_dist<-c(edge_dist,sqrt((edge_mod[i,1]-64)^2+(edge_mod[i,2]-64)^2))
}

edge_dist_srt<-sort(edge_dist)

#############################################################################
######index of 20% edge-points having least distance from the origin################
#############################################################################
indx<-NULL
for (i in 1:243)
{
  indx<-c(indx,match(edge_dist_srt[i],edge_dist))
}  

edge_mod_eff<-edge_mod[indx,]    ##Effective edge points


#################################################
###########Registration Algorithm################
#################################################


#Here we are generating the registered image for the `Aral Sea' image for r1=8 and r2=16.

#For any other image and and any other r1 and r2 value, we can generate the regitered images accordingly.
#For that we need to change the `im1, im2, r1 and r2' accordingly.

w1<-30  ##pad-width


img<- im1  ##Reference image
img_zoom<- im2  ##Zoomed image

#Change the parameter values of r1 and r2 below.
r1=8
r2=16

m=proc.time()

img1<-img  
img_zoom1<-img_zoom  


mat1<- padarray(img_zoom1,c(w1,w1),"symmetric","both") ##Zoomed image (padded)
mat2<- padarray(img1,c(w1,w1),"symmetric","both")   ##Reference image (padded)


edg<-edge_mod_eff +w1
edge_reg_L1<-matrix(NA,nrow = nrow(edg),ncol = 2)
edge_reg_L2<-matrix(NA,nrow = nrow(edg),ncol = 2)

edge_x<-edg[,1]
edge_y<-edg[,2]

for(i in 1:nrow(edg))
{
  
  N1<-c()
  for(k1 in (edge_x[i]-r1):(edge_x[i]+r1))
  {
    for(l1 in (edge_y[i]-r1):(edge_y[i]+r1))
    {
      if( (k1-edge_x[i])^2 + (l1-edge_y[i])^2 <= (r1)^2)
      {
        N1<-c(N1,mat1[k1,l1])
      }
    }
  }
  msd<-c()
  mad<-c()
  MI<-c()
  # r<-c()
  xz<-c()
  yz<-c()
  # cov1<-0
  # sd1<-0
  # sd2<-0
  for(m1 in (edge_x[i]-r2):(edge_x[i]+ r2))
  {
    for(n1 in (edge_y[i]-r2):(edge_y[i]+ r2))
    {
      if( ((m1-edge_x[i])^2 + (n1-edge_y[i])^2) <= (r2)^2)
      {
        N2<-c()
        N<-0
        for(s1 in (m1-r1):(m1+r1))                              
        {
          for(t1 in (n1-r1):(n1+r1))
          {
            if( ((s1-m1)^2+(t1-n1)^2) <= (r1)^2)
            {
              N2<-c(N2,mat2[s1,t1])
              N<-N+1
            }
          }
        }

        msd<-c(msd,sum((N1-N2)^2)/N)
        mad<-c(mad,sum(abs(N1-N2))/N)
        xz<-c(xz,round(m1))
        yz<-c(yz,round(n1))
      }
    }
  }
  
  index_L2<-which.min(msd)
  index_L1<-which.min(mad)
  
  
  edge_reg_L1[i,1] <- xz[index_L1]
  edge_reg_L1[i,2] <- yz[index_L1]
  
  edge_reg_L2[i,1] <- xz[index_L2]
  edge_reg_L2[i,2] <- yz[index_L2]
  
  
}
edge_x_L1<-c()
edge_x_L1<-edge_reg_L1[,1]
edge_y_L1<-c()
edge_y_L1<-edge_reg_L1[,2]

edge_x_L2<-c()
edge_x_L2<-edge_reg_L2[,1]
edge_y_L2<-c()
edge_y_L2<-edge_reg_L2[,2]


################################################################
#########Estimation of rotation angle###########################
################################################################
a1<-matrix(NA,nrow = nrow(edg),ncol = 2)
a1<- edg-matrix(rep(c(mean(edg[,1]),mean(edg[,2])),nrow(edg)),ncol=2,byrow = TRUE)+ matrix(rep(c(mean(edge_reg_L1[,1]),mean(edge_reg_L1[,2])),nrow(edg)),ncol=2,byrow = TRUE) -w1

b1<-matrix(NA,nrow = nrow(edg),ncol = 2)
b1<-  edge_reg_L1 -w1

a2<-matrix(NA,nrow = nrow(edg),ncol = 2)
a2<- edg-matrix(rep(c(mean(edg[,1]),mean(edg[,2])),nrow(edg)),ncol=2,byrow = TRUE)+ matrix(rep(c(mean(edge_reg_L2[,1]),mean(edge_reg_L2[,2])),nrow(edg)),ncol=2,byrow = TRUE) -w1

b2<-matrix(NA,nrow = nrow(edg),ncol = 2)
b2<- edge_reg_L2 -w1


theta_L1<-NULL
theta_L2<-NULL
th1<-0
th2<-0
x1<-0 
x2<-0
y11<-0
y12<-0
y2<-0
y3<-0
for(i in 1:nrow(edg))
{
  x1<- sum(c(a1[i,1]-64,a1[i,2]-64)* c(b1[i,1]-64,b1[i,2]-64))
  x2<- sum(c(a2[i,1]-64,a2[i,2]-64)* c(b2[i,1]-64,b2[i,2]-64))
  
  y11<- sum(c(a1[i,1]-64,a1[i,2]-64)* c(a1[i,1]-64,a1[i,2]-64))
  y12<- sum(c(a2[i,1]-64,a2[i,2]-64)* c(a2[i,1]-64,a2[i,2]-64))
  
  y2<- sum(c(b1[i,1]-64,b1[i,2]-64)* c(b1[i,1]-64,b1[i,2]-64))
  y3<- sum(c(b2[i,1]-64,b2[i,2]-64)* c(b2[i,1]-64,b2[i,2]-64))
  
  theta_L1<- c(theta_L1,acos(x1/sqrt(y11*y2))*180/pi)
  theta_L2<- c(theta_L2,acos(x2/sqrt(y12*y3))*180/pi)
}
th1<- mean(theta_L1,na.rm = TRUE)
th2<- mean(theta_L2,na.rm = TRUE)


s_L1<-0
h_L1<-0
m_L1<-0


ex_L1<-0
ey_L1<-0
ex_L2<-0
ey_L2<-0
ex_L1<-(edge_x_L1-64-w1)
ey_L1<-(edge_y_L1-64-w1)
ex_L2<-(edge_x_L2-64-w1)
ey_L2<-(edge_y_L2-64-w1)


z1<-0
z2<-0

z1<- cos(th1*pi/180)*ex_L1+  sin(th1*pi/180)*ey_L1

z2<- (-1)*sin(th1*pi/180)*ex_L1+ cos(th1*pi/180)*ey_L1

A1<-matrix(rep(0,9),nrow = 3,ncol=3)

B1<-rep(0,3)
for(i in 1:nrow(edge_mod_eff))
{
  A1[1,1]<- A1[1,1]+ z1[i]^2+ z2[i]^2
  A1[1,2]<- A1[1,2]+z1[i]
  A1[1,3]<- A1[1,3]+z2[i]
  A1[2,2]<- A1[2,2]+1
  B1[1]<- B1[1]+ z1[i]*(edg[i,1]-64-w1)+ z2[i]*(edg[i,2]-64-w1)
  B1[2]<- B1[2]+ (edg[i,1]-64-w1)
  B1[3]<- B1[3]+ (edg[i,2]-64-w1)
}

A1[3,3]<-A1[2,2]
A1[2,1]<-A1[1,2]
A1[3,1]<-A1[1,3]

sol1<- rep(0,3)
sol1<- solve(A1)%*%B1


s_L1= sol1[1]
h_L1= sol1[2]
m_L1= sol1[3]



reg_L1<- matrix(NA, nrow = nrow(mat1),ncol = ncol(mat1))
i1<-0
j1<-0
for(i in 1:(nrow(mat1)))
{
  for(j in 1:(nrow(mat1)))
  {
    i1<- ceiling(((i-64-h_L1-w1)/s_L1)*cos(th1*pi/180) - ((j-m_L1-64-w1)/s_L1)*sin(th1*pi/180))+64+w1
    j1<- ceiling(((i-64-h_L1-w1)/s_L1)*sin(th1*pi/180) + ((j-m_L1-64-w1)/s_L1)*cos(th1*pi/180))+64+w1
    if(is.na(i1)==TRUE)
      i1<-i
    else if(i1<1)
      i1<-1
    else if(i1>nrow(mat1))
      i1<-nrow(mat1)
    
    if(is.na(j1)==TRUE)
      j1<-j
    else if(j1<1)
      j1<-1
    else if(j1>nrow(mat1))
      j1<-nrow(mat1)
    reg_L1[i1,j1]<-mat1[i,j]
  }
}

t1<-0

n1<-0

for(i in (w1+r2):(nrow(mat2)-w1-r2+1))
{
  for(j in (w1+r2):(nrow(mat2)-w1-r2+1))
  {
    if(is.na(reg_L1[i,j])==FALSE)
    {
      t1<-t1+(mat2[i,j]-reg_L1[i,j])^2
      n1<-n1+1
    }
    
  }
}
msd_L1<-t1/n1  ##MSE under L1-norm

s_L2<-0
h_L2<-0
m_L2<-0

z3<-0
z4<-0


z3<- cos(th2*pi/180)*ex_L2+  sin(th2*pi/180)*ey_L2

z4<- (-1)*sin(th2*pi/180)*ex_L2+ cos(th2*pi/180)*ey_L2

A2<-matrix(rep(0,9),nrow = 3,ncol=3)

B2<-rep(0,3)
for(i in 1:nrow(edge_mod_eff))
{
  A2[1,1]<- A2[1,1]+ z3[i]^2+ z4[i]^2
  A2[1,2]<- A2[1,2]+z3[i]
  A2[1,3]<- A2[1,3]+z4[i]
  A2[2,2]<- A2[2,2]+1
  B2[1]<- B2[1]+ z3[i]*(edg[i,1]-64-w1)+ z4[i]*(edg[i,2]-64-w1)
  B2[2]<- B2[2]+ (edg[i,1]-64-w1)
  B2[3]<- B2[3]+ (edg[i,2]-64-w1)
}

A2[3,3]<-A2[2,2]
A2[2,1]<-A2[1,2]
A2[3,1]<-A2[1,3]

sol2<- rep(0,3)
sol2<- solve(A2)%*%B2

s_L2= sol2[1]
h_L2= sol2[2]
m_L2= sol2[3]


reg_L2<- matrix(NA, nrow = nrow(mat1),ncol = ncol(mat1))
i2<-0
j2<-0
for(i in 1:(nrow(mat1)))
{
  for(j in 1:(nrow(mat1)))
  {
    i2<- ceiling(((i-64-h_L2-w1)/s_L2)*cos(th2*pi/180) - ((j-64-m_L2-w1)/s_L2)*sin(th2*pi/180))+64+w1
    j2<- ceiling(((i-64-h_L2-w1)/s_L2)*sin(th2*pi/180) + ((j-64-m_L2-w1)/s_L2)*cos(th2*pi/180))+64+w1
    if(is.na(i2)==TRUE)
      i2<-i
    if(i2<1)
      i2<-1
    else if(i2>nrow(mat1))
      i2<-nrow(mat1)
    
    if(is.na(j2)==TRUE)
      j2<-j
    if(j2<1)
      j2<-1
    else if(j2>nrow(mat1))
      j2<-nrow(mat1)
    
    
    reg_L2[i2,j2]<-mat1[i,j]
  }
}


t2<-0

n2<-0

for(i in (w1+r2):(nrow(mat2)-w1-r2+1))
{
  for(j in (w1+r2):(nrow(mat2)-w1-r2+1))
  {
    if(is.na(reg_L2[i,j])==FALSE)
    {
      t2<-t2+(mat2[i,j]-reg_L2[i,j])^2
      n2<-n2+1
    }
  }
}
msd_L2<-t2/n2 ###MSE under L2-norm

n=proc.time()

n-m


#############################################
######Creation of registered image#########
#############################################


L1_img<-matrix(NA,nrow = nrow(mat1)-2*w1,ncol = ncol(mat1)-2*w1)  ##Registered image under L1-norm
L1_img[(r2+1):((nrow(L1_img)-r2)),(r2+1):((ncol(L1_img)-r2))]<-reg_L1[(w1+r2+1):((nrow(mat1)-w1-r2)),(w1+r2+1):((ncol(mat1)-w1-r2))]

L2_img<--matrix(NA,nrow = nrow(mat1)-2*w1,ncol = ncol(mat1)-2*w1) ####Registered image under L2-norm
L2_img[(r2+1):((nrow(L2_img)-r2)),(r1+1):((ncol(L2_img)-r2))]<-reg_L2[(w1+r2+1):((nrow(mat1)-w1-r2)),(w1+r2+1):((ncol(mat1)-w1-r2))]


############################################
###########Residual Image###################
############################################
res_L1<- img - L1_img  ##residual image under L1-norm

res_L2<- img - L2_img  ##residual image under L2-norm


############################################
#############Visualization##################
############################################
image(rot90(L1_img,3),col = grey(seq(0,1,length=256)))

#################################
######### MSE ###################
#################################
msd_L1
msd_L2

################################
#### Parameter Estimates #######
################################

####### For L1-norm ############
# th1      ##Estimate of rotation angle   
# s_L1     ##Estimate of zooming factor
# h_L1     ##Estimate of x-coordinate of the translation parameter
# m_L1     ##Estimate of y-coordinate of the translation parameter
