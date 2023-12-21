
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.io import fits
import scipy.misc 
import math as mt
from jdcal import gcal2jd,jd2gcal
import os
import timeit
import scipy as sp
import cv2
import pathlib
import numpy.ma as ms

from skimage.morphology import erosion, dilation, opening, closing, white_tophat
from skimage.morphology import black_tophat, skeletonize, convex_hull_image
from skimage.morphology import disk




#from skimage import data
from skimage.morphology import  closing, square
from scipy.ndimage import label, generate_binary_structure,find_objects,measurements,map_coordinates

startTime = timeit.default_timer()
totelIm=0

#fol = "/media/adithya/Sunku/sundata/1990/01";  #paste your foleder location
#spl=fol.split('/')

pathlib.Path('image_data').mkdir(parents=True, exist_ok=True) 
pathlib.Path('plage_data').mkdir(parents=True, exist_ok=True) 
pathlib.Path("Extracted_images").mkdir(parents=True, exist_ok=True) 

for x in range(1,2,1):
 p=str('%02.0f' %x)
 newlist = []
 thArry=[]
 lonArry=[]
 area_no=[]
 Area=[]
 DelTh=[]
 DelLn=[]
 DelA=[]
 DateArray=[]
 pixArray=[]
 lonCArray=[]
 
 Da_=[]
 h_=[]
 k_=[]
 JD_=[]
 P_=[]
 Bo_=[]
 Lo_=[]
 Rad_=[]
 Pno_=[]
 S3_=[]
 folderTosearch ="/mnt/A6DA8B10DA8ADBC5/PROJECT/1984/{}".format(p);  
 for root, dirs, files in os.walk(folderTosearch):
    for file in files:
        if file.endswith(".fits.gz"):
          newlist.append(file)
        
          
 newlist.sort()
 Length=len(newlist)  
 print('______Month:',p,'______')      
 print('No. of Images found :',Length)
 print('')
 print('_____')
 a =np.zeros((4096))
 for sir in range(Length):
  img=fits.open(os.path.join(root,newlist[sir]))
  scidata=img[0].data #for image matrix
  
  #xs=np.linspace(1,4096,4096) #to plot graph
  #ac=np.zeros((4096)) #to graph
  #xi=[]
  #yi=[]  #to add bright pixels
  #ec=np.zeros((4096,4096))
  
  m=scidata.mean()
  halfmean=m*0.5 
  clip=np.clip(scidata,0,halfmean)
  ec=np.absolute(sp.ndimage.filters.sobel(clip))
  edge=np.clip(ec,m,2*m) #
  c=np.nonzero(edge==2*m)#(row,col)
  yi=c[0]
  xi=c[1]
  
  #scipy.misc.imsave('grad{}.jpg'.format(sir),edge)
  #plt.plot(xs,ac)
  #plt.plot(xi,yi,'r*')
  #plt.show()
  
  x=xi
  y=yi
  N=len(x)
 
  X=np.sum(x)
  Y=np.sum(y)
  x_ =X/N
  y_ =Y/N
  
  u=x-x_
  v=y-y_
  
  ui2=np.sum(u**2)
  vi2=np.sum(v**2)
  uvi=np.sum(u*v)
  ui3=np.sum(u**3)
  vi3=np.sum(v**3)
  uiv2=np.sum(u*(v**2))
  viu2=np.sum(v*(u**2))
  
  c1=0.5*(ui3+uiv2)
  c2=0.5*(vi3+viu2)
  
  a=np.zeros(shape=(2,2))

  a[0][0]=ui2
  a[0][1]=uvi
  a[1][0]=uvi
  a[1][1]=vi2
  b=[c1,c2]
  cen=np.linalg.solve(a,b)
  uc=cen[0]
  vc=cen[1]
  centre=uc+x_,vc+y_ #x,y
  R=((uc)**2+(vc)**2+(1/N)*(ui2+vi2))**0.5
  
  #print(R,centre)
  
  Date=img[0].header['DATE-OBS']
  DOB=[]
  DOB.append(Date)
  d=[]
  d+=Date
  print('DOB :',Date)
  print("centre:",centre)
  print("Radius:",R)
  
  yr=[''.join(d[0:4])]
  mo=[''.join(d[4:6])]
  Da=[''.join(d[6:8])]
  t=[''.join(d[9:11])]
  mi=[''.join(d[11:13])]

  DObStr=[''.join(d[0:8]),''.join(d[9:13])]
  DnTstr=[''.join(DObStr[0:2])]
  DnMStr=[''.join(d[0:6])]
  #DObStr=[''.join(DObStr[0:2])] #to include time
  #timeStr=[''.join(d[9:13])]
  DOb=int(DObStr[0])
  DnM=int(DnMStr[0])
  DnT=int(DnTstr[0])
  #print(DnT)
  #time=int(timeStr[0])

  M=int(mo[0])
  Y=int(yr[0])
  D=int(Da[0])
  Hr=int(t[0])
  Min=int(mi[0])

  d=((Hr/24)+(Min/(24*60)))   #+sec/(24*60*60)  #its in ut..already corrected           

  c=gcal2jd(Y,M,(D))
  JD=c[0]+c[1]+d
  #print('Julian date:',JD)
  
  
  T=((JD)-(2415020))/(36525)
#This is for before 2000.
  L_=np.deg2rad(279.69668)+np.deg2rad(36000.76892)*T+np.deg2rad(.0003025)*(T**2)
  L_=np.remainder(np.rad2deg(L_),360) 
  L_=np.deg2rad(L_)	
  g=np.deg2rad(358.47583)+np.deg2rad(35999.04975)*T-np.deg2rad(0.00015)*(T**2)-np.deg2rad(.0000033)*(T**3) 
  g=np.remainder((np.rad2deg(g)),360)
  g=np.deg2rad(g)
  Ii=np.deg2rad(7.25)
  omg=np.deg2rad(259.18)-np.deg2rad(1934.142)*T
  C=(np.deg2rad(1.91946)-np.deg2rad(.004789)*T-np.deg2rad(0.000014)*(T**2))*mt.sin(g) + (np.deg2rad(0.020094)-np.deg2rad(.0001)*T)*mt.sin(2*g)  +np.deg2rad(0.000293)*mt.sin(3*g)
  Ltr=L_+C
  Ltr=np.remainder(np.rad2deg(Ltr),360)
  Ltr=np.deg2rad(Ltr)	
  La=Ltr-np.deg2rad(.00569)-np.deg2rad(.00479)*mt.sin(omg)
  fi=np.remainder((JD-2398220)*(360/25.38),360)
  K=np.deg2rad(74.3646)+np.deg2rad(1.395833)*T
  epo=np.deg2rad(23.452295)-np.deg2rad(0.0130125)*T-np.deg2rad(.00000164)*(T**2)+np.deg2rad(0.000000503)*(T**3)
  ep=epo+np.deg2rad(.00256)*mt.cos(omg)
  X=mt.atan(-mt.cos(La)*mt.tan(ep))
  Y=mt.atan(-mt.cos(Ltr-K)*mt.tan(Ii))

#position of poles
  P=X+Y
  Bo=mt.asin(mt.sin(Ltr-K)*mt.sin(Ii))
  M=(360-fi)
  M=(np.deg2rad(M))
  y=mt.sin(K-Ltr)*mt.cos(Ii)
  x=-(mt.cos(K-Ltr))
  Lo=mt.atan2(y,x)+M

#pre defined
  R_=1.00014-0.01671*mt.cos(g)-0.00014*mt.cos(2*g)
  Rad=(0.2684/R_)*3600                      #0.2684-rad of sun plus 1800km chromosphere
  T=(Rad/15)
  Ro=29.5953*mt.cos((mt.acos(-.00629*T)/3+np.deg2rad(240)))	
 	
  h=int(round(centre[0])) #x                #centre coordinates
  k=int(round(centre[1])) #y
  Rr=int(round(R)-50)                       #R-50,  -50 is to avoid edge 
  aa=np.zeros((4096,4096))                  #mtx for img
  circ=cv2.circle(aa,(h,k),Rr,(255,0,0),-1) #filed
  circle=circ.astype(np.bool)
  Circle=np.invert(circle)                  #empty inside
  mask=ms.array(scidata,mask=Circle)        #with disk
  disk=ms.array(scidata,mask=circle)        #masked disk
  DD=disk*0
  sun=DD.data                               #only solar disk
  index=np.nonzero(sun)
  mtxData=sun[index[0],index[1]]
  
  #print(mtxData.mean(),np.median(mtxData))
  
  s=np.std(mask)
  d=np.median(mtxData)                      #np.mean(mask)
  print('median',d)
  s4=4*s
  s3=3*s
  S4=d+s4
  S3=d+s3      
  print('threshold intensity :',S3) 	#midean+3sigma
  thresh =S3
  im = closing(sun > thresh, square(3)) 
  
  #scipy.misc.imsave('thresh.jpg',im)

#first labeling
  S = generate_binary_structure(2,2)
  lba, num_features = label(im,structure=S)
  print(num_features)

  no_of_plage=0
  Roo=np.deg2rad(Ro)
  sundata=np.zeros((4096,4096))
  for lbl in range (num_features):
    loc=find_objects(lba)[lbl]
    area=np.count_nonzero(lba==lbl+1)
    if area>311 :
      sundata[loc]=lba[loc]
      plage=ms.masked_where(lba[loc]!=lbl+1,lba[loc])
      yy,xx=np.nonzero(plage)
      x=np.array(loc)
      Col=(x[1].start)
      Row=(x[0].start) 
      i=(xx+Col) #cartissian x
      j=(yy+Row) # y  
      r=((j-k)**2+(i-h)**2)**0.5
      y=(j-k)
      x=(i-h)     
      th_ =np.arctan2(y,x)
      ro_=np.deg2rad(Ro*(r/R))  
      ro=np.arcsin(np.sin(ro_)/np.sin(Roo))-(ro_)
      th=np.arcsin(np.cos(ro)*np.sin(Bo)+np.sin(ro)*np.cos(Bo)*np.sin(th_))
      l= np.arcsin(np.cos(th_)*np.sin(ro)/np.cos(th))
      L= l+Lo
      ld=np.rad2deg(l)
      th=np.rad2deg(th)
      Ld=np.rad2deg(L)
      LD=np.mod((Ld),360)
      numl=np.sum(ld*sun[j,i])
      Nu=np.sum(LD*sun[j,i])
      nu=np.sum(th*sun[j,i])
      den=np.sum(sun[j,i])
      Lt=(nu/den)
      Ln=Nu/den
      ln=numl/den
      Del=abs(np.sin(Bo)*np.sin(Lt)+np.cos(Bo)*np.cos(Lt)*np.cos(ln))
      Carea=(area/Del)
      HArea=2*np.pi*(R)**2/10e5
      Ccarea=(Carea/ HArea)  #3.04e6) 
      Sdth=np.std(th)
      Sdln=np.std(LD)
      Dlth=Sdth/(N**0.5)
      Dlln=Sdln/(N**0.5)
      area_no.append(lbl+1)  
      no_of_plage+=1
      area_no.append(lbl+1)
      pixArray.append(area)
      Area.append(Ccarea)#Ccarea)
      thArry.append(Lt)
      lonArry.append(Ln)
      DelTh.append(Dlth)
      DelLn.append(Dlln)
      lonCArray.append(ln)
      DateArray.append(DnT)
  ##
  Edge=closing(edge > m, square(3))
  data=np.add(Edge,sundata)
  dataImage=np.clip(data,0,1)
  scipy.misc.imsave('Extracted_images/sun{}.jpg'.format(DOb),dataImage[::-1])    #[::-1] to flip it in y axis
  ##
  P=np.rad2deg(P)
  Bo=np.rad2deg(Bo)
  Lo=np.rad2deg(Lo)
  Lo=np.remainder((Lo),360)
  
  print('P :',P,'Bo :',Bo,'Lo :',Lo)
  Da_.append(DnT)
  JD_.append(JD)
  h_.append(h)
  k_.append(k)
  Rad_.append(R)
  S3_.append(S3)
  P_.append(P)
  Bo_.append(Bo)
  Lo_.append(Lo)
  Pno_.append(no_of_plage)  
  
  print('No. of Plages:',no_of_plage)
  f=open('image_data/Ca-II_image_data{}.dat'.format(DnM),'a')
  np.savetxt('image_data/Ca-II_image_data{}.dat'.format(DnM),np.c_[Da_,JD_,h_,k_,Rad_,S3_,P_,Bo_,Lo_,Pno_], fmt='%9.2f',header='     DOB         JD        centre-X  centre-Y  Radius   Thresh Int    P          Bo       Lo     No. of Plages')
  f.close()
  print(sir+1,"done out of", Length)
  f=open('plage_data/ca-II_Plage_data{}.dat'.format(DnM),'a')
  np.savetxt('plage_data/ca-II_Plage_data{}.dat'.format(DnM),np.c_[DateArray,Area,thArry,lonArry,DelTh,DelLn, lonCArray,pixArray], fmt='%9.2f' ,header='    DOB           Area      Hlat      Hlon        dth      dln       lonC   ' )
  f.close()
    
 
  print('------:::::######:::::------')
 totelIm+=Length
 print('______________________________________________________')  
stopTime = timeit.default_timer()
runtime=stopTime-startTime
print('')
print('.......COMPLEETED........')
print("Totel No. of Images",totelIm)
print('Run Time',runtime/60 ,'Min')
print("Avg. time per image",runtime/totelIm,'Sec')


