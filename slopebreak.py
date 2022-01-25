# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 13:02:13 2020

@author: tonian
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from linearstatistics import linestatistics as lstats

pckfile = '308.txt'
pckdata = pd.read_csv(pckfile, sep="\s", header=None)
show_slopebreaks = 'true'

pckdata.columns = ["dist", "time"]

x = np.array(pckdata.dist)
y = np.array(pckdata.time)

nPts = len(x); #number of Points in each time series

#Plotting 
plt.plot(x,y, 'b-o', label='1st Arrival Pcks')
#plt.title("Displacment for PS Point %i" %(in_ID))
plt.ylabel('Time (ms)',fontsize=15)
plt.xlabel('Distance (m)',fontsize=15)


# Define range of "kinks" to test for a slope break in best fit line
firstkink = 3; #allowing the kinks to start after the first 3 data Points (remember python indices start at 1 so the first kink is at #4); this is the location of the kink as it relates to the extracted 56 data Pointss
lastkink = nPts-3;#:nPts; #leaving 4 data points for the last kink
kinks = np.arange(firstkink, lastkink+1);
nKinks = len(kinks);

#------------CREATING ARRAYS TO STORE INFORMATION----------------------#
b = np.zeros((nKinks,4)); #stores best-fitting parameters from 2-line fit for all positions
yfitk = np.zeros((nPts,1)); #empty array for all kinked lines
yfitminslopeuncertaintyk = np.zeros(nKinks); #empty array for y values of minslopeuncertaintyk (shows optimum kink point)
minslopeuncertaintyk = np.zeros(2); #creating empty array to store max Rsquared values for kinked lines
kinkpos = np.zeros(1); #creating empty array to store corresponding time values for max Rsquared values
yfitkinkline = np.zeros((1,nPts)); #empty array for saving y fit values caluclated with kink location (basically plots kink line) 
bminslopeuncertaintyk = np.zeros(4); #b values that correspond to the min slopeuncertaintyk value #b(1) intercept of line 1 for nth position
#b(2)#slope of line 1, b(3)intercept of second line added past kink Points, b(4) slope of second line added past kink Points
A = np.zeros((nPts,4))
D = []

#--------Best Fit Two Lines With Kink------------------------
bk = np.zeros((1,4)); # make bk array to store best-fit coefficients for each kink time
statistics_k1 = np.zeros((nKinks,5)) # make statistics array to for each kink time for 1st slope
statistics_k2 = np.zeros((nKinks,5)) # make statistics array to for each kink time for 2nd slope


#-----------Begin finding slope breaks------------------------------------
nk = -1; # make counter for kinks
for ik in range(firstkink,lastkink+1):
    nk = nk+1;
    D = np.zeros(nPts);
    D[ik:nPts+1] = 1;  # for Pointss with D = 1 second line is added; letting the first 3 points = 0
    
    # build A array , for which y = A * b and b = A\y 
    A[:,0] = np.ones(nPts)
    A[:,1] = x
    A[:,2] = D
    A[:,3] = np.multiply(D,x)
    
    bk = np.linalg.lstsq(A,y)[0]; #calculating b array using the linear algebra - least squares from numpy (has to be least squares becasue the A matrix is not a perfect square)
    #input('Press ENTER to continue');
    
    #store b (array of 4) into the bk row for this kink time
    b[nk,:] = bk; #using the corresponding b array to calculate the yfit value for each Points
    
    #compute the yfit best-fitting line with kink through the data
    yfitk = bk[0] + bk[1]*x + bk[2]*D + bk[3]*D*x
    
    if show_slopebreaks == 'true':
        # -------Plot kinks------
        #random colors for fun
        rc = np.random.rand(1)
        bc = np.random.rand(1)
        gc = np.random.rand(1)
        color = (rc[0], gc[0], bc[0])
        plt.plot(x,yfitk, c=color, linewidth=1,)
    elif show_slopebreaks == 'false':
        pass
    
    # Compute line statistics for 1st slope, function returns [m, b, m_uncertainty, b_uncertainty,  R2]
    statistics_k1[nk,:] = lstats(x[0:ik],y[0:ik]) 
    
    # Compute slope uncertainty for 2nd slope, function returns [m, b, m_uncertainty, b_uncertainty,  R2]
    statistics_k2[nk,:] = lstats(x[ik-1:nPts],y[ik-1:nPts]) # want points to cross over the "kink point"
    

#-------------end of for loop for for finkding slope breaks---------------------    

    
#Finding minimum combination of slope R2 for all kink lines
sq_combo = np.square(statistics_k1[:,2]) + np.square(statistics_k2[:,2])
[mincombo_min,mincombo_index] = [np.amin(sq_combo),np.where(sq_combo == sq_combo.min())];
#np.where(a == a.min())

#Storing minimum slope uncertainties
minslopeuncertaintyk[0] = statistics_k1[mincombo_index,2];
minslopeuncertaintyk[1] = statistics_k2[mincombo_index,2];

bminslopeuncertaintyk[0:4] = b[mincombo_index,0:4]; #finding
#corresponding b values for each minslopeuncertaintyk; this is done before finding
#the kink times because the b array is 48+/-2 (based on site) cells long

posminslopeuncertaintyk = mincombo_index[0]+2; #adding the 3 points
#that were excluded from the first kink, to find positions on whole
#array

kinkpos = x[posminslopeuncertaintyk];#finding
#corresponding time location for min slopeuncertaintyk values

#----------RECALULATING Y VALUES FOR BEST FITTED KINK LINE----------#

#Recreating the D vector to calculate the yfit for the minslopeuncertaintyk
Dminslopeuncertaintyk = np.zeros(nPts);
Dminslopeuncertaintyk[posminslopeuncertaintyk[0]:nPts] = 1;

yfitminslopeuncertaintyk = ((bminslopeuncertaintyk[0]) + (bminslopeuncertaintyk[1])*x + (bminslopeuncertaintyk[2])*Dminslopeuncertaintyk + (bminslopeuncertaintyk[3])*Dminslopeuncertaintyk*x);
yfitkinkline = yfitminslopeuncertaintyk; #saving yfit line with kink for future plotting

#---------Plotting final Slope Break
#plt.plot(x,yfitminslopeuncertaintyk, 'y-')
#plt.title("Displacment for PS Point %i" %(in_ID))


x_direct = x[0:posminslopeuncertaintyk[0]+1]
y_direct = y[0:posminslopeuncertaintyk[0]+1]

x_refract = x[posminslopeuncertaintyk[0]:nPts]
y_refract = y[posminslopeuncertaintyk[0]:nPts]


#------- Computing final slope statistics and labeling direct and refreacted waves
# Compute line statistics for 1st slope, function returns [m, b, m_uncertainty, b_uncertainty,  R2], for final slope break
statistics_s1 = lstats(x_direct,y_direct) 
    
# Compute slope uncertainty for 2nd slope, function returns [m, b, m_uncertainty, b_uncertainty,  R2]
statistics_s2 = lstats(x_refract,y_refract) # want points to cross over the "kink point"

if np.abs(statistics_s1[0]) < np.abs(statistics_s2[0]) :
         print('\n>>> !ERROR! You may have an hidden layer problem, your direct wave has a much faster velocity (higher absolute slope) OR you did not start picking from time (0)!')
      
if statistics_s1[0] and statistics_s2[0] < 0: #if slopes neg they are reverse shot
    rotn = -180
    vanchor = 'baseline' #'top', 'bottom', 'center', 'baseline', 'center_baseline'
    hanchor = 'left' #'center', 'right', 'left'
    
elif statistics_s1[0] and statistics_s2[0] > 0:
    rotn = 0
    vanchor = 'baseline'
    hanchor = 'center'

plt.title('SEG2 First Arrivals Pcks from: %s'%(pckfile),fontweight ="bold")

label_s1 = 'Direct Wave'
label_s2 = 'Refracted Wave'
color1 = (0.1, 0.7, 0.1)
color2 = (0.9, 0.1, 0.1)
m_dir = statistics_s1[0]
m_refract = statistics_s2[0]
b_dir = statistics_s1[1]
b_refract = statistics_s2[1]

#Adding trendline for slope 1
z1 = np.polyfit(x_direct, y_direct, 1)
p1 = np.poly1d(z1)
plt.plot(x_direct,p1(x_direct),c = color1,linewidth=3, label=label_s1)

#Adding trendline for slope 2
z2 = np.polyfit(x_refract,y_refract, 1)
p2 = np.poly1d(z2)
plt.plot(x_refract,p2(x_refract),c = color2,linewidth=3,label=label_s2)

plt.axvline(x= kinkpos,c=(1,.2,1), color = 'r')

plt.legend(fontsize=15)

# =============================================================================
#         #Just to move eqns up and down
#         delta_ydirect = y_direct[3] - y_direct[0]
#         delta_xdirect = x_direct[3] - x_direct[0]
#         delta_yrefract = y_refract[3] - y_refract[0]
#         delta_xrefract = x_refract[3] - x_refract[0]
# =============================================================================

# Tranformations & Locations to plot eqn text
T1_direct = plt.gca().transData.transform_point((x_direct[0], y_direct[0]))
T2_direct = plt.gca().transData.transform_point((x_direct[-1], y_direct[-1]))
dy_direct = (T2_direct[1] - T1_direct[1])
dx_direct = (T2_direct[0] - T1_direct[0])
angle_direct = np.degrees(np.arctan2(dy_direct, dx_direct))
txt1 = np.array((np.median(x_direct), np.median(y_direct)))

T1_refract = plt.gca().transData.transform_point((x_refract[0], y_refract[0]))
T2_refract = plt.gca().transData.transform_point((x_refract[-1], y_refract[-1]))
dy_refract = (T2_refract[1] - T1_refract[1])
dx_refract = (T2_refract[0] - T1_refract[0])
angle_refract = np.degrees(np.arctan2(dy_refract, dx_refract))
txt2 = np.array((np.median(x_refract), np.median(y_refract)))

   # print(delta_ydirect, delta_xrefract, delta_ydirect, delta_xdirect)

#Plot text
plt.text(txt1[0], txt1[1], 'y = %0.3fx + %0.3f'%(m_dir,b_dir), fontsize=16, ha=hanchor, va=vanchor,
               rotation=angle_direct+rotn)

plt.text(txt2[0], txt2[1], 'y = %0.3fx + %0.3f'%(m_refract,b_refract), fontsize=16, ha=hanchor, va=vanchor,
               rotation=angle_refract+rotn)
