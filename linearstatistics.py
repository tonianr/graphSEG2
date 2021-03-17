# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 19:25:25 2020

@author: tonian
"""
import sys
import numpy as np


def linestatistics(x,y):
       # print('This functions using input (x,y) datasets to compute various statistical parameters for a line.')
       # print('It returns [m, b, m_uncertainty, b_uncertainty,  R2]')

        # STATISTICS
        sx = sum(x);
        sy = sum(y);
        sxx = sum(x * x);
        syy = sum(y * y);
        sxy = sum(x * y);
        sx2 = sx * sx;
        sy2 = sy * sy;
    
        n=len(x);
    
        m = (n * sxy - sx * sy)/(n * sxx - sx2); #slope
    
        b = (1/n) * (sy - m * sx); #intercept
    
        s = np.sqrt((syy - sy2/n - m * (sxy - sx * sy/n))/(n-2));
    
        m_uncertainty = s * np.sqrt(n / (n * sxx - sx2)); #uncertainty in slope
    
        b_uncertainty = s * np.sqrt(sxx / (n * sxx - sx2)); #uncertainty in intercept
    
        R2 = ((sxy - sx * sy /n) / np.sqrt((sxx - sx2/n) * (syy - sy2/n)))**2; #R^2
        
        statistics_out = [m, b, m_uncertainty, b_uncertainty,  R2]
        
    # =============================================================================
    #     #Plotting 
    #     plt.plot(x,y, 'b-o')
    #     #plt.title("Displacment for PS Point %i" %(in_ID))
    #     plt.ylabel('Displacement (mm)')
    #     plt.xlabel('Time (yr)')
    #     
    #     #Adding trendline
    #     z = np.polyfit(x, y, 1)
    #     p = np.poly1d(z)
    #     plt.plot(x,p(x),"r--")
    # =============================================================================
        return statistics_out
    
 #To make this script callable as a function add the line below and do a "practice run" code below   
if __name__ == "__main__":
   print("File %s executed as a program, so run some tests on the functions"  % sys.argv[0])
   x = np.array([2,7,8,9])
   y = np.array([20, 93, 120, 349])
   stats = linestatistics(x, y)
   print(stats)
