# -*- coding: utf-8 -*-
"""
=====================================================
SEG2grapher: Mini User Guide
=====================================================
This module focuses on using the obspy framework for creating an interactive seismic refraction plot.
\n Created on Mon Nov 30 15:28:30 2020
\n @author: Tonian Robinson; tonianr@usf.edu

\n=====================================================
\n`IMPORT`::
    >>> from graphSEG2 import SEG2grapher

\n=====================================================
\n`PLOT` 
    * Plots data for viewing."
    * Use this step for understanding the gain best for your picking.
    * Format for input::
        >>> SEG2grapher.pick(infile,gain).
    * e.g.
        >>> SEG2grapher.plot('305.dat',6,)
           
\n===================================================
\n`PICK`
    * Plots data for first arrival point picking.
    * Start picks from the (0,0) point on your plot.
    * Format for input:
        >>> SEG2grapher.pick(infile,gain, originalname).
    * Picking:
        >>> selecting = LEFT MOUSE
        >>> removing point = RIGHT MOUSE
        >>> end selection = ENTER or CENTER MOUSE
    * e.g.
        >>> SEG2grapher.plot('305.dat',6,'true')
  
\n====================================================
\n`ANALYZE`
    * Plots pick data and automatically finds crossover point.
    * SEG2grapher.analyze requires a simple x, y dataset (no headers), where the distance is in the 1st collumn and the time in the 2nd collumn.
    * SEG2grapher.analyze uses a slope break method to find the the intercept point of the direct and refracted wave arrivals, with the assumption that the user has experience picking arrivals for refraction datasets.
    * Format for input:
        >>> SEG2grapher.analyze(pckfile,show_slopebreaks)         
    * e.g. 
        >>> SEG2grapher.analyze('304.txt','true')
    * NB: You can reshape plot window to allow text annotations to align with slopes. Text annotations align best when plot windows are square.
\n====================================================
\n`CROSSOVER`
    * Plots pick data for user to manually select crossover point.
    * SEG2grapher.crossover requires a simple x, y dataset (no headers), where the distance is in the 1st collumn and the time in the 2nd collumn.
    * In SEG2grapher.crossover the user is allowed to select desired crossover point on plotted refracted dataset.
    * Format for input:
        >>> SEG2grapher.crossover(pckfile,show_slopebreaks)     
    * e.g. 
        >>> SEG2grapher.crossover('304.txt')   
    * NB: You can reshape plot window to allow text annotations to align with slopes. Text annotations align best when plot windows are square.
\n====================================================
                 

"""

class SEG2grapher():
    
    @staticmethod
    def info():
        infomsg = "\n>> Welcome graphSEG2 user! For a mini user guide highlight graphSEG2 and select CTRL + i. \n"\

                
        return infomsg
    
    @staticmethod
    def plot(infile,gain):
        from matplotlib import pyplot as plt
        import obspy
        import numpy as np
        

        print(SEG2grapher.info()) # print simple user guide information 
        
        
        #Import SEG2 using Obspy
        obspySEG2_input = obspy.read(infile, format="seg2")
        
        geophone1st_pos = np.abs(float(obspySEG2_input[0].stats.seg2.RECEIVER_LOCATION)) #position of 1st geophone
        shot_location = float(obspySEG2_input[0].stats.seg2.SOURCE_LOCATION)
        geophone_spacing =  np.abs(float(obspySEG2_input[0].stats.seg2.RECEIVER_LOCATION) - float(obspySEG2_input[1].stats.seg2.RECEIVER_LOCATION)) #geophone interval/ spacing
        record_start = 0
        record_end = obspySEG2_input[0].stats.delta*len(obspySEG2_input[0])*1000
        
        for i in range(len(obspySEG2_input)):
                
                dist = (obspySEG2_input[i].data)/max(obspySEG2_input[i].data)*gain+geophone1st_pos+geophone_spacing*i
                time = ([obspySEG2_input[i].stats.delta*k*1000 for k in range(len(obspySEG2_input[i]))])
                
                #-----Plot traces
                trace, = plt.plot(dist,time, color=(0.2, 0.1, 0.3)) 
                #----- Fill traces
                plt.fill_betweenx(trace.get_ydata(),geophone1st_pos+geophone_spacing*i,trace.get_xdata(),
                                        where = trace.get_xdata() >= geophone1st_pos+i*geophone_spacing,color='black')
                #-----Clip traces
                trace.get_xdata()[trace.get_xdata() < geophone1st_pos+i*geophone_spacing-((geophone_spacing)*1.5)] = geophone1st_pos+i*geophone_spacing-((geophone_spacing)*1.5)
                trace.get_xdata()[trace.get_xdata() > geophone1st_pos+i*geophone_spacing+((geophone_spacing)*1.5)] = geophone1st_pos+i*geophone_spacing+((geophone_spacing)*1.5)
                trace.set_xdata(trace.get_xdata())
                
                plt.title('SEG2 file: %s.dat'%(obspySEG2_input[0].stats.seg2.SHOT_SEQUENCE_NUMBER),fontweight ="bold",fontsize=20)
                #plt.ylim(0,ymax)
                
        plt.xlabel('Distance (m)',fontsize=20)
        plt.xticks(fontsize=18)
        plt.ylabel('Time (ms)',fontsize=20)
        plt.yticks(fontsize=18)
        plt.ylim(record_start-5,record_end)
        plt.xlim(geophone1st_pos-geophone_spacing,geophone1st_pos+geophone_spacing*len(obspySEG2_input))
        plt.gca().invert_yaxis()  

        plt.text(max(dist)-5, -8, 'shot location: %0.1f m'%(shot_location), fontsize=18, ha='right')
        # mng = plt.get_current_fig_manager()
        # mng.window.showMaximized()
        plt.show() 
   
        return 
    
    @staticmethod
    def pick(infile,gain, originalname):
        import os
        from matplotlib import pyplot as plt
        import numpy as np
        import obspy

        print(SEG2grapher.info()) # print simple user guide information 
        
        input("\n<<<Press ENTER to continue OR CTRL+C to end >>>")  # add press enter to continue?
        
        #Import SEG2 using Obspy
        obspySEG2_input = obspy.read(infile, format="seg2")
        
        nclicks = 96 #number of clicks per image
        
        geophone1st_pos = np.abs(float(obspySEG2_input[0].stats.seg2.RECEIVER_LOCATION)) #position of 1st geophone
        shot_location = float(obspySEG2_input[0].stats.seg2.SOURCE_LOCATION)
        geophone_spacing =  np.abs(float(obspySEG2_input[0].stats.seg2.RECEIVER_LOCATION) - float(obspySEG2_input[1].stats.seg2.RECEIVER_LOCATION)) #geophone interval/ spacing
        record_start = 0
        record_end = obspySEG2_input[0].stats.delta*len(obspySEG2_input[0])*1000
        
        for i in range(len(obspySEG2_input)):
                
                dist = (obspySEG2_input[i].data)/max(obspySEG2_input[i].data)*gain+geophone1st_pos+geophone_spacing*i
                time = ([obspySEG2_input[i].stats.delta*k*1000 for k in range(len(obspySEG2_input[i]))])
                
                #-----Plot traces
                trace, = plt.plot(dist,time, color=(0.2, 0.1, 0.3)) 
                #----- Fill traces
                plt.fill_betweenx(trace.get_ydata(),geophone1st_pos+geophone_spacing*i,trace.get_xdata(),
                                        where = trace.get_xdata() >= geophone1st_pos+i*geophone_spacing,color='black')
                #-----Clip traces
                trace.get_xdata()[trace.get_xdata() < geophone1st_pos+i*geophone_spacing-((geophone_spacing)*1.5)] = geophone1st_pos+i*geophone_spacing-((geophone_spacing)*1.5)
                trace.get_xdata()[trace.get_xdata() > geophone1st_pos+i*geophone_spacing+((geophone_spacing)*1.5)] = geophone1st_pos+i*geophone_spacing+((geophone_spacing)*1.5)
                trace.set_xdata(trace.get_xdata())
                
                plt.title('SEG2 file: %s.dat'%(obspySEG2_input[0].stats.seg2.SHOT_SEQUENCE_NUMBER),fontweight ="bold",fontsize=20)
                #plt.ylim(0,ymax)
                
        plt.xlabel('Distance (m)',fontsize=20)
        plt.xticks(fontsize=18)
        plt.ylabel('Time (ms)',fontsize=20)
        plt.yticks(fontsize=18)
        plt.ylim(record_start-5,record_end)
        plt.xlim(geophone1st_pos-geophone_spacing,geophone1st_pos+geophone_spacing*len(obspySEG2_input))
        plt.gca().invert_yaxis()  
        plt.text(max(dist)-5, -8, 'shot location: %0.1f m'%(shot_location), fontsize=18, ha='right')

        # mng = plt.get_current_fig_manager()
        # mng.window.showMaximized()
 

        pcks = plt.ginput(n=nclicks, timeout=300, show_clicks=True) #using ginput to select values from plot
        print('\n>> %d picks made, output pick values:' %(len(pcks)))
        print(pcks)
        
        if originalname == 'true':
            outpckfilename = os.path.splitext(infile)[0] #to keep original in file name
            outpckfilename = outpckfilename+'.txt'
            print(outpckfilename)
        elif originalname == 'false':
            outpckfilename = input("Write the output file's name as filename.txt OR file location as ./foldername/filename.txt:")
        np.savetxt(outpckfilename, pcks)


        plt.show() 

        return 
    
    @staticmethod
    def analyze(pckfile,show_autocrossover):
        import pandas as pd
        import numpy as np
        from matplotlib import pyplot as plt
        from linearstatistics import linestatistics as lstats
        import os
        
        print(SEG2grapher.info()) # print simple user guide information
        
        pckdata = pd.read_csv(pckfile, sep="\s", header=None, engine='python')
        
        pckdata.columns = ["dist", "time"]
        
        x = np.array(pckdata.dist)
        y = np.array(pckdata.time)
        y = y-np.min(y) #making directwave start at (0,0)
        x = x-np.min(x) #making directwave start at (0,0)
        nPts = len(x) #number of Points in each time series
        
        # Compute overall stats for data, function returns [m, b, m_uncertainty, b_uncertainty,  R2]
        initial_stats = lstats(x,y) 

        if initial_stats[0] < 0: #if slopes neg they are reverse shot
            #Text rotation based on reverse vs forward shot
            vanchor = 'baseline' #'top', 'bottom', 'center', 'baseline', 'center_baseline'
            hanchor = 'left' #'center', 'right', 'left'
            print("\n>> Overall slope is (-) this dataset is likely a reverse shot.")
        elif initial_stats[0] > 0:
            #Text rotation based on reverse vs forward shot
            vanchor = 'baseline'
            hanchor = 'right'
        
        
        
        #Plotting 
        plt.figure() #remove if you want points to be plotted over wiggle window immage
        plt.plot(x,y, 'b-o', label='1st Arrival Picks', markersize=8)
        plt.ylabel('Time (ms)',fontsize=20)
        plt.yticks(fontsize=18)
        plt.xlabel('Distance (m)',fontsize=20)
        plt.xticks(fontsize=18)
        plt.title('SEG2 First Arrivals Picks from: %s'%(os.path.basename(pckfile)),fontweight ="bold",fontsize=20)
        plt.xlim(-0.1,max(x))
        plt.ylim(-0.1,max(y))
        
        # Define range of "kinks" to test for a slope break in best fit line
        firstkink = 3 #allowing the kinks to start after the first 3 data Points (remember python indices start at 1 so the first kink is at #4); this is the location of the kink as it relates to the extracted 56 data Pointss
        lastkink = nPts-3 #:nPts; #leaving 4 data points for the last kink
        kinks = np.arange(firstkink, lastkink+1)
        nKinks = len(kinks)
        
        #------------CREATING ARRAYS TO STORE INFORMATION----------------------#
        b = np.zeros((nKinks,4)) #stores best-fitting parameters from 2-line fit for all positions
        yfitk = np.zeros((nPts,1)) #empty array for all kinked lines
        #yfitminslopeuncertaintyk = np.zeros(nKinks); #empty array for y values of minslopeuncertaintyk (shows optimum kink point)
        minslopeuncertaintyk = np.zeros(2) #creating empty array to store max Rsquared values for kinked lines
        kinkpos = np.zeros(1) #creating empty array to store corresponding time values for max Rsquared values
        #yfitkinkline = np.zeros((1,nPts)); #empty array for saving y fit values caluclated with kink location (basically plots kink line) 
        bminslopeuncertaintyk = np.zeros(4) #b values that correspond to the min slopeuncertaintyk value #b(1) intercept of line 1 for nth position
        #b(2)#slope of line 1, b(3)intercept of second line added past kink Points, b(4) slope of second line added past kink Points
        A = np.zeros((nPts,4))
        D = []
        
        #--------Best Fit Two Lines With Kink------------------------
        bk = np.zeros((1,4)) # make bk array to store best-fit coefficients for each kink time
        statistics_k1 = np.zeros((nKinks,5)) # make statistics array to for each kink time for 1st slope
        statistics_k2 = np.zeros((nKinks,5)) # make statistics array to for each kink time for 2nd slope
        
        
        #-----------Begin finding slope breaks------------------------------------
        nk = -1 # make counter for kinks
        for ik in range(firstkink,lastkink+1):
            nk = nk+1
            D = np.zeros(nPts)
            D[ik:nPts+1] = 1  # for Pointss with D = 1 second line is added; letting the first 3 points = 0
            
            # build A array , for which y = A * b and b = A\y 
            A[:,0] = np.ones(nPts)
            A[:,1] = x
            A[:,2] = D
            A[:,3] = np.multiply(D,x)
            
            bk = np.linalg.lstsq(A,y,rcond=None)[0] #calculating b array using the linear algebra - least squares from numpy (has to be least squares becasue the A matrix is not a perfect square)
            #input('Press ENTER to continue');
            
            #store b (array of 4) into the bk row for this kink time
            b[nk,:] = bk #using the corresponding b array to calculate the yfit value for each Points
            
            #compute the yfit best-fitting line with kink through the data
            yfitk = bk[0] + bk[1]*x + bk[2]*D + bk[3]*D*x

            # Compute line statistics for 1st slope, function returns [m, b, m_uncertainty, b_uncertainty,  R2]
            statistics_k1[nk,:] = lstats(x[0:ik],y[0:ik]) 
            
            # Compute slope uncertainty for 2nd slope, function returns [m, b, m_uncertainty, b_uncertainty,  R2]
            statistics_k2[nk,:] = lstats(x[ik-1:nPts],y[ik-1:nPts]) # want points to cross over the "kink point"
        

            # -------Plot kinks/slope breaks------
            # #random colors for fun
            # rc = np.random.randint(1,9)/10
            # bc = np.random.randint(1,9)/10
            # gc = np.random.randint(1,9)/10
            # color = (rc, gc, bc)
            # plt.plot(x,yfitk, c=color, linewidth=1,)
    
        
        #-------------end of for loop for for finding slope breaks---------------------    
        
            
        #Finding minimum combination of slope R2 for all kink lines
        sq_combo = np.square(statistics_k1[:,2]) + np.square(statistics_k2[:,2])
        [mincombo_min,mincombo_index] = [np.amin(sq_combo),np.where(sq_combo == sq_combo.min())]
        #np.where(a == a.min())
        
        #Storing minimum slope uncertainties
        minslopeuncertaintyk[0] = statistics_k1[mincombo_index,2]
        minslopeuncertaintyk[1] = statistics_k2[mincombo_index,2]
        
        bminslopeuncertaintyk[0:4] = b[mincombo_index,0:4] #finding
        #corresponding b values for each minslopeuncertaintyk; this is done before finding
        #the kink times because the b array is 48+/-2 (based on site) cells long
        
        posminslopeuncertaintyk = mincombo_index[0]+2 #adding the 3 points
        #that were excluded from the first kink, to find positions on whole
        #array
        
        kinkpos = x[posminslopeuncertaintyk]#finding
        #corresponding time location for min slopeuncertaintyk values
        
        #----------RECALULATING Y VALUES FOR BEST FITTED KINK LINE----------#
        
        #Recreating the D vector to calculate the yfit for the minslopeuncertaintyk
        Dminslopeuncertaintyk = np.zeros(nPts)
        Dminslopeuncertaintyk[posminslopeuncertaintyk[0]:nPts] = 1
        
        #yfitminslopeuncertaintyk = ((bminslopeuncertaintyk[0]) + (bminslopeuncertaintyk[1])*x + (bminslopeuncertaintyk[2])*Dminslopeuncertaintyk + (bminslopeuncertaintyk[3])*Dminslopeuncertaintyk*x);
        #yfitkinkline = yfitminslopeuncertaintyk; #saving yfit line with kink for future plotting
        
        #---------Plotting final Slope Break
        #plt.plot(x,yfitminslopeuncertaintyk, 'y-')
        #plt.title("Displacment for PS Point %i" %(in_ID))
    
    
        #---------Separating Direct and Refracted Waves
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
    

        color1 = (0.1, 0.7, 0.1)
        color2 = (0.9, 0.1, 0.1)
        m_dir = statistics_s1[0]
        m_refract = statistics_s2[0]
        b_dir = statistics_s1[1]
        b_refract = statistics_s2[1]
        label_s1 = ('Direct Wave [y = %0.2fx + %0.2f]'%(m_dir,b_dir))
        label_s2 = ('Refracted Wave [y = %0.2fx + %0.2f]'%(m_refract,b_refract))
        
        #Adding trendline for slope 1
        z1 = np.polyfit(x_direct, y_direct, 1)
        p1 = np.poly1d(z1)
        plt.plot(x_direct,p1(x_direct),c = color1,linewidth=3, label=label_s1)
        
        #Adding trendline for slope 2
        z2 = np.polyfit(x_refract,y_refract, 1)
        p2 = np.poly1d(z2)
        plt.plot(x_refract,p2(x_refract),c = color2,linewidth=3,label=label_s2)
        
        
        #plotting Crossover point
        if show_autocrossover == 'true':
            #random colors for fun
            rc = (np.random.randint(1,9))/10
            bc = (np.random.randint(1,9))/10
            gc = (np.random.randint(1,9))/10
            color = [rc, bc, gc]   
            label_s3 = ('Crossover [%0.2f m]'%kinkpos)
            #plt.axvline(x= kinkpos,c=(1,.2,1), color = 'r',label=label_s3)
            plt.axvline(x= kinkpos, c=color,label=label_s3)
        elif show_autocrossover == 'false':
            pass  
        print('\n >> Crossover point (x,y): %s'%kinkpos)
 
        plt.legend(fontsize=18)
        
        # # Tranformations & Locations to plot eqn text
        # txt1 = np.array((np.median(x_direct), np.median(y_direct))) 
        # txt2 = np.array((np.median(x_refract), np.median(y_refract)))

        # #Plot equation text
        # plt.text(txt1[0], txt1[1], 'y = %0.3fx + %0.3f'%(m_dir,b_dir), fontsize=18, ha=hanchor, va=vanchor)
        
        # plt.text(txt2[0], txt2[1], 'y = %0.3fx + %0.3f'%(m_refract,b_refract), fontsize=18, ha=hanchor, va=vanchor)
        plt.show()
        
        return
    
    @staticmethod
    def crossover(pckfile, show_autocrossover): #pick cross over point
        import pandas as pd
        import numpy as np
        from matplotlib import pyplot as plt
        from linearstatistics import linestatistics as lstats
        from findclosestpoint import closestpoint
        import os
        
        print(SEG2grapher.info()) # print simple user guide information
        
        pckdata = pd.read_csv(pckfile, sep="\s", header=None, engine='python')
        
        pckdata.columns = ["dist", "time"]
        
        x = np.array(pckdata.dist)
        y = np.array(pckdata.time)
        y = y-np.min(y) #making directwave start at (0,0)
        x = x-np.min(x) #making directwave start at (0,0)
        
        # Compute overall stats for data, function returns [m, b, m_uncertainty, b_uncertainty,  R2]
        initial_stats = lstats(x,y) 

        if initial_stats[0] < 0: #if slopes neg they are reverse shot
            #y = y-np.min(y) #making directwave start at (0,0) for reverse 
            #Text rotation based on reverse vs forward shot
            vanchor = 'baseline' #'top', 'bottom', 'center', 'baseline', 'center_baseline'
            hanchor = 'left' #'center', 'right', 'left'
            print("\n>> Overall slope is (-) this dataset is likely a reverse shot.")
        elif initial_stats[0] > 0:
            x = x-np.min(x) #making directwave start at (0,0) for forward
            #Text rotation based on reverse vs forward shot
            vanchor = 'baseline'
            hanchor = 'right'
        
        
        nPts = len(x) #number of Points in each time series
        
        #---------Plotting Input
        plt.figure()
        plt.title('SEG2 First Arrivals Picks from: %s'%(os.path.basename(pckfile)),fontweight ="bold", fontsize=20)
        plt.plot(x,y, 'b-o', label='1st Arrival Picks', markersize=8)
        plt.ylabel('Time (ms)',fontsize=20)
        plt.yticks(fontsize=18)
        plt.xlabel('Distance (m)',fontsize=20)
        plt.xticks(fontsize=18)
        plt.xlim(-0.1,max(x))
        plt.ylim(-0.1,max(y))
        
        
        pck = plt.ginput(n=1, timeout=30, show_clicks=True) #using ginput to select desired crossover point
        print('\n >> Crossover point (x,y): %s'%pck)
        
        #---------Finding Position of Crossover point
        pck = np.array(pck,ndmin=2)
        pt = np.stack(pck)
        pts = np.stack((x,y))
        pts = pts.T
        
        pospck = closestpoint(pt, pts)
        #print(pospck)
        
        #---------Separating Direct and Refracted Waves

        x_direct = x[0:pospck+1]
        y_direct = y[0:pospck+1]
        
        x_refract = x[pospck:nPts]
        y_refract = y[pospck:nPts]

        #------- Computing final slope statistics and labeling direct and refreacted waves
        # Compute line statistics for 1st slope, function returns [m, b, m_uncertainty, b_uncertainty,  R2], for final slope break
        statistics_s1 = lstats(x_direct,y_direct) 
            
        # Compute slope uncertainty for 2nd slope, function returns [m, b, m_uncertainty, b_uncertainty,  R2]
        statistics_s2 = lstats(x_refract,y_refract) # want points to cross over the "kink point"
        
        if np.abs(statistics_s1[0]) < np.abs(statistics_s2[0]) :
                 print('\n>>> !ERROR! You may have an hidden layer problem, your direct wave has a much faster velocity (higher absolute slope) OR you did not start picking from time (0)!')
        

       
        color1 = (0.1, 0.7, 0.1)
        color2 = (0.9, 0.1, 0.1)
        m_dir = statistics_s1[0]
        m_refract = statistics_s2[0]
        b_dir = statistics_s1[1]
        b_refract = statistics_s2[1]
        label_s1 = ('Direct Wave [y = %0.2fx + %0.2f]'%(m_dir,b_dir))
        label_s2 = ('Refracted Wave [y = %0.2fx + %0.2f]'%(m_refract,b_refract))
        
        #Adding trendline for slope 1
        z1 = np.polyfit(x_direct, y_direct, 1)
        p1 = np.poly1d(z1)
        plt.plot(x_direct,p1(x_direct),c = color1,linewidth=3, label=label_s1)
        
        
        #Adding trendline for slope 2
        z2 = np.polyfit(x_refract,y_refract, 1)
        p2 = np.poly1d(z2)
        plt.plot(x_refract,p2(x_refract),c = color2,linewidth=3,label=label_s2)
        
        #plotting Crossover point
        if show_autocrossover == 'true':
            #random colors for fun
            rc = (np.random.randint(1,9))/10
            bc = (np.random.randint(1,9))/10
            gc = (np.random.randint(1,9))/10
            color = [rc, bc, gc]   
            label_s3 = ('Crossover [%0.2f m]'%x[pospck])
            #plt.axvline(x= x[pospck],c=(1,.2,1), color = 'r',label=label_s3)
            plt.axvline(x= x[pospck], c=color,label=label_s3)
        elif show_autocrossover == 'false':
            pass  

        plt.legend(fontsize=18)

        
        # # Tranformations & Locations to plot eqn text
        # txt1 = np.array((np.median(x_direct), np.median(y_direct))) 
        # txt2 = np.array((np.median(x_refract), np.median(y_refract)))

        # #Plot equation text
        # plt.text(txt1[0], txt1[1], 'y = %0.3fx + %0.3f'%(m_dir,b_dir), fontsize=18, ha=hanchor, va=vanchor)
        
        # plt.text(txt2[0], txt2[1], 'y = %0.3fx + %0.3f'%(m_refract,b_refract), fontsize=18, ha=hanchor, va=vanchor)

        plt.show()
        return
        
        
#To make this script callable as a function add the line below and do a "practice run" code below   
if __name__ == "__main__":
    #SEG2grapher.plot("206.dat",6, 100,3,'true')
    print(SEG2grapher.info()) # print simple user guide information 
    SEG2grapher.analyze('./exampledata/206.txt','false')

