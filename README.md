# graphSEG2: a Python based interactive Seismic Refraction Graphing Toolset 

Seismic refraction datasets are crucial to Near-surface geophysical explorations and are often used to delineate depths to bedrock. This Python 3 module utilizes the Obspy framework to read SEG2 refraction datasets and plot interactive stream datasets using the matplotlib.pyplot.ginput function. It was created as a simple way for students to visualize and analyze SEG2 datasets.

Note: The user is advised that this toolset was created with SEG2 files from Geometrics, and they may need to adjust some data extract parameters within the graphSEG2.plot and graphSEG2.pick sections for use. According to obspy "UserWarning: Many companies use custom defined SEG2 header variables. This might cause basic header information reflected in the single traces' stats to be wrong (e.g. recording delays, first sample number, station code names, ..). Please check the complete list of additional unmapped header fields that gets stored in Trace.stats.seg2 and/or the manual of the source of the SEG2 files for fields that might influence e.g. trace start times."

Table of contents
-----------------
- [Recommended ObsPy Installation](#Recommended-ObsPy-Installation)
- [Importing](#Importing)
- [SEG2grapher.plot](#seg2grapherplot)
- [SEG2grapher.pick](#seg2grapherpick)
  * [Making First Arrival Picks](#making-first-arrival-picks)
- [SEG2grapher.analyze](#seg2grapheranalyze)
- [SEG2grapher.crossover](#seg2graphercrossover)
- [Other Information](#other-information)
- [Possible Hiccups](#possible-hiccups)
  * [non-GUI Backend Error](#non-GUI-Backend-Error)

## Recommended ObsPy Installation

1.	Add obspy to separate conda environment:
```bash 
conda create -n obspy python=3.8
```
2.	Initialize the environment using the command while in the anaconda command prompt: 
```bash
conda activate obspy
```
3.	Install obspy package: 
```bash
conda install -c conda-forge obspy
```
4.	Since this module works best with spyder, while in activated obspy environment install spyder using:
```bash 
conda install -c conda-forge spyder
```
5.	For missing dependencies e.g. pandas use:
``` bash
conda install pandas
```
6.	To get started, activate spyder by typing spyder in open terminal.
``` bash
spyder
```


## Importing

To import the function, use the code:`from graphSEG2_v2 import graphSEG2` in python script.
This allows the user script in use to have access to the graphSEG2 function and its executables. While in spyder highlight graphSEG2 in the script and press *CTRL + i* for a mini user guide.

## SEG2grapher.plot
SEG2grapher.plot allows the user to plot Seismic Refraction datasets by utilizing the obspy.read SEG2 support function. To execute the SEG2grapher.plot, the user will need to input a filename, gain. The user will also use this function to make estimates on the best gain value to use for the input file.
```bash
SEG2grapher.plot(filename, gain)
```
```bash
SEG2grapher.plot(‘304.dat’,15)
```
Where:
1.	Filename (*str*) – should contain a string, e.g., ‘103.dat.’
2.	Gain (*int*) – should be an integer, e.g., 1 for no gain or 10 to amplify the dataset by a factor of 10

![image](https://user-images.githubusercontent.com/78233913/151029942-4e87d072-44fe-4dc8-9fe8-33e4557505b4.png)



## SEG2grapher.pick
SEG2grapher.pick allows the user to pick first arrivals for Seismic Refraction datasets by utilizing the obspy.read SEG2 support function. To execute the SEG2grapher.pick, the user will need to input a filename, gain, nclicks, original name and startpicks. 
```bash
SEG2grapher.pick(filename, gain, originalname)
```
```bash
SEG2grapher.pick(‘300.dat’,15,'false')
```
Where:
1.	Filename (*str*) – should contain a string, e.g., ‘103.dat.’
2.	Gain (*int*) – should be an integer, e.g., 1 for no gain or 10 to amplify the dataset by a factor of 10
3.	Original name (*str*) – program asks if the user wants to save the file’s original name after completing the first arrival picks. ‘true’ indicates that the user would like to save picks with the files original name while ‘false’ prompts in the python console the user to enter the desirable name after completing picks. The input filename should be filename.txt, with no quotations.

**NB**	Important notice: Zoom tool in the interactive figure window allows user to clearly see first arrivals for picking, however the user risks added an unwanted selection point. If the user does decides to use the zoom tool, the user should click the RIGHT-MOUSE after zooming to remove any unwanted point selections.
After running the above code, the user is prompted to read a mini guide and informed to press ENTER to continue or CTRL+C to end.

#### Making First Arrival Picks
First arrival picking in the module utilizes matplotlib.pyplot.ginput for interactive selections. When picking points: 
   * selecting = LEFT MOUSE
   * removing point = RIGHT MOUSE or CENTER MOUSE
   *	end selection = ENTER

**NB** Note that the user defined the maximum number of selections the user allowed is 96 (within the max time of 5 minutes) but one can exit the selection mode by simply pressing ENTER or CENTER MOUSE. 

## SEG2grapher.analyze
SEG2grapher.analyze allows the user to plot and analyze Seismic Refraction first arrival picks using a Slope Break method to find possible kinks in the plotted first arrivals that would represent the crossover point. To execute the SEG2grapher.analyze, the user will need to input a pick file and show_autocrossover.
```bash
SEG2grapher.analyze(pckfilename, show_autocrossover)
```
````bash
SEG2grapher.analyze(‘300.txt’,'false')
````
Where:
1.	Pckfilename (*str*) – should contain a string of the text file, e.g., ‘103.txt.’ This text file should have x,y dataset space or tab-separated. Where the x column is the distance, and the y column is the time.
2.	Show_autocrossover (*str*) – indicates whether the user would like to plot all possible crossover points on the analysis plot; ‘true’ shows slope breaks while ‘false’ plots only the analysis and final crossover point. 

## SEG2grapher.crossover
SEG2grapher.crossover allows the user to plot Seismic Refraction first arrival picks from text file and select desired crossover point. To execute the SEG2grapher.crossover, the user will need to input a filename.
`````bash
SEG2grapher.crossover(filename, Show_autocrossover )
`````
````
SEG2grapher.crossover (‘300.txt’,'true')
````
Where:
1.	Filename (*str*) – should contain a string, e.g., ‘103.txt.’ This text file should have x,y dataset space or tab-separated. Where the x column is the distance, and the y column is the time.
2.	Show_autocrossover (*str*) – indicates whether the user would like to plot all possible crossover points on the analysis plot; ‘true’ shows slope breaks while ‘false’ plots only the analysis and final crossover point.

## Other Information
**Loading file from another folder**
When loading file from folder make sure to source folder.
Example:
```
SEG2grapher.pick(‘./foldername/300.dat’,15, 'false')
```
OR
```
SEG2grapher.analyze('./foldername/102.txt','true')
```
OR
```
SEG2grapher.crossover('./foldername/102.txt','true')
```

**Saving pick txt file to specified folder**
When prompted to input output file name after picking first arrivals (SEG2grapher.pick original name set to ‘false’) user can add specific file location.
Example when prompted in console:
```
Write the output file's name as filename.txt: ./foldername/103.txt
```
## Possible Hiccups
#### non-GUI Backend Error
Error message in console: *Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.* 
OR 
There is no figure window pop-up after running the module.

Fix by changing IPython Console graphics backend to Qt5 in the Preferences window. Source: https://stackoverflow.com/questions/36700083/in-spyder-plot-using-matplotlib-with-interactive-zoom-etc	


