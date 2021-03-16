# graphSEG2: a Python based interactive Seismic Refraction Graphing Toolset 

The use of seismic refraction datasets is crucial to Near-surface geophysical explorations and are often used to delineate depths to bedrock. This project will utilize the Obspy python framework to read SEG2 refraction datasets and plot interactive stream datasets using the matplotlib.pyplot.ginput function.


Table of contents
-----------------
- [Installation](#Installation)
- [Importing](#Importing)
- [SEG2grapher.plot](#seg2grapherplot)
- [SEG2grapher.pick](#seg2grapherpick)
  *  [Making First Arrival Picks](#making-first-arrivals-picks)
- [SEG2grapher.analyze](#seg2grapheranalyze)
  *  [Slope Break Method](#slope-break-method)
  *  [Linestatistics.py](#linestatistics)
- [SEG2grapher.crossover](#seg2graphercrossover)
  *  [Findclosestpoint.py](#findclosestpoint)
  *  [Running Findclosestpoint.py](#runningFindclosestpoint.py)


## Installation

1.	Add obspy to separate conda environment:
```bash 
conda create -n obspy python=3.6
```
2.	Initialize the environment using the command while in the anaconda command prompt: 
```bash
conda activate obspy
```
3.	Since this module works best with spyder, while in activated obspy environment install spyder using:
```bash 
conda install spyder 
```
then update all packages using:
```bash 
conda update --all 
```
4.	For missing dependencies e.g. pandas use:
``` bash
conda install pandas
```
6.	To get started, activate spyder by typing spyder in open terminal.


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

## SEG2grapher.pick
SEG2grapher.pick allows the user to pick first arrivals for Seismic Refraction datasets by utilizing the obspy.read SEG2 support function. To execute the SEG2grapher.pick, the user will need to input a filename, gain, nclicks, original name and startpicks. 
SEG2grapher.pick(filename, gain, nclicks, originalname, startpicks)
SEG2grapher.pick(‘300.dat’,15, 100,24,'false', ‘true’)
Where:
1.	Filename (*str*) – should contain a string, e.g., ‘103.dat.’
2.	Gain (*int*) – should be an integer, e.g., 1 for no gain or 10 to amplify the dataset by a factor of 10
3.	Nclicks (*int*) – typically the number of geophones used within the survey, e.g., 24
4.	Original name (*str*) – program asks if the user wants to save the file’s original name after completing the first arrival picks. ‘true’ indicates that the user would like to save picks with the files original name while ‘false’ prompts the user to enter the desirable name after completing picks. The input filename should be filename.txt, with no quotations.
5.	Startpicks (*str*) – user will need to specific if they will be picking first arrivals or viewing the plot by setting startpicking to ‘true’ or ‘false’, where ‘true’ indicates picking first arrivals.

**NB**	Important notice: Zoom tool   in interactive figure window allows user to clearly see first arrivals for picking, however the user risks added an unwanted selection point. If the user does decides to use the zoom tool, the user should click the right mouse after zooming to remove any unwanted point selections.
After running the above code, the user is prompted to read a mini guide and informed to press ENTER to continue or CTRL+C to end.

   ###### Making First Arrival Picks
   First arrival picking in the module utilizes matplotlib.pyplot.ginput for interactive selections. When picking points: 
   * selecting = LEFT MOUSE
   * removing point = RIGHT MOUSE
   *	end selection = ENTER

**NB** Note that the user defined Nclicks is the maximum number of selections the user is allowed but one can exit the selection mode by simply pressing ENTER.

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
SEG2grapher.crossover(filename)
`````
````
SEG2grapher.crossover (‘300.txt’)
````
Where:
1.	Filename (*str*) – should contain a string, e.g., ‘103.txt.’ This text file should have x,y dataset space or tab-separated. Where the x column is the distance, and the y column is the time.

