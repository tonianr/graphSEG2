# graphSEG2: a Python based interactive Seismic Refraction Graphing Toolset. 

The use of seismic refraction datasets is crucial to Near-surface geophysical explorations and are often used to delineate depths to bedrock. This project will utilize the Obspy python framework to read SEG2 refraction datasets and plot interactive stream datasets using the matplotlib.pyplot.ginput function.


Table of contents
-----------------
- [Installation](#Installation)
- [Importing](#Importing-function)
- [graphSEG2.pick](#graphSEG2.pick)
  *  [Making First Arrival Picks](#first-arrivals)
- [graphSEG2.analyze](#graphSEG2.analyze)
  *  [Slope Break Method](#slope-break)
  *  [Linestatistics.py](#linestatistics)
- [graphSEG2.crossover](#graphSEG2.crossover)
  *  [Findclosestpoint.py](#Findclosestpoint)
  *  [Running Findclosestpoint.py](#RunningFindclosestpoint.py)


## Installation

- (1)	Add obspy to separate conda environment:
```bash 
conda create -n obspy python=3.6.
```
ii.	Initialize the environment using the command while in the anaconda command prompt: conda activate obspy.
iii.	Since this module works best with spyder in activated obspy environment install spyder using conda install spyder then conda update --all. 
iv.	For missing dependencies e.g. pandas use conda install pandas. 
v.	To get started, activate spyder by typing spyder in open terminal.
