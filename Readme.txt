
JeongYeon Han and Wonsuck Kim
11/09/21

This folder contains code to run numerical simulations for the manuscript, 'Linking levee-building process with river avulsion: geomorphic analysis for assessing avulsion frequency and style'.
Data was produced by running this code repeatedly, changing the input overflow velocity and median grainsize parameters as described in the manuscript. Below is a brief outline of what each file does in this folder.

1. LBM.m 
This is the Test1, 4, and 5 script, responsible for running repeated levee aggradation under a given set of input conditions. (Title is shorthand for "Levee Building Model")

2. LBM_T2.m 
This is the Test2 script, responsible for running repeated levee aggradation under a given set of input conditions. (Title is shorthand for "Levee Building Model Test2")

3. LBM_T3.m 
This is the Test3 script, responsible for running repeated levee aggradation under a given set of input conditions. (Title is shorthand for "Levee Building Model Test3")

4. LvSlopeAvFreq.m 
This is the levee slope and avulsion frequency script, responsible for running repeated levee aggradation with ranges of input overflow and suspended sediment conditions. (Title is shorthand for "Levee Slope and Avulsion Frequency")
Here, for diversing colour pattern for plotting Figure 7 contour maps, download multigradient.m file from https://github.com/lrkrol/multigradient, version 1.5.6 (46.3 KB) by Laurens R Krol if it is necessary. 
