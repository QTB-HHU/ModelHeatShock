# ModelHeatShock

## Credits

The python code present in this directory has been written by Stefano Magni and it is based on scientific work performed by Stefano Magni, Antonella Succurro, Alexander Skupin and Oliver Ebenhoeh, which can be found in the following paper: XXX. This code can be used according to what specifyed in the license file that you can find in this folder.


## What is it?

This code represent the implementation of our mode lon the Heat Shock Response in Chlamydomonas reinhardtii which we used to perform all the simulations (and additional ones) described in the above mentioned paper, which we suggest to read in order to understand the scientific content of this code.


## How to run the code

To run our code:

1. Download the full content of the directory containing this file.
2. Install the requirements listed in the corresponding file "requirements.txt" (this step might be unnecessary if you already have installed on your system all necessary python libraries).
3. open a terminal and type "python3 Model.py", press enter
4. the code will execute the various simulations contained in the paper one by one, showing 2 or 3 plots of results each time. To proceed to the next simulation it is necessary to close the plots of the previous one (they will be automatically saved).


## Brief description of how the code function

The file "Model.py" is tha main file of the project, and contains all the variables that one might wish to modify when starting to use the code. In particular, among the first few lines you can find a list of tags that can be set either to "No" or to "Yes". This determines if the corresponding simulation will or not be run when executing the script Model.py.
We initially set all the tags to "Yes", except those related to the calibration procedure (MC, gradient search, RMS etc., wich can be extremely time consuming and we already runned them to determine once and for all the parameter set to be used for all the other simulations, which is contained in the file "OutputFileBestParametersSet.csv" in the folder containing the code) and those generating 3D plots (a little bit time consuming).


## Content of the file Model.py

Scrolling down the file Model.py, you will encounter the following content:

1. Tags ("Yes" or "No"), to chose which simulations to lunch.
2. Parameters sets and initial conditions. Notew that the values hard coded in the python script for the rate constants k are those called "fiducial set" in the paper, and are used as starting point for the calibration. Those that resulted from the calibration and are then used for all the simulations of point 4 below are contained in a file named "OutputFileBestParametersSet.csv" in the folder containing the code.
3. Functions necessary to the calibration of the model (MC and gradient search, described in the Supplementary Material of the paper)
4. Main simulations, including HSR, comparisons with experimental data, double HS, long term HS, trade off between duration and temperature of the HS, study of the steady state, hot day behaviour etc (they represent the simulations of the main body of the paper)

The fourth part is subdivided in several sections, one corresponding to roughly each simulation of the paper. the main parameters that one may wish to change for each simulation are contained in the corresponding section of part 4 of the file Model.py. In principle they can be changed and simulation relaunched without opening any of the other python files. These files contain the implementation of all the classes, functions etc. of the code and might need to be changed for more advanced developements of the code.

Let us note that the simulation studying the HSP accumulation as a function of duration and max temperature of the HS is initially set to have only 4x4=16 points, while the one of the paper has 30x30=900 (but then it is ways longer to produce).
