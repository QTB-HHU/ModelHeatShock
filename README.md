# ModelHeatShock

To run our implementation of our Heat shock model:
1) Download the full content of the directory containing this file.
2) open a terminal and type "python3 Model.py"
3) the code will execute the various simulations contained in the paper one by one, showing 2 or 3 plots of result each time. To proceed to the next simulation it is necessary to close the plots of the previous one (they will be automatically saved).

The file "Model.py" is tha main file of the project, and contain all the variables that one might wish to modify when starting to use the code.
In particular, among the first few lines you can find a list of tags that can be set either to "No" or to "Yes". This determines if the corresponding simulation will or not be run when executing the script Model.py.
We initially set all the tags to "Yes", except those related to the calibration procedure (MC, gradient search, RMS etc., wich can extremely time consuming and we already runned them to determine once and for all the parameter set to be used for all the other simulations) and those generating 3D plots (a little bit time consuming).
In addition, the simulation studying the HSP accumulation as a function of duration and max temperature of the HS is initially set to have only 4x4=16 points, while the one of the paper has 30x30=900 (but then it is ways longer to produce).
