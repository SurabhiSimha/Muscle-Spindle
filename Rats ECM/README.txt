StarterCode README
Author: KPB
Date: 08/01/2017

This folder contains functions for curve fitting data collected in rat electrophysiology 
experiments. In total, there are 5 '.m' files to be used by MATLAB by default.

1) rat_operate.m
	- This script is the driver for running curve-fitting for different models and
	different trials. Currently, the user supplies which model the code should use for 
	curve fitting and which trial of the data file the actions should be performed on.
	This script will return a new data structure called fit_data, which contains relevant
	data from the optimization. 
	
2) findGains.m
	- This function finds the optimal gains that satisfy the objective function in 
	cost.m for the specified trial and model. For example, if the user specifies model = 1, 
	and trial = 5, findGains.m will take the raw data from the proc_data.mat file (loaded 
	beforehand by the user) from trial 5, and will attempt to fit gains and offsets using 
	fmincon.m to the data in the prescribed procedure (this depends on which model is 
	defined as "model 1". By default, model 1 is the force-related model, and the 
	procedure for estimating IFR in this case is defined in 'kinetics.m', though this can 
	be changed by the user). 
	
3) cost.m
	- This is the objective function defined to estimate the parameters of a model. By
	default, this function takes the current estimated firing rate by the model and the 
	recorded firing rate the model output is being fit to, and calculates the sum of 
	squared errors (SSE) divided by the sum of squares about the mean (SSM). fmincon.m 
	will attempt to find the minimum value of the cost, by changing model parameters until
	a minimum in the objective function is found (important to note this may not be a
	global minimum). 
	
4) kinetics.m
	- This is an example of a data-driven model that the curve-fitting code uses to
	predict firing data. This model estimates the components of force due to contractile
	and noncontractile mechanisms, and uses the former to predict muscle spindle firing.  
	
5) kinematics.m
	- This is another example of a data-driven model that the curve-fitting code uses to
	predict firing data. This model takes recorded length, velocity, and acceleration and
	creates a pseudolinear combination of them to predict muscle spindle firing. 