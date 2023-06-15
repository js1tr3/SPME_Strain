# README #

Lithium Ion Battery SPMe with Mechanical, Thermal and Degradation

### What is this repository for? ###

* This model runs SPMe added submodels for Temperature, Expansion and Degradation.
* Degradation mechanisms considered are: Lithium Plating, SEI formation and Mechanical Degradation
* This model is created using Matlab 2020b and requires Simulink 2020b to run.
* Version - 1

### How do I run the model? ###

* Set current and current mode in run_simulink_stress.m file
* Update model parameters in the parameters.m file (if necessary)
* Run run_simulink_stress.m 
* Run plot_plating.m to get Basic plots (needs to be updated for temperature, mechanical plots, concentration plots and SEI)
* Run plot_slider_pl.m to plot the side reaction overpotentials and current across the electrode as a function of time



### Who do I talk to? ###

* Sravan Pannala (spannala@umich.edu)
* Dr. Jason B. Siegel (siegeljb@umich.edu)