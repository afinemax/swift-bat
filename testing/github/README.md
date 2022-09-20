### swift-bat pipeline 

### By Maxwell A. Fine
### Contact: maxwell.fine@mail.utoronto.ca


To do:
* fix bug with generation of .RSP files
* add terminal interface
* change subproccess comands for mkdir to OS commands
* add proper spectrum modeling function that users can adjust
	* understand spectrum modeling, and why we are using which models 
	* pick "best" ones for default
* improve SNR search method
	* compare against making n lightcurves with the heafsoft software
* add additonal plots /  group together plots for a nice readable diagnostic plot
	* SNR vs time, (for each energy band), 
	* SNR vs fluence limit
	* SNR vs energy
	* SNR vs time bin size 
	* info about CHIME/FRB target and SWIFT/BAT data
* for SWIFT/BAT targets out of view / lacking the full 200s add a JSON entry saying so
* add SNR cutoff ?
* add back / add image search method
	* can only find sources 3.5 >= SNR on lowest setting
	* need to fix bug with input catalog not being read in correctly with heasoft
	* Do this as an inital search or on the best target after light curve?
		* needs event window, and background window
* make an output catalog function
	* outputs at the moment are stored individually in their output dirs	
	 	
 
Wish list to do:
* add additonal plots /  group together plots for a nice readable diagnostic plot
	* SNR vs time, (for each energy band), 
	* SNR vs fluence limit
	* SNR vs energy
	* SNR vs time bin size 
	* info about CHIME/FRB target and SWIFT/BAT data
* add SNR cutoff ?
* add back / add image search method
	* can only find sources 3.5 >= SNR on lowest setting
	* need to fix bug with input catalog not being read in correctly with heasoft
	* Do this as an inital search or on the best target after light curve?
		* needs event window, and background window
* make an output catalog function
* outputs at the moment are stored individually in their output dirs	
