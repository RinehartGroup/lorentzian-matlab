# lorentzian-matlab
A MATLAB script for Lorentzian fitting of the derivative of magnetic hysteresis data from Quantum Design instrumentation.

## Setup
The script should be added to your MATLAB path. This can be done by clicking the **Set Path** button under the **Home** tab. See [this](https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html) for more details.

## Usage
This script uses MPMS3 .dat files and will differentiate between measurement types. Typically, DC mode is used for nanoparticle samples, while VSM mode is used for molecular samples.

1. The user is prompted for the sample identity. Do not include the ".dat" portion of your datafile when typing it into the command window.
2. The script will display the available temperatures and prompt the user to select one. 
3. The user is prompted to fit either the forward or reverse scan
4. The user is then prompted to input how many peaks to fit. This will take some manual inspection.

The script will output three plots and save each as a .png:
  Figure 1: Full hysteresis loop
  Figure 2: Full hysteresis loop and the derivative of the forward or reverse scan with corresponding Lorentzian fit
  Figure 3: Full hysteresis loop and the corresponding Lorentzian fit, accounting for multiple peaks (if included)
  
The script will save each plot as a .fig and save the fit parameters as .mat for later use.

## Example
**_Data import:_**
```
>> Lorentzian
Sample: IronOxide_14nm_MvsH

```
**_Data selection:_**
```
Available Temps: 
     5
   300
   
Which temperature do you want to fit? 5

Forward (1) or reverse (2) scan? (input 1 or 2): 1

How many peaks do you want to fit? 1
```

