The EuLag-PATS is an Eulerian-Lagrangian model to solve transient state equations describing gravitational sinking of particulate multicomponent 
aggregates with organic (POC) and mineral (opals and CaCO3) components. The aggregate components degrade under different factors. 
In our model POC component is influenced by the age of aggregate, temperature and oxygen concentration. 
Opal component degradation is influenced by temperature only (the user can optionally include oxygen and age-dependence, if he finds it necessary). 
CaCO3 degradation rate is dependent on calcite saturation only or 0. 
The aggregate can be considered consisting 1, 2 or 3 different components. The user can also include other components, with same parameters.

EuLag-PATS.f90 is the source file. It was tested under Windows 11 using Intel Fortran Compiler with VS2008 shell.

Input files: 
Input1.nml - number of components, output path and filenames are defined here,
Input2.nml - model parameters are defined here. 
The model parameters should be specified in the order. 
We recommend 1 -POC component, 2 - opals, 3 - CaCO3.
We use this order for all parameters, but the user can define another order. 

*temperature profile file, ex. "Atlantic Temperature.dat"* - temperature profile (C), 
*oxygen concentration file, ex. "Atlantic Oxygen.dat"* - oxygen concentration profile (kg/m^3).
Optionally, the user can include a file for CO3 concentration and CO3 saturation concentration, if calcite degradation is included.

Output files:
The names of the output files are set in the format *name of variable*_*time in days*.dat. 
Files starting with Fd(z) are the total POM flux profile (kg/m^2/s), and Sp(z) is the total POM concentration profile (kg/m^3).
For example
Fd(z)_10.dat - total POM flux profile (kg/m^2/s) at t=10 days, 
Sp(z)_10.dat - total POM concentration profile (kg/m^3),  at t= 10 days, 
Fd(z)_norm_10.dat, Sp(z)_norm_10.dat - normalized total POM flux and concentration profiles at t=10 days. 

The filenames for concentrations and fluxes of separate components start as
Sp1(z)_..., Sp2(z)_.., Sp3(z)_.., 
Fd1(z)_..., Fd2(z)_.., Fd3(z)_.., 
and for normalized concentrations and fluxes of separate components
Sp1(z)_norm_.., Sp2(z)_norm_.., Sp3(z)_norm_..,
Fd1(z)_norm_.., Fd2(z)_norm_.., Fd3(z)_norm_..,
Here 1 - is a POC component, 2 - opals, 3 - CaCO3. 
We use this order for all parameters, but a user can define another order. 
For example Sp1(z)_norm_20.dat is normalized POC concentration at t=20 days;
Fd3(z)_norm_365.dat is normalized CaCO3 flux at t=365 days;


How to use:
Step 1: input the parameters into input-files. First, input the number of fractions and filenames for temperature and oxygen concentration in the Input1.nml file. 
Also set filenames for CO3 concentration and CO3 saturated concentration if calcite degradation is considered.
Follow the instructions to input the data correctly. In this release, parameters, temperature and oxygen profiles are set to those that were used in our simulations for the Atlantic Ocean.


Step 2: input the model parameters into Input2.nml. Follow the instructions to input the data correctly. 
The aggregate can be considered to consist of 1, 2, or 3 different components. 
The user can also include other components with the same parameters.
The model parameters should be specified in order. 
We recommend 1 -POC component, 2 - opals, 3 - CaCO3.
We use this order for all parameters, but the user can define another order. 

Note that the program has different working modes.
a) The mode parameter is a velocity calculation mode. 

If mode=0, then velocity is calculated with the general nonlinear Stokes formula with an iterative calculation algorithm;

If mode=1, then the simplified power-law formula for Re<<1 is used, with parameters Cw and η calculated in the program;

If mode=2, then the simplified power-law formula for Re<<1 is used, with Cw and η predefined by the user.

b) The modecaco3 parameter is a calcite degradation parameter;

If modecaco3=1, then calcite degradation dependent on calcite saturation is included

If modecaco3=0, then calcite degradation is not included.

c) The mode_eps is a parameter to define the epsilon function. 
The basic epsilon value is epsilon = const, defined by the epsilon parameter. 
If the user wants to define epsilon=const, he can set mode_eps=0 or any other value, but not 1 or 2.

If  mode_eps=1, epsilon=eps_min if time <= time1 parameter, epsilon increase
linearly from eps_min to eps_max if time1 <= time <= time2, and epsilon=eps_max if time >= time2.

If  mode_eps=2, epsilon=eps_max if time <= time1 and if time>=time1+delta_t2; 
epsilon decreases linearly from eps_max to eps_min during the period of time 
time1 <= time <= time1+delta_t1, and increases linearly from eps_min to eps_max 
during the period of time time1+delta_t1<= time <= time1+delta_t2. 
This mode is used to simulate seasonal spectrum variations caused by algae bloom. 

d) The mode_init parameter is used to define the starting concentration and flux profile. 
In numerical simulations with a 1-component model, an analytical solution can be used as a starting profile.

If mode_init=0, or any other value then simulations start with 0 profile (default setting).

If mode_init=1, and number of components nf=1, an analytical solution is used as starting profile.

e) The datamode parameter specifies the way the boundary conditions are defined;

If datamode=0, then the boundary conditions are defined by total POM concentration (Sp_measur);

If datamode=1, then the boundary conditions are defined by total POM flux (Fd_measur); 

If datamode=2, then the boundary conditions are defined by separate component concentrations (conc array); 

If datamode=3, then the boundary conditions are defined by separate component fluxes (flux array); 
For datamode = 2 and 3, the component mass fractions pM are also defined by separate component concentrations or fluxes. 

Step 3: If you want to model some specific ocean and change temperature and oxygen concentration profiles, replace the corresponding temperature and oxygen data files and corresponding measurements data (Z_measur and Sp_measur, or Fd_measur, or flux, or conc..) in the input files. You can found measurements data in the https://github.com/KKovalets/EuLag_DataSet repository or in the EuLag_DataSet release https://github.com/KKovalets/EuLag_DataSet/releases/tag/v0.0.0. 
Change the corresponding parameters n_data1 and n_data2 in the input file. n_data1 is the number of lines in the temperature and n_data2 is the number of lines in the oxygen concentration file. n_data3 is the number of lines in calcite degradation parameter files, such as CO3 concentration and CO3 saturation concentration, if calcite degradation is included.

Step 4: Compile an executable file using EuLag-PATS.f90 and run it.
