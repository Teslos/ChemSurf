This documents gives some description of the
ChemSurf.for program.
This program is copied from paper Tanaka et al. "Calculation of Surface Tension of Liquid Bi-Sn Alloy using 
Thermochemical Application Library ChemApp", Calphad, Vol24. No.4, pp 485-474,2000.
It works with ThermoCalc ver. S TQ-interface. Be aware TQ interface is buggy piece of software!
So you allways need to check your results. There are many hacks in this program to keep it running with TQ.
I suppose that with new version a new hacks will be needed and old abolished :)

So beware you have been warn!

What is needed to run the program:
- Put your myGES.GES5 data file in the same directory where the program is.
- Create input file which look something like this:

The parameters are put in the 
# Name of the phase diagram in this
# case it is name of thermocalc GES5 workspace.
CECO.GES5
# Name of the liquid phase
IONIC_LIQ
# Beta-mix is parameter which controls the correction
# of the excess Gibbs energy of the solution.
0.83
# Name of the solid phase 1
CO
# Coefficients of the surface energy equation sigmma = a+b*T and melting temperature Tm
2.82  0.00007  1200.0
# Molar volume of liquid and solid phase1
0.0000208   0.000117
# Name of the solid phase 2
CE
# Coefficients of the surface energy equation sigmma = a+b*T and melting temperature Tm
1.00  0.00009  1200.0
# Molar volume of liquid and solid phase2
0.000017  0.00087
# Max temperature
2500
# step for concentration
0.05

- run the program using run.bat batch script

