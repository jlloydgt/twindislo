## Synopsis

Twindislo is a python code for calculating the growth of twins within and across grain boundaries using dislocation statics. The code accompanies the paper "A dislocation-based model for twin growth within and across grains", submitted to Proceedings of the Royal Society: A, by Jeffrey T Lloyd. 

## Code Examples

The easiest, fastest example is running "python test.py", which is the same file was paper-scripts/1um-tau10.py. This will generate a twin in a 1um grain subjected to a shear stress of 10 MPa, using material properties for magnesium and quarter symmetry in the folder 1um/tau10/. The anticipated output to terminal is: 

Step 1 used 836 iterations__
Step 2 used 1599 iterations__
Step 3 used 1513 iterations__
Step 4 used 1433 iterations__
Step 5 used 1413 iterations__
Step 6 used 245 iterations

One more example is the same set of simulations, but with quarter symmetry turned off, which can be evaluated by running the file "python test-nosym.py". The anticipated output to the terminal is:

Step 1 part 1 used 339 iterations__
Step 1 part 2 used 132 iterations__
Step 2 part 1 used 474 iterations__
Step 2 part 2 used 330 iterations__
Step 3 part 1 used 190 iterations__
Step 3 part 2 used 279 iterations__
Step 4 part 1 used 494 iterations__
Step 4 part 2 used 335 iterations__
Step 5 part 1 used 803 iterations__
Step 5 part 2 used 422 iterations__
Step 6 part 1 used 44 iterations__
Step 6 part 2 used 1279 iterations

For some of the simulations, the code does not converge properly once annihilation occurs. At this point the user may terminate the program manually or let it run and visualize when this occurs. 

## Motivation

The code is used to calculate the growth of twins within and across grain boundaries using a dislocation statics approach. Right now the code is only used for a single twin with a shear stress applied directly to the twin direction, but it would be very easy to extend to multiple twins. The twin stops growing when the resolved shear stress at the nucleation location becomes negative.  

## Installation

The code was only tested on python 2.7.9 with matplotlib version 1.5.1 and numpy version 1.11.3. The code saves images as .png, but if you want to use .jpeg, you will need PIL installed. The code was testing on linux and mac - the file saving will need to be altered for windows. 

## Tests

All of the files needed to generate the data for figures in the manuscript are located in paper-scripts/ under the appropriate folders. These need to be placed in the main directory to be run using python filename.py, just as was done for the test problems. 

## Contributors

Jeff Lloyd - jeffrey.t.lloyd.civ@mail.mil. The step algorithm is better explained by Richard LeSar's notes "Mesoscale Dynamics: Dislocation Dynamics examples", which you can find on his website. If you would like to contribute or collaborate on this project, please contact me and I'd be glad to give some potential directions.

## License
    The program twindislo is licensed under the Creative Commons Zero 1.0 Universal (CC0 1.0) Public Domain Dedication license.

