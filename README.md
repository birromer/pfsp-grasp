# Permutational Flowshop Scheduling Problem - GRASP

In order to run the meta-heuristic you'll need the Julia programming language installed,
as well the libraries Random and Statistics

You can install those libraries with the following commands:
* Open the julia console by typing `julia` in a terminal
* Press `]` and write the command `add Random Statistics`
* After that you can close the REPL by typing a backspace followed by `exit()`

To run the GRASP implementation you'll need to supply the following parameters:
* filename: the txt file with the t matrix.
* ALPHA: the percentage of the candidate set to be considered for selecting a
  random initial solution.
* NUMBER_CANDIDATES: the amount of random solutions that will be generated in
  the constructor.
* STOP_GRASP: how many iterations of GRASP to run.
* STOP_HILL_CLIMBING: how many iterations without imprevement in hill climbing
  are needed for it to stop.
* SEED: a seed for the random functions calls

The parameters used in our assigment are:
* ALPHA = 0.1
* NUMBER_CANDIDATES = 270
* STOP_GRASP = 1500
* STOP_HILL_CLIMBING = 50
* SEED = 10 (and increased by 10 every new run)

And the execution of that with the basic instance would be
`julia grasp.jl VFR10_15_1_Gap.txt 0.1 270 1500 50 10`
