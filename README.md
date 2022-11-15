# HIRES/KPF Automatic Scheduling/Script Generating Software
Generate scripts for HIRES/KPF observers for the next upcoming observing night by processing and optimizing a semester of data

## Installation

Clone the git repository onto your local system with
```
$ git clone https://github.com/lukehandley/KPFautoscheduler.git
```

#### Basic Requirements
* numpy
* astropy
* astroplan
* pandas
* gurobipy



### Install Gurobi and Obtain a License
  1. Download Gurobi for your OS here: https://www.gurobi.com/documentation/9.5/quickstart_windows/software_installation_guid.html
  2. Obtain an academic (non trial) license for the ILP Solver https://www.gurobi.com/solutions/licensing/
  3. You can learn more about Gurobi here https://www.gurobi.com/documentation/quickstart.html

For pip, you can use 
```
python -m pip install gurobipy
```
 
On Conda, first set up the channel, then install gurobipy using:
```
conda config --add channels https://conda.anaconda.org/gurobi
```
```
conda install -c gurobi gurobi
```

### Using the scheduler
Gurobi tutorials are included in /tutorials.

View the available commands with:
```
python single_night_optimizer.py -h
```

You can run the scheduler with included example data for the 2022B semester valid through 10/04/22. These are set as the default parameters, so just run:
```
python single_night_optimizer.py
```

#### Data Inputs
* *observers_sheet* https://docs.google.com/spreadsheets/d/18r_xWaz26ya6sI0BQ6xpD3ibZLCwgkbjjAbwNmoxIUM/edit#gid=508439574
  * Requests/Cadence/Nobs/Exposure/Coordinates/Program included in optimization
  * Other columns needed for script formatting
* *twilight_times*
  * .csv with columns for twilight times for at least every day in the semester
  * Its much easier to query these than calculate them, a function to create one of these with astroplan is included
* *allocated_nights* https://github.com/California-Planet-Search/jump-config/tree/master/allocations/hires_j
  * Dates and start/stop markers for each quarter night
* *marked_scripts*
  * Directory of text files of format 'YYYY--MM-DD.txt' from all previous dates in observing semester
  * Copy the marked starlist targets relevant to the observers sheet (i.e. above the line of X's)
* *current_day*
  * Date string in same format as above
 
## The Scheduler

Organizes all observing requests from *observers_sheet* into a large set of possible 'reservations' and attempts to schedule as many as possible 
into the *allocated_nights* based on target accessibility. It uses Integer Linear Programming (https://en.wikipedia.org/wiki/Integer_programming)
to construct a mathematical formulation of the scheduling problem. This problem is then optimized to maximize observations at the desired cadence
using a third party solver (Gurobi) in real time.

This framework adapts throughout the semester based upon current progress. Completed observations are read from *marked_scripts*, and the
scheduler will construct an optimal semester plan moving forward. It generates three potential plans for the possibility of different weather conditions
(nominal, weathered, poor) in the upcoming night for observers to swap between at their convenience. 

Targets are binned into discretized quarter nights, and those from the next upcoming night are run through an additional optimization. The scheduler
solves the traveling salesman problem (see https://www.gurobi.com/resource/traveling-salesman-problem/) to minimize the slews for the next night. The 
ordered list is then automatically formatted into scripts.

#### Outputs
* A csv is created for each weather condition for *current_day* and saved to the current directory for inspection. Each row represents a quarter night
slot, and listed values are the assigned request numbers (indeces in *observers_sheet*) to that quarter night. These are purely for inspection/debugging.
* A text file for *current_day* from each of the semester schedules formatted as a starlist using the column values from *observers_sheet*.

Note that *current_day* is the allocated night you are planning for and must exist in *allocated_nights*

The process for updating data for each new day:
* Download the newest observers sheet (see inputs)
* Copy the marked script from the previous night (or nights) and save them to text files with appropriate dates
* Update the day you are scheduling for in the generator function

### Features

* There are additional columns that can be edited in *observers_sheet* to customize the output of the generator.

  * **Comment**: Text is imported directly into the scriptlines. This is better used for observing tips than making notes about the semester 
  progress of a target. For example, 'companion at 0.7", observe brighter component (SW)'.
  * **Include**: This column allows you to forcibly constrain what targets should be observed in the upcoming night. An 'N' in this column will
  block the target for the upcoming night. Equivalently a 'Y' would force a target to be observed, however this is broken in the current version.
  * **Done**: Similar to **Include**, but for the *entire* semester. If column is marked 'Done' target will no longer be considered.
  * **Cadence**: Specify the amount of days between observations as an integer. If one is not specified, the scheduler will construct one that 
  provides the most uniform sampling throughout the semester, calculated using target accessibility and the 'Nobs' column in *observers_sheet*.

* Weather
  * Takes potential weather losses into account by randomly sampling out 30% of quarter nights in the remaining semester
  * Produces three scripts for the upcoming night based on conditions. Accounts for increased exposure time and in the future will enforce a
  maximum magnitude limit
  
* Program Equity
  * The scheduler can also enforce uniformity of program completion using the optional argument *equalize_programs* when set to True. This is
  generally discouraged as it can reduce total scheduled observations

* Script Processing
  * Reads and understands current semester progress, then reacts with optimality




Written by Luke Handley (lukehandley@g.ucla.edu)
