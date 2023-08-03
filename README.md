# HIRES/KPF Automatic Scheduling/Script Generating Software
Generate scripts for HIRES/KPF observers for the next upcoming observing night by processing and optimizing a semester of data

## Installation

Clone the git repository onto your local system with:
```
$ git clone https://github.com/lukehandley/KPFautoscheduler.git
```

#### Basic Requirements
* numpy
* matplotlib
* imageio
* astropy
* astroplan
* pandas
* gurobipy


### Install Gurobi and Obtain a License
  1. Download Gurobi for your OS here: https://www.gurobi.com/documentation/9.5/quickstart_windows/software_installation_guid.html
  2. Obtain an academic (non trial) license for the ILP Solver https://www.gurobi.com/solutions/licensing/
  3. You can learn more about Gurobi here https://www.gurobi.com/documentation/quickstart.html

You can install everything at once with pip using:
```
pip install -r requirements.txt
```
 
For Conda users, start by setting up a new environment in python 3.9:
```
conda create -n kpfauto python=3.9
```
Installing gurobipy will require first setting up the channel, then installing using:
```
conda config --add channels https://conda.anaconda.org/gurobi
```
```
conda install -c gurobi gurobi
```
Note that some external files may be downloaded while installing system requirements (i.e. ephemerides for the astroplan module)

### Using the scheduler
Custom gurobi tutorials are included in /tutorials. A reduced form of the HIRES semester data in 2023A (for easy tinkering) also lives here.

View the available commands with:
```
python scheduler.py -h
```

You can run the scheduler with data for the 2023B semester. These are set as the default parameters, so to generate a schedule for an allocated night on HIRES or KPF (simulating the given date and onwards), run:
```
python scheduler.py -i 'INSTRUMENT' -d 'YYYY-MM-DD'
```
Where INSTRUMENT is replaced by HIRES or KPF. The proper data directories will be automatically selected.

To initiate the 'high production mode' and generate fully formatted schedules for a sequence of nights, set the date parameter multiple times 
to append them all.
```
python scheduler.py -i 'INSTRUMENT' -d 'YYYY-MM-DD' -d 'YYYY-MM-DD'
```
And so on until every date is listed. Be sure to check there are no gaps between these dates (in the CPS observing schedule, they need not 
be consecutive dates in the year).


#### Data Inputs
* *observers_sheet* https://docs.google.com/spreadsheets/d/18r_xWaz26ya6sI0BQ6xpD3ibZLCwgkbjjAbwNmoxIUM/edit#gid=508439574
  * Requests/Cadence/Nobs/Exposure/Coordinates/Program included in optimization
  * Other columns needed for script formatting
* *twilight_times*
  * .csv with columns for twilight times for (at least) every day in the semester
  * Its much easier to query these than calculate them, the script 'twilight_calculator' is included to create one of these with astroplan
* *allocated_nights* https://github.com/California-Planet-Search/jump-config/tree/master/allocations/hires_j
  * Dates and start/stop markers for each quarter night
* *marked_scripts*
  * Directory of text files of format 'YYYY--MM-DD.txt' from all previous dates in observing semester
  * Includes the marked starlist targets relevant to the observers sheet (i.e. above the line of X's)
* *schedule_dates*
  * Date strings in same format as above, which will be slew optimized and script formatted
 
## The Scheduler

Organizes all observing requests from *observers_sheet* into a large set of possible observations and attempts to schedule as many as possible 
into the *allocated_nights* based on target accessibility. It uses Integer Linear Programming (https://en.wikipedia.org/wiki/Integer_programming)
to construct a mathematical formulation of the scheduling problem. This problem is then optimized to maximize observations at the desired cadence
using a third party solver (Gurobi) in real time.

This framework adapts throughout the semester based upon current progress. Completed observations are read from *marked_scripts*, and the
scheduler will construct an optimal semester plan moving forward. It generates three potential plans for the possibility of different weather conditions (nominal, weathered, poor) in the upcoming night for observers to swap between at their convenience. 

Targets are binned into discretized quarter nights, and those from the next upcoming night are run through an additional optimization. The scheduler solves a complex formulation of the traveling salesman problem to minimize the slews for the next night. The ordered list is then automatically formatted into scripts.

#### Outputs
* A csv is created for each weather condition for *current_day* and saved to the current directory for inspection. Each row represents a quarter night slot, and listed values are the assigned request numbers (indeces in *observers_sheet*) to that quarter night. These are purely for inspection/debugging.
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
  * Reads and understands current semester progress, then sources a new optimal schedule for the remaining dates




Written by Luke Handley (lukehandley@g.ucla.edu)