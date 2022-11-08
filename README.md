# HIRES/KPF Automatic Scheduling/Script Generating Software
Generate scripts for HIRES/KPF observers for the next upcoming observing night by processing and optimizing a semester of data

#### Basic Requirements
* numpy
* astropy
* astroplan
* pandas

#### Computing Requirement
* Gurobi Python API (gurobipy)
  * Obtain a license for the ILP Solver https://www.gurobi.com/solutions/licensing/
    * There are both free trial licenses and free academic licenses
  * You can learn more about Gurobi here https://www.gurobi.com/documentation/quickstart.html

## Installation
Since the scheduer is not yet packaged, clone the git repository and open the file single_night_optimizer.py

Here you can adjust your data directories, but they are linked to appropriate example files included in the repo by default. To make your first
schedule, just run the file through command line
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
  * Directory of text files of format 'MM-DD-YYYY.txt' from all previous dates in observing semester
  * Copy the marked starlist targets relevant to the observers sheet (i.e. above the line of X's)
* *current_day*
 
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

### Example Usage

This repository comes equipped with data for the 2022B Semester up to October 4th. To generate a schedule, you would use the function:

```
semester_schedule(observers_sheet = 'Data/HIRES_observers - PI-requests2022B.csv',
                            twilight_times = 'Data/twilight_times.csv',
                            allocated_nights = 'Data/hires_schedule_2022B.csv',
                            marked_scripts = 'Data/MarkedScripts', 
                            current_day = '2022-10-04'
                            )
```

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
  
## Notes on the Formalism

Integer Linear Programming (ILP) is an abstraction of a logical problem using extraordinarily large matrices. The solver maximizes a linear objective
function subject to a series of constraints (inequality vectors). By keeping the objective simple, you can handle problems with millions of parameters!
 
### What does a parameter look like?

Parameters take many forms in the scheduler. The most efficient algorithms utilize binary variables to denote whether an action will take place or not.
Variables can take the form of matrices, such that we can index them to find the 1 or 0 corresponding to a certain decision.

For example, consider the variables assosciated with assigning a list of targets to a set of nights with no other complications. Let the matrix Y contain the
parameters that define what targets are assigned to what nights. The shape of the matrix is (num_targets,num_nights). If Y[target,night] == 1, the target is 
assigned to that night.

Coding this parameter in Gurobi might look like this:
```
import gurobipy as gp
from gurobipy import GRB

#The model object contains all of our information
m = gp.Model('Example')
Y = m.addVars(num_targets,num_nights, vtype = GRB.BINARY, name='YesOrNoMatrix')
``` 
See here: https://www.gurobi.com/documentation/9.5/refman/py_model_addvars.html

### How do constraints represent logic?
 
Once you create a framework of variables, you can constrain their values as linear inequalities, and enforce what form your solution can take. 

For example, lets say we want to enforce that each target can only be observed a single time across all our nights. The sum of the Y matrix across any row
should be less than or equal to one (scheduled once or less than once). 

We can constrain our Y matrix in Gurobi with:
```
for r in targets:
    m.addConstr((gp.quicksum(Y[r,t] for t in nights) <= 1), name='observe_once') 
```
More info here: https://www.gurobi.com/documentation/9.5/refman/py_model_addconstrs.html

### What is the objective?

In the scheduler, we maximize the sum of successfully scheduled observations, while reducing the objective for long periods between observations of the same 
target. We set a minimum duration, such as one specified in the **Cadence** column of the *observers_sheet*, and penalize for anything greater.
 
For our test problem here, let's do something simple. We'll tell Gurobi to maximize the number of scheduled observations like this:
```
#Maximize the number of 1's in our matrix
m.setObjective(gp.quicksum(Y[r,t] for r in reservations for t in nights),GRB.MAXIMIZE)
``` 
More here: https://www.gurobi.com/documentation/9.5/refman/cs_model_setobjective.html

Clearly this is not a realistic nor complete model, but Gurobi can solve it (for one of its countless solutions in this case). The solved model contains
parameters that now hold distinct values, and we can query these decision variables and convert them to a schedule.



Written by Luke Handley (lukehandley@g.ucla.edu)
