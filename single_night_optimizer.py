from collections import defaultdict
from gurobipy import *
import numpy as np
import astropy.units as u
import pandas as pd
import gurobipy as gp
import math
import astropy as apy
import astroplan as apl
from astropy.time import Time
from astropy.time import TimeDelta
from itertools import combinations
import csv

def generate_twilight_sheet():

    #Initialize the astroplan observer object
    keck = apl.Observer.at_site('W. M. Keck Observatory')

    #Create a range of dates to calculate (some before and some after semester in this case)
    twilight_frame = pd.date_range('2021-06-01','2023-06-01').to_frame()
    twilight_frame = twilight_frame.rename(columns={0:'time_utc'})
    eighteen_deg_evening = []
    twelve_deg_evening = []
    six_deg_evening = []
    eighteen_deg_morning = []
    twelve_deg_morning = []
    six_deg_morning = []
    for day in twilight_frame.index.strftime(date_format='%Y-%m-%d').tolist():
        as_day = Time(day,format='iso',scale='utc')
        eighteen_deg_evening.append(keck.twilight_evening_astronomical(as_day,which='next'))
        twelve_deg_evening.append(keck.twilight_evening_nautical(as_day,which='next'))
        six_deg_evening.append(keck.twilight_evening_civil(as_day,which='next'))
        eighteen_deg_morning.append(keck.twilight_morning_astronomical(as_day,which='next'))
        twelve_deg_morning.append(keck.twilight_morning_nautical(as_day,which='next'))
        six_deg_morning.append(keck.twilight_morning_civil(as_day,which='next'))

    twilight_frame['18_evening'] = eighteen_deg_evening
    twilight_frame['12_evening'] = twelve_deg_evening
    twilight_frame['6_evening'] = six_deg_evening
    twilight_frame['18_morning'] = eighteen_deg_morning
    twilight_frame['12_morning'] = twelve_deg_morning
    twilight_frame['6_morning'] = six_deg_morning

    #Create simple csv with date and twilight evening/morning times
    twilight_frame.to_csv('twilight_times.csv')

def generate_reservation_list(all_targets_frame,plan,twilight_frame,observatory):
    dates = plan.Date.to_list()

    #Reservations is a list of tuples of the form (r,t), where r is the request id, and t is a qn slot
    #in which the request target is accessible
    reservations = []

    #The minimium separation list will contain either the PI specified cadence, or one generated if necessary
    min_separations = []

    #To quickly calculate observability, we use the start date and calculate minute by minute to create a 1D binary
    #array of accessibility for one day
    start = Time(dates[0],format='iso',scale='utc')
    stop = start + TimeDelta(1,format='jd')
    step = TimeDelta(60,format='sec')
    t = np.arange(start.jd,stop.jd,step.jd)
    t = Time(t,format='jd')

    #These coordinates designate observability at Keck
    min_az = 5.3
    max_az = 146.2
    min_alt = 33.3
    else_min_alt = 18
    max_alt = 85

    #Go request by request, and find which quarter nights we can throw each target into
    for index,row in all_targets_frame.iterrows():
        name = row['Starname']
        ra = row['ra']
        dec = row['dec']
        dur = row['discretized_duration']

        #Standard astroplan procedure
        coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
        target = apl.FixedTarget(name=name, coord=coords)
        AZ = observatory.altaz(t, target)
        alt=AZ.alt.deg
        az=AZ.az.deg

        #Filter by observability constraints, assign '1' to good minutes out of the day, otherwise '0'
        deck = np.where((az >= min_az) & (az <= max_az))
        deck_height = np.where((alt <= max_alt) & (alt >= min_alt))
        first = np.intersect1d(deck,deck_height)
        not_deck_1 = np.where((az < min_az))
        not_deck_2 = np.where((az > max_az))
        not_deck_height = np.where((alt <= max_alt) & (alt >= else_min_alt))
        second = np.intersect1d(not_deck_1,not_deck_height)
        third = np.intersect1d(not_deck_2,not_deck_height)
        good = np.concatenate((first,second,third))
        a = np.zeros(len(t),dtype=int)
        a[good] = 1

        #Now that we have a starting observability array, we can quickly do the same for any day we want by
        #permuting this array by the change in set times due to sidereal time
        for ind,r in plan.iterrows():
            date = dates[ind]
            day_diff = (Time(date)-start).value

            #This shift is easy to compute and keeps us accurate within a minute at any time
            offset = math.floor(day_diff/15)
            shift = day_diff * 4 - offset

            #Find the starting/ending minute for each interval in our plan to use as accessibility bounds
            b = int(np.round((Time(r['qn_start'],format='jd')-
                Time(date,format='iso',scale='utc')).jd * 24*60,0) + shift)
            c = int(np.round((Time(r['qn_stop'],format='jd')-
                Time(date,format='iso',scale='utc')).jd * 24*60,0) + 1 + shift)

            #Permute the accessibility array
            if c > len(a):
                one = a[b:]
                remainder = c % len(a)
                two = a[:remainder]
                d = np.concatenate((one,two))
            else:
                d = a[b:c]

            #Here I've decided to make it so that a target must be accessible for at least half of a quarter night
            #to make the reservation schedulable. This is crucial for time independent traveling salesman later
            if (len(d) - np.bincount(d)[0]) >= len(d)/2:
                reservations.append((index,ind))

        #This generates a minimum separation for targets that don't have one specified by the PI
        if math.isnan(row['Cadence']):

            #Only applicable for cadenced observations
            if row['N_obs(full_semester)'] > 1:

                #Find the number of total days in the semester that the target is up
                days_observable = 0
                for day in pd.date_range(dates[0],dates[-1]).strftime(date_format='%Y-%m-%d').tolist():
                    
                    #Similar process as above
                    day_start = Time(day,format='iso',scale='utc')
                    day_diff = (day_start-start).value
                    offset = math.floor(day_diff/15)
                    shift = day_diff * 4 - offset
                    morning_twilight = Time(twilight_frame.loc[day,'12_morning'],format='jd')
                    evening_twilight = Time(twilight_frame.loc[day,'12_evening'],format='jd')
                    b = int(np.round((evening_twilight - day_start).jd * 24*60,0) + shift)
                    c = int(np.round((morning_twilight - day_start).jd * 24*60,0) + 1 + shift)
                    if c > len(a):
                        one = a[b:]
                        remainder = c % len(a)
                        two = a[:remainder]
                        d = np.concatenate((one,two))
                    else:
                        d = a[b:c]
                    #For this calculation the observability constraint is relaxed. We just look to see if it could be
                    #observed at any point at all, not for a major chunk
                    if (len(d) - np.bincount(d)[0]) >= dur:
                        days_observable += 1

                #This is the formula to calculate separation based on accessible days and Nobs
                if int(np.round(0.8 * days_observable/row['N_obs(full_semester)'],0)) - 3 > 1:
                    min_separations.append(int(np.round(0.8 * days_observable/row['N_obs(full_semester)'],0))-3)
                    row['Cadence'] = min_separations[-1]
                else:
                    min_separations.append(1)
                    row['Cadence'] = min_separations[-1]
            
            #For single shots
            else:
                min_separations.append(0)
                row['Cadence'] = min_separations[-1]

        #For PI specified cadences    
        else:
            min_separations.append(row['Cadence'])
    return reservations,min_separations

def calculate_intervals(plan,twilight_frame):
    dates = plan.Date.to_list()

    #starts and stops are the fractional beginning and ending to each quarter night in the plan
    starts = plan.start.to_numpy()
    stops = plan.stop.to_numpy()
    night_starts = []
    night_ends = []

    #Move through each date in the plan, calculate the start and end times for each qn using
    for date in dates:
        day_start = Time(date,format='iso',scale='utc')

        #Queary the twilight times to find the start and end to each total night
        morning_twilight = Time(twilight_frame.loc[date,'12_morning'],format='jd')
        evening_twilight = Time(twilight_frame.loc[date,'12_evening'],format='jd')

        #Round to nearest minute
        length1 = np.round((evening_twilight-day_start).jd*24*60,0)*60
        length2 = np.round((morning_twilight-day_start).jd*24*60,0)*60
        night_starts.append((day_start+TimeDelta(length1,format='sec')).jd)
        night_ends.append((day_start+TimeDelta(length2,format='sec')).jd)
    
    #Find the duration of each night
    durations = []
    for i in range(len(dates)):

        #diff is the minute length of each total night assigned to each date in the plan
        diff = night_ends[i]-night_starts[i]
        observing_starts = night_starts[i] + diff*starts[i]
        observing_ends = night_starts[i] + diff*stops[i]
        durations.append(math.floor(((observing_ends-observing_starts)*24)*60))

    #Finally assign the rouned beginning and ending of each quarter night in the plan
    interval_starts = []
    interval_stops = []
    for i in range(len(dates)):
        diff = night_ends[i]-night_starts[i]
        start = night_starts[i] + diff*starts[i]
        end = (Time(start,format='jd') + TimeDelta(60*(durations[i]),format='sec')).jd
        interval_starts.append(start)
        interval_stops.append(end) 
    for i in range(1,len(interval_stops)):
        start_min = int(np.round((Time(interval_starts[i],format='jd')-
            Time(plan.Date.values[i],format='iso',scale='utc')).jd * 24*60,0))
        end_min = int(np.round((Time(interval_stops[i-1],format='jd')-
            Time(plan.Date.values[i-1],format='iso',scale='utc')).jd * 24*60,0))

        #Make sure rounding doesn't lead to any overlap between consecutive quarter nights
        if start_min == end_min:
            interval_stops[i-1] = ((Time(interval_stops[i-1],format='jd') - TimeDelta(60,format='sec')).jd)

    plan['qn_start'] = interval_starts
    plan['qn_stop'] = interval_stops
    return plan

def discretize(seconds):
    return np.round(seconds/60,1) + 1

def can_force(forced_slots,reserved_slots):
    for x in forced_slots:
        for y in reserved_slots:
            if x == y:
                return True
    return False

def relaxed_minimum(minimum,relaxation):
    return math.ceil(minimum * relaxation)

def process_scripts(all_targets_frame,plan,marked_scripts):
    import os
    observed_dict = defaultdict(list)

    #The Marked Scipts directory includes the lines from the previous nights in the semester organized into files
    #titled by date. i.e 'marked_scripts/2022-08-07.txt'
    #Note: these files should include only targets for that night, and none of the text/notes that follow
    directory = marked_scripts
    for filename in os.listdir(directory)[1:]:

        fpath = os.path.join(directory, filename)
        with open(fpath) as f:

            #Read in each processed script line
            lines = f.readlines()

            #Parse the corresponding date
            day = f.name[:-4][-10:]

            #The observation is always assigned to the first qn of that date for simplicity. The date of observation
            #all that matters for the optimization, not the specific qn
            slot = plan[plan['Date'] == day].index[0]
            good_lines = []

            #lines that are marked are completed observations
            for line in lines:
                if (line[0] == 'x') or (line[0] == 'X'):
                    good_lines.append(line)

            #Search for corresponding request by cross referencing name and program
            #This is a uniqueness problem that still needs to be solved in Jump
            for line in good_lines:
                name = 'default'
                program = 'default'
                for item in line.split(' '):
                    
                    #Sift through the text to find valid name/program identifiers
                    if item in all_targets_frame['Starname'].tolist():
                        name = item
                    if item in all_targets_frame['Program code'].tolist():
                        program = item
                
                #See if this target appears in our sheet
                targ = all_targets_frame.loc[(all_targets_frame['Starname'] == name) 
                                & (all_targets_frame['Program code'] == program)]['request_number'].values
                
                #If we can find it, append the corresponding qn identifier to the targets
                #past observation dictionary
                if targ.size > 0:
                    observed_dict[targ[0]].append(slot)

    return observed_dict

def conflicting_constraints(model,time_limit):

    #Utilize Gurobi's built in tools for finding model infeasibility
    model.Params.TimeLimit = time_limit
    search = model.computeIIS()

    #Print the names of all conflicting constraints for debugging
    for c in model.getConstrs():
        if c.IISConstr:
            print('%s' % c.ConstrName)
    for c in model.getGenConstrs():
        if c.IISGenConstr:
            print('%s' % c.GenConstrName)

def coordinates(frame,plan,targ,qn):

    #Target coordinates are estimated at halfway through the quarter night
    time = Time((plan.loc[qn,'qn_start']+plan.loc[qn,'qn_stop'])/2,format='jd')
    ra = frame.loc[targ,'ra']
    dec = frame.loc[targ,'dec']
    observatory = apl.Observer.at_site('W. M. Keck Observatory')
    name = frame.loc[targ,'Starname']
    coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
    target = apl.FixedTarget(name=name, coord=coords)
    AZ = observatory.altaz(time, target)
    
    return AZ

def distance(frame, plan, targ1, targ2, qn):
    t1 = coordinates(frame,plan,targ1,qn)
    t2 = coordinates(frame,plan,targ2,qn)

    #If you wanted to calculate pairwise distances including elevation
    #diff = (t1.az.deg-t2.az.deg, t1.alt.deg-t2.alt.deg)
    #return math.sqrt(diff[0]*diff[0]+diff[1]*diff[1])

    return np.abs(t1.az.deg-t2.az.deg)

def create_coordinates(name,ra,dec,time,observatory):
    
    #I've included hour angle here for future implementation, however it's not currently used
    coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
    target = apl.FixedTarget(name=name, coord=coords)
    AZ = observatory.altaz(time, target)
    HA = observatory.target_hour_angle(time, target)
    
    return AZ,HA

def write_starlist(frame,requests,condition,current_day):

    #Retrieve columns from observers sheet that are relevant to the starlist
    columns = ['Starname','RAH','RAM','RAS','DECD','DECM','DECS','epoch','vmag=','Vmag',
        'T_exp(sec)','T_max(sec)','Exp_meter_counts(k)','decker','N_obs','Iodine',
        'Priority','Program code','Telescope Propsoal Code','Comment']

    formatted_frame = frame.loc[requests][columns]

    lines = []
    for index,row in formatted_frame.iterrows():

        #Just a bunch of string formatting. This prints standard starlists as ordered by the salesman optimization
        namestring = ' '*(16-len(row['Starname'])) + row['Starname']
        
        rastring = ('0'*(2-len(str(int(row['RAH'])))) + str(int(row['RAH'])) + ' '
                        + '0'*(2-len(str(int(row['RAM'])))) + str(int(row['RAM'])) + ' '
                            + '0'*(4-len(str(row['RAS']))) + str(row['RAS']))
        
        starter = '+'
        if int(row['DECD']) < 0:
            starter = '-'
        decstring = (starter + '0'*(2-len(str(abs(int(row['DECD']))))) + str(abs(int(row['DECD']))) + ' '
                        + '0'*(2-len(str(int(row['DECM'])))) + str(int(row['DECM'])) + ' '
                            + '0'*(4-len(str(row['DECS']))) + str(int(row['DECS'])))
        
        magstring = row['vmag='] + str(row['Vmag']) + ' '*(5-len(str(row['Vmag'])))
        
        exposurestring = (' '*(4-len(str(int(row['T_exp(sec)'])))) + str(int(row['T_exp(sec)'])) + '/' 
                            + str(int(row['T_max(sec)'])) + ' '*(4-len(str(int(row['T_max(sec)'])))))
        
        countstring = (' '*(3-len(str(int(row['Exp_meter_counts(k)']))))
                            + str(int(row['Exp_meter_counts(k)'])) + 'k')
        
        line = (namestring + ' ' + rastring + ' ' + decstring + ' ' + str(int(row['epoch'])) + ' '
                    + magstring + ' ' + exposurestring + ' ' + countstring + ' ' + row['decker'] + ' ' + 
                    str(int(row['N_obs'])) + 'x' + ' ' + ' '*(3-len(row['Iodine'])) + row['Iodine']
                    + ' ' + row['Priority'] + ' ' + 
                    row['Program code'] + ' ' + row['Telescope Propsoal Code'])
        
        if not pd.isnull(row['Comment']):
            line += (' ' + str(row['Comment']))
        
        lines.append(line)

        #Formatted starlists appear as text files in the directory
        with open('{}_{}.txt'.format(current_day,condition), 'w') as f:
            f.write('\n'.join(lines))

def salesman_scheduler(all_targets_frame,plan,current_day):
    for condition in ['nominal','weathered','poor']:

        #Retrieve the conditions from respective csv's
        file = open("2022B_{}.csv".format(condition), "r")
        blocked_targets = list(csv.reader(file, delimiter=","))
        file.close()
        for i in range(len(blocked_targets)):
            for j in range(len(blocked_targets[i])):
                blocked_targets[i][j] = int(blocked_targets[i][j])


        #Functions unique to Gurobi's salesman implementation. I've only touched these enough to make them work
        def subtourelim(model, where):
            if where == GRB.Callback.MIPSOL:
                # make a list of edges selected in the solution
                vals = model.cbGetSolution(model._vars)
                selected = gp.tuplelist((i, j) for i, j in model._vars.keys()
                                    if vals[i, j] > 0.5)
                # find the shortest cycle in the selected edge list
                tour = subtour(selected)
                if len(tour) < len(blocked_targets[qn]):
                    # add subtour elimination constr. for every pair of cities in subtour
                    model.cbLazy(gp.quicksum(model._vars[i, j] for i, j in combinations(tour, 2))
                                <= len(tour)-1)

        def subtour(edges):
            unvisited = blocked_targets[qn][:]
            cycle = blocked_targets[qn][:] # Dummy - guaranteed to be replaced
            while unvisited:  # true if list is non-empty
                thiscycle = []
                neighbors = unvisited
                while neighbors:
                    current = neighbors[0]
                    thiscycle.append(current)
                    unvisited.remove(current)
                    neighbors = [j for i, j in edges.select(current, '*')
                                if j in unvisited]
                if len(thiscycle) <= len(cycle):
                    cycle = thiscycle # New shortest subtour
            return cycle    


        ordered_requests = []

        #The path for each quarter night is computed independently
        for qn in plan[plan['Date'] == current_day].index.tolist(): 

            #Traveling Salesman https://www.gurobi.com/resource/traveling-salesman-problem/

            dist = {(t1, t2): distance(all_targets_frame,plan,t1, t2, qn) for t1, t2 in combinations(blocked_targets[qn], 2)}
            # tested with Python 3.7 & Gurobi 9.0.0

            m = gp.Model()
            m.Params.OutputFlag = 0
            # Variables: is city 'i' adjacent to city 'j' on the tour?
            vars = m.addVars(dist.keys(), obj=dist, vtype=GRB.BINARY, name='x')

            # Symmetric direction: Copy the object
            for i, j in vars.keys():
                vars[j, i] = vars[i, j]  # edge in opposite direction

            # Constraints: two edges incident to each city
            cons = m.addConstrs(vars.sum(t, '*') == 2 for t in blocked_targets[qn])
            m._vars = vars
            m.Params.lazyConstraints = 1
            m.optimize(subtourelim)

            # Retrieve solution

            vals = m.getAttr('x', vars)
            selected = gp.tuplelist((i, j) for i, j in vals.keys() if vals[i, j] > 0.5)

            tour = subtour(selected)
            assert len(tour) == len(blocked_targets[qn])

            for i in range(len(tour)):
                #Stitch targets from each qn together
                ordered_requests.append(tour[i])
            
        #Turn all the nights targets into starlists   
        write_starlist(all_targets_frame,ordered_requests,condition,current_day)

def semester_schedule(observers_sheet,twilight_times,allocated_nights,marked_scripts,current_day):

    ############Load Observer information and data files, Pandas is an easy way to manage this############
    keck = apl.Observer.at_site('W. M. Keck Observatory')

    #Retrieve the night allocations as csv from jump-config, drop the RM observation nights not part of Community Cadence
    obs_plan = pd.read_csv(allocated_nights)
    obs_plan = obs_plan.drop([25,28,35])

    #Retrieve the generated twilight times information and the HIRES observers sheet as csv's
    twilight_frame = pd.read_csv(twilight_times, parse_dates=True, index_col=0)
    all_targets_frame = pd.read_csv(observers_sheet)
    
    #Turn the observers sheet into a dataframe that is simplified and useful for caluculations. This assumes the 
    #sheet is left downloaded directly
    all_targets_frame = all_targets_frame.drop([0,1])


    ############Compute coordinates in decimal RA hours and decimal DEC deg############
    all_targets_frame['ra'] = all_targets_frame['RAH'] + (1/60)*all_targets_frame['RAM'] + (1/3600)*all_targets_frame['RAS']
    all_targets_frame['dec'] = all_targets_frame['DECD'] + (1/60)*all_targets_frame['DECM'] + (1/3600)*all_targets_frame['DECS']

    #For triple shots, we assume the total exposure time is triple the listed nominal for our calculations
    for index,row in all_targets_frame.iterrows():
        if row['N_obs'] == 3.0:
            all_targets_frame.at[index,'T_exp(sec)'] = row['T_exp(sec)'] * 3
    all_targets_frame['discretized_duration'] = all_targets_frame['T_exp(sec)'].apply(discretize)

    #PI's may choose to input a minimum cadence for each target. If not, one will be generated
    for index,row in all_targets_frame.iterrows():
        if type(row['Cadence']) == str:
            all_targets_frame.at[index,'Cadence'] = int(row['Cadence'])
    all_targets_frame['N_obs(full_semester)'] = all_targets_frame['N_obs(full_semester)'].astype('int')


    #Remove the RM targets corresponding to the RM nights
    all_targets_frame = all_targets_frame[all_targets_frame['Starname'] != 'HIP41378']
    all_targets_frame = all_targets_frame[all_targets_frame['Starname'] != 'TIC256722647']
    all_targets_frame = all_targets_frame[all_targets_frame['Starname'] != 'TIC419523962']

    #Remove the duplicated targets from our list that are calibration shots not relevant to community cadence
    dup = all_targets_frame[all_targets_frame.duplicated(subset='Starname',keep=False)]
    all_targets_frame = all_targets_frame.drop(dup[dup['N_obs(full_semester)'] == 1].index.tolist())
    all_targets_frame = all_targets_frame[all_targets_frame['Done'] != 'done']

    #Finally assign a unique identifier to each request to be used in lookup tables
    all_targets_frame.reset_index(inplace=True,drop=True)
    all_targets_frame['request_number'] = all_targets_frame.index.tolist()

    #Create simplified 'plan' dataframe of quarter nights with columns date, start, and stop
    #These are the discretized time 'buckets' that the top level optimization fills
    sta = []
    sto = []
    dat = []
    for index,row in obs_plan.iterrows():
        num = int((row['stop'] - row['start'])/.25)
        n = row['start']
        for i in range(num):
            sta.append(n + i*.25)
            sto.append(n + (i+1)*.25)
            dat.append(row['Date'])
    plan_dict = {'Date' : dat,'start' : sta, 'stop' : sto}
    plan = pd.DataFrame.from_dict(plan_dict)

    #Call the interval calculator function
    plan = calculate_intervals(plan,twilight_frame)

    #Create a list of reservations and the minimum separation for each target
    reservations,min_separations = generate_reservation_list(all_targets_frame,plan,twilight_frame,keck)

    #Group targets by program in case we want to use statistics or equalizing algorithm
    program_dict = defaultdict(list)
    for code in all_targets_frame['Program code'].unique():

        #'Ex' is the identifier for targets we schedule outside of community cadence, not an actual internal program
        if code != 'Ex':
            requests = all_targets_frame[all_targets_frame['Program code'] == code]['request_number'].tolist()
            program_dict[code] = requests

    #Formally make lists of targets,slots,and observations to communicate to Gurobi
    target_ids = all_targets_frame['request_number'].tolist()
    slots = plan.index.tolist()
    Nobs = all_targets_frame['N_obs(full_semester)'].tolist()

    #Turn our reservation list into a dictionary for when we want to easily access reservations by target
    reservation_dict = defaultdict(list)
    for targ, slot in reservations:
        reservation_dict[targ].append(slot)

    ############Process completed observations############
    observed_dict = process_scripts(all_targets_frame,plan,marked_scripts)

    #Remove reservations that occur before the current day unless actually took place
    starting_slot = plan[plan['Date'] == current_day].index[0]
    for targ in reservation_dict.keys():
        reservation_dict[targ] = [item for item in reservation_dict[targ] 
                            if item >= starting_slot or item in observed_dict[targ]]

    #Account for observations that took place that we may not have calculated to be possible
    #i.e the observer caught a target that was only available for a brief moment that night
    for targ in observed_dict.keys():
        for slot in observed_dict[targ]:
            if slot not in reservation_dict[targ]:
                reservation_dict[targ].append(slot)

    ############Account for weather losses############
    np.random.seed(20)
    day_buffer = 3
    ind = obs_plan[obs_plan['Date'] == current_day].index.values[0] + day_buffer

    #Move forward a bit in time before sampling out quarter nights
    protected = plan[plan['Date'] == obs_plan.loc[ind,'Date']].index.values[0]
    weathered_slots = []

    #Randomly remove 30% of total quarter nights to force the optimizer to work harder
    for slot in slots[protected:]:
        if np.random.randint(10) < 3:
            weathered_slots.append(slot)
    
    #Remove possible reservations from sampled out nights
    for targ in reservation_dict.keys():
        reservation_dict[targ] = [item for item in reservation_dict[targ] if item not in weathered_slots]

    #Cadenced observations are much higher priority when it comes to placement in the semester
    priority_dict = defaultdict(int)
    for targ in all_targets_frame['request_number'].tolist():
        if Nobs[targ] > 1:
            priority_dict[targ] = 10
        else:
            priority_dict[targ] = 1

    #Generate four lists for various conditions for observers to bounce between
    for condition_type in ['nominal','weathered','poor']:

        #Initialize our Gurobi Model
        m = Model('Semester_Plan')

        #Determine the size of our qn buckets
        interval_dict = defaultdict(int)
        for index, row in plan.iterrows():
            interval_dict[index] = (int(np.round((Time(row['qn_stop'],format='jd') - 
                                            Time(row['qn_start'],format='jd')).jd * 24 * 60,0)))

        #Group together our community cadence targets for the optimizer
        cadenced_obs = all_targets_frame[(all_targets_frame['N_obs(full_semester)'] > 1) 
                                    & (all_targets_frame['Program code'] != 'Ex')].index.tolist()


        ############Variable format adopted by ZTF's scheduler############

        #yrt = 1 if target r is scheduled to quarter night t  
        yrt = m.addVars(target_ids,slots, vtype = GRB.BINARY, name = 'Yrt')
        for r in target_ids:
            for t in slots:
                if t not in reservation_dict[r]:
                    m.addConstr((yrt[r,t] == 0),'constr_observable_{}_{}'.format(r,t))

        #If the observers make the mistake of overobserving, we need to increase Nobs to avoid model breakdown
        #by constraint 'constr_num_obs'
        for r in observed_dict.keys():
            for t in observed_dict[r]:
                m.addConstr(yrt[r,t] == 1,'already_sched_{}_{}'.format(r,t))
            if len(observed_dict[r]) > Nobs[r]:
                Nobs[r] = len(observed_dict[r])

        #Less than or equal to sum of durations assigned to each quarter night
        constr_lim_dur = m.addConstrs(
                (gp.quicksum(yrt[r,t] * all_targets_frame.loc[r,'discretized_duration'] 
                            for r in target_ids) <= interval_dict[t] for t in slots[starting_slot:]),'constr_lim_dur')

        #Don't overschedule the number of observations for a target to be greater than those requested
        constr_lim_obs = m.addConstrs(
                (yrt.sum(r,'*') <= Nobs[r] for r in target_ids), 'constr_num_obs')
        
        #Special rules will be applied to the observing night we want to retrieve a schedule for
        next_slots = plan[plan['Date'] == current_day].index.values

        #Force the next night to contain an amount of stars necessary for different conditions
        if condition_type == 'nominal':
            #These bound the size of the upcoming nights 'bin'
            lb = 0.9
            ub = 1.0
            mag_lim = np.inf
            outpfile = '2022B_nominal.csv'
        if condition_type == 'weathered':
            lb = 0.6
            ub = 0.8
            mag_lim = 12
            outpfile = '2022B_weathered.csv'
        if condition_type == 'poor':
            lb = 0.3
            ub = 0.5
            mag_lim = 11
            outpfile = '2022B_poor.csv'
        fill_above = m.addConstr((gp.quicksum(yrt[r,t] * all_targets_frame.loc[r,'discretized_duration'] 
                            for r in target_ids for t in next_slots) >= 
                                lb * sum(interval_dict[t] for t in next_slots))
                                ,'constr_fill_above')
        fill_below = m.addConstr((gp.quicksum(yrt[r,t] * all_targets_frame.loc[r,'discretized_duration'] 
                            for r in target_ids for t in next_slots) <= 
                                ub * sum(interval_dict[t] for t in next_slots))
                                ,'constr_fill_below')
        
        #Magnitude limitation can be turned on using this code
        #for r in all_targets_frame[all_targets_frame['Vmag'] > mag_lim].index.values:
        #    m.addConstr(gp.quicksum(yrt[r,t] for t in next_slots) == 0,
        #                    'constr_lim_magnitude_{}'.format(r))


        #yrtt = 1 if target r is scheduled to both slot t and t2. This is how we build a cadence framework
        yrtt = m.addVars(cadenced_obs,slots,slots,vtype=GRB.BINARY)
        for r in cadenced_obs:
            for t in reservation_dict[r]:
                for t2 in reservation_dict[r]:
                    if ((t2 < t) and (t2 >= starting_slot or t >= starting_slot)):
                        m.addGenConstrAnd(yrtt[r,t,t2],[yrt[r,t],yrt[r,t2]],
                                    "slotdiff_and_{}_{}_{}".format(r,t,t2))
                        #Variabled activated if and only if both corresponding slots are occupied by the target
        
        #Assign a day duration between any two given quarter nights
        dtdict = defaultdict(list)
        for t in slots:
            for t2 in slots:
                if ((t2 <= t) and (t2 >= starting_slot or t >= starting_slot)): 

                    #simplified to days between the date the quarter nights occur
                    dt = np.abs(int((Time(plan.loc[t2,'Date'],format='iso',scale='utc')
                                    -Time(plan.loc[t,'Date'],format='iso',scale='utc')).jd))
                    dtdict[dt].append((t,t2))

        #yrdt = 1 if target r is scheduled to two slots with interval dt between them           
        yrdt = m.addVars(cadenced_obs,dtdict.keys(), vtype=GRB.BINARY)
        for r in cadenced_obs:
            for dt in dtdict.keys():

                #This is activated for ANY two slots that produce this difference dt
                m.addGenConstrOr(yrdt[r,dt], 
                [yrtt[r,t,t2] for (t,t2) in dtdict[dt]],
                    "slot_dt_indicator_{}_{}".format(r,dt))

        #These constraints are built on interface from the human schedulers
        #They may enforce that a target must or must not be scheduled in the next night
        dont_sched = all_targets_frame[(all_targets_frame['Include'] == 'N') | 
                                        (all_targets_frame['Include'] == 'n')].index.tolist()
        dont_include = m.addConstrs((yrt[r,t] == 0 for r in dont_sched for t in next_slots)
                                    ,'dont_sched_{}_{}'.format(r,t))
        force_sched = all_targets_frame[(all_targets_frame['Include'] == 'Y') | 
                                        (all_targets_frame['Include'] == 'y')].index.tolist()
        for r in force_sched:

            #Require that the request can go through only if it is also observable
            #TO DO: find a solution for communicating this error
            '''if can_force(next_slots,reservation_dict[r]):
                m.addConstr((gp.quicksum(yrt[r,t] for t in next_slots) == 1),
                            'forced_sched_{}_{}'.format(r,next_slots[0]))'''

        #This activates the minimum slot separation constraint
        relax_coeff = 1
        constr_min_slotsep = m.addConstrs((yrdt[r,dt] == 0 for r in cadenced_obs for dt in dtdict.keys() 
                                       if dt < relaxed_minimum(min_separations[r],relax_coeff)), 'constr_min_slot_sep')


        #We can also choose to enforce program equality. Research suggests this is often harmful to the end result,
        #for now we've set it to false
        equalize_programs = False
        if equalize_programs:
            equity_list = [item for item in program_dict.keys()]
            yp = m.addVars(program_dict.keys(),vtype=GRB.INTEGER)
            ytot = m.addVar(vtype=GRB.INTEGER)
            m.addConstrs((yp[program] == gp.quicksum(yrt[r,t] for r in program_dict[program] for t in reservation_dict[r])
                        for program in program_dict.keys()),'N_per_program')
            m.addConstr((ytot == gp.quicksum(yrt[r,t] for r,t in reservations)),'N_total')
            dev = m.addVars(program_dict.keys(),lb=-np.inf,vtype=GRB.CONTINUOUS)
            abs_dev = m.addVars(program_dict.keys(),vtype=GRB.CONTINUOUS)
            t_divisor = sum(Nobs)
            for program in program_dict.keys():
                p_divisor = sum(Nobs[r] for r in program_dict[program])
                m.addConstr((dev[program] == ((yp[program]/p_divisor) - (ytot/t_divisor))),'calc_dev')
                m.addGenConstrAbs(abs_dev[program],dev[program],'dev_to_abs')
        
        #These scalars should be adjusted through trial and error, and are dependent on model size.
        cadence_scalar = 1/800
        program_scalar = 600

        #Reward more observations, closer cadence to minimum, and also program equity
        if equalize_programs:
            m.setObjective(gp.quicksum(yrt[r,t] * priority_dict[r] for r,t in reservations)
                    - cadence_scalar * gp.quicksum(yrdt[r,dt] * dt * 1/(Nobs[r]) for r in cadenced_obs for dt in dtdict.keys())
                    - program_scalar * gp.quicksum(abs_dev[program] for program in equity_list)
                    , GRB.MAXIMIZE)

        else:
            m.setObjective(gp.quicksum(yrt[r,t] for r,t in reservations)
                    - cadence_scalar * gp.quicksum(yrdt[r,dt] * dt * 1/(Nobs[r]) for r in cadenced_obs for dt in dtdict.keys())
                    , GRB.MAXIMIZE)

        #Optimization almost always complete or plateaued within 6 minutes
        m.Params.TimeLimit = 360
        m.update()
        m.optimize()

        #Often, for good conditions, with our hard cadence cap it may not be easy to schedule a full night
        while m.Status == GRB.INFEASIBLE:
            #Remove the constraint so it can be redone
            m.remove(constr_min_slotsep)
            relax_coeff -= 0.1
            constr_min_slotsep = m.addConstrs((yrdt[r,dt] == 0 for r in cadenced_obs for dt in dtdict.keys() 
                                        if dt < relaxed_minimum(min_separations[r],relax_coeff)), 'constr_min_slot_sep')

            #It's possible the solve issue isn't with our fill constraints
            #If so, this allows us to view our limiting constraints that break the model
            if relax_coeff <= 0:
                print('Issue with solve for condition {}. Proceeding to IIS computation'.format(condition_type))
                break
            m.update()
            m.optimize()
        
        #Search out model issues
        if m.Status == GRB.INFEASIBLE:
            conflicting_constraints(m,300)
        
        if m.Status != GRB.INFEASIBLE:
            #Create lists of target ids that correspond to each quarter night slot
            scheduled_targets = []
            for v in yrt.values():
                #If decision variable = 1, append the id to that slot
                if v.X == 1:
                    #Perhaps theres a better way than parsing the names, I haven't found it!
                    scheduled_targets.append(v.VarName[4:][:-1].split(','))

            unordered_times = []

            #Reorder the quarter nights
            for i in range(len(scheduled_targets)):
                unordered_times.append(int(scheduled_targets[i][1]))
            order = np.argsort(unordered_times)
            scheduled_targets = [scheduled_targets[i] for i in order]
            
            #Create nightly lists of star requests to be assigned there
            starlists = []
            for i in range(len(plan)):
                s = []
                for j in range(len(scheduled_targets)):
                    if int(scheduled_targets[j][1]) == i:
                        obs = int(scheduled_targets[j][0])
                        s.append(all_targets_frame.loc[obs,'request_number'])
                starlists.append(s)


            #Write four individual files. I like to print them so it's easy to inspect what the top level is producing
            with open(outpfile, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerows(starlists)
    
    salesman_scheduler(all_targets_frame,plan,current_day)
    
    
#Input the filepaths for the relevant data and specify the current date for calculation
semester_schedule('Data/HIRES_observers - PI-requests2022B.csv','Data/twilight_times.csv',
                            'Data/hires_schedule_2022B.csv','Data/MarkedScripts','2022-10-04')


