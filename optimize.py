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
from itertools import permutations
import csv
import logging
import plotting
import accounting


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

def generate_reservation_list(all_targets_frame,plan,twilight_frame):
    keck = apl.Observer.at_site('W. M. Keck Observatory')

    #Astroplan observers don't interface well with astropy functions
    keckapy = apy.coordinates.EarthLocation.of_site('Keck Observatory')

    dates = plan.Date.to_list()
    
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
    else_min_alt = 25
    #min_alt = 40
    #else_min_alt = 40
    max_alt = 85

    logger.info('Generating Reservation List...')

    target_list = []
    ra_dec_list = []
    for index,row in all_targets_frame.iterrows():
        name = row['Starname']
        ra = row['ra']
        dec = row['dec']
        dur = row['discretized_duration']

        #Standard astroplan procedure
        coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
        target = apl.FixedTarget(name=name, coord=coords)
        target_list.append(target)

        #Convert everything to radians
        ra_dec_list.append((ra*np.pi/12,dec*np.pi/180))

        if math.isnan(row['Cadence']):
            min_separations.append(1)
            row['Cadence'] = 1
        else:
            min_separations.append(row['Cadence'])

    AZ = keck.altaz(t, target_list, grid_times_targets=True)
    observability_matrix = []
    for i in range(len(AZ)):
        alt=AZ[i].alt.deg
        az=AZ[i].az.deg
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
        observability_matrix.append(a)

    reservation_matrix = np.zeros((len(all_targets_frame),len(plan)))
    for ind,r in plan.iterrows():
        date = dates[ind]
        day_diff = (Time(date)-start).value

        #This shift is easy to compute and keeps us accurate within a minute at any time
        offset = math.floor(day_diff/15)
        shift = day_diff * 4 - offset

        #Need the time list for calculating the moon observability array
        moon = apy.coordinates.get_moon(Time(r['qn_start'],format='jd'), keckapy)

        #Find the starting/ending minute for each interval in our plan to use as accessibility bounds
        b = int(np.round((Time(r['qn_start'],format='jd')-
            Time(date,format='iso',scale='utc')).jd * 24*60,0) + shift)
        c = int(np.round((Time(r['qn_stop'],format='jd')-
            Time(date,format='iso',scale='utc')).jd * 24*60,0) + 1 + shift)

        for i in range(len(observability_matrix)):
            a = observability_matrix[i]
            if c > len(a):
                one = a[b:]
                remainder = c % len(a)
                two = a[:remainder]
                d = np.concatenate((one,two))
            else:
                d = a[b:c]
            if (len(d) - np.bincount(d)[0]) >= len(d)/2 and moon_safe(moon,ra_dec_list[i]):
                reservation_matrix[i,ind] = 1


    logger.info('Reservations done')
    
    return reservation_matrix,min_separations

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

    #Finally assign the rounded beginning and ending of each quarter night in the plan
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
    return np.round(seconds/60,1)

def can_force(forced_slots,reserved_slots):
    for x in forced_slots:
        for y in reserved_slots:
            if x == y:
                return True
    return False

def relaxed_minimum(minimum,relaxation):
    return math.ceil(minimum * relaxation)

def process_scripts(all_targets_frame,plan,marked_scripts,current_day):

    logger.info('Processing Scripts...')

    import os
    observed_dict = defaultdict(list)

    #The Marked Scipts directory includes the lines from the previous nights in the semester organized into files
    #titled by date. i.e 'marked_scripts/2022-08-07.txt'
    #Note: these files should include only targets for that night, and none of the text/notes that follow
    directory = marked_scripts
    dates = []
    for filename in os.listdir(directory):
        
        targ_per_script = 0
        fpath = os.path.join(directory, filename)
        with open(fpath) as f:

            #Read in each processed script line
            lines = f.readlines()

            #Parse the corresponding date
            day = f.name[:-4][-10:]
            dates.append(day)

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
                    targ_per_script += 1
                    observed_dict[targ[0]].append(slot)

        if targ_per_script == 0:
            logger.warning('No targets recognized in {}'.format(filename))

    current_slot = plan[plan['Date'] == current_day].index[0]

    for script_date in list(dict.fromkeys(plan[:current_slot].Date.tolist())):
        if script_date not in dates:
            logger.warning('Script not found for {}'.format(script_date))


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

def get_ra_dec(all_targets_frame,id_num):
    ra = all_targets_frame.loc[id_num,'ra']
    dec = all_targets_frame.loc[id_num,'dec']
    
    return ra,dec

def write_starlist(frame,requests,condition,current_day):

    #Retrieve columns from observers sheet that are relevant to the starlist
    columns = ['Starname','RAH','RAM','RAS','DECD','DECM','DECS','epoch','vmag=','Vmag',
        'T_exp(sec)','T_max(sec)','Exp_meter_counts(k)','decker','N_obs','Iodine',
        'Priority','Program code','Telescope Propsoal Code','Comment']

    formatted_frame = frame.loc[requests][columns]

    logger.info('Writing starlist to {}_{}.txt'.format(current_day,condition))
    lines = []
    for index,row in formatted_frame.iterrows():

        #Just a bunch of string formatting. This prints standard starlists as ordered by the salesman optimization
        namestring = ' '*(16-len(row['Starname'])) + row['Starname']
        
        rastring = ('0'*(2-len(str(int(row['RAH'])))) + str(int(row['RAH'])) + ' '
                        + '0'*(2-len(str(int(row['RAM'])))) + str(int(row['RAM'])) + ' '
                            + '0'*(4-len(str(row['RAS']))) + str(row['RAS']))
        
        starter = '+'
        if row['DECD'][0] == '-' or row['DECD'][0] == '−':
            starter = '-'
            decstring = (starter + '0'*(2-len(str(abs(int(row['DECD'][1:]))))) + str(abs(int(row['DECD'][1:]))) + ' '
                            + '0'*(2-len(str(int(row['DECM'])))) + str(int(row['DECM'])) + ' '
                                + '0'*(4-len(str(row['DECS']))) + str(int(row['DECS'])))
        else:
            decstring = (starter + '0'*(2-len(str(abs(int(row['DECD']))))) + str(abs(int(row['DECD']))) + ' '
                            + '0'*(2-len(str(int(row['DECM'])))) + str(int(row['DECM'])) + ' '
                                + '0'*(4-len(str(row['DECS']))) + str(int(row['DECS'])))
        
        magstring = row['vmag='] + str(row['Vmag']) + ' '*(4-len(str(row['Vmag'])))
        
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

def get_alt_az(frame,id_num,time,observatory):
    ra,dec = get_ra_dec(frame,id_num)
    
    coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
    target = apl.FixedTarget(name=id_num, coord=coords)
    AZ = observatory.altaz(time, target)
    
    return AZ.alt.deg,AZ.az.deg

def moon_safe(moon,target_tuple):
    ang_dist = apy.coordinates.angular_separation(moon.ra.rad,moon.dec.rad,
                                                  target_tuple[0],target_tuple[1])
    if ang_dist*180/(np.pi) >= 30:
        return True
    else:
        return False
            
def salesman_scheduler(all_targets_frame,plan,current_day,output_flag,plot_results):
    
    keck = apl.Observer.at_site('W. M. Keck Observatory')
    for condition in ['nominal']:

        logger.info('Building traveling salesman for condition: {}'.format(condition))

        #Retrieve the conditions from respective csv's
        file = open("2023A_{}.csv".format(condition), "r")
        blocked_targets = list(csv.reader(file, delimiter=","))
        file.close()
        for i in range(len(blocked_targets)):
            blocked_targets[i] = blocked_targets[i][1:]
            for j in range(len(blocked_targets[i])):
                blocked_targets[i][j] = int(blocked_targets[i][j])

        ordered_requests = []       
        
        #Plot folder
        plotpath = '{}_plots'.format(current_day)
        
        #The path for each quarter night is computed independently
        for qn in plan[plan['Date'] == current_day].index.tolist(): 

            nightly_targets = blocked_targets[qn]

            #Find reason for duplicates
            nightly_targets = list(set(nightly_targets))

            ind_to_id = defaultdict(int)
            i = 1
            for targ in nightly_targets:
                ind_to_id[i] = targ
                i += 1

            start = Time(plan.loc[qn,'qn_start'],format='jd')
            end = Time(plan.loc[qn,'qn_stop'],format='jd')
            dur = np.round((end-start).jd*24*60,0)*60
            stop = start + TimeDelta(dur,format='sec')
            step = TimeDelta(60,format='sec')
            t = np.arange(start.jd,stop.jd,step.jd)
            t = Time(t,format='jd')

            R = len(nightly_targets)+2
            T = math.ceil(dur/(60*15))
            slots = range(T)

            min_az = 5.3
            max_az = 146.2
            min_alt = 33.3
            else_min_alt = 25
            #min_alt = 40
            #else_min_alt = 40
            max_alt = 90

            e_i = []
            l_i = []
            s_i = []
            for i in range(R):
                if i == 0 or i == R-1:
                    e_i.append(0)
                    l_i.append((stop-t[0]).jd*24*60)
                    s_i.append(0)
                else:
                    targ = ind_to_id[i]
                    exp = all_targets_frame.loc[targ,'discretized_duration']
                    s_i.append(exp)
                    ra,dec = get_ra_dec(all_targets_frame,targ)
                    coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
                    target = apl.FixedTarget(name=targ, coord=coords)
                    AZ = keck.altaz(t, target)
                    alt=AZ.alt.deg
                    az=AZ.az.deg
                    deck = np.where((az >= min_az) & (az <= max_az))
                    deck_height = np.where((alt <= max_alt) & (alt >= min_alt))
                    first = np.intersect1d(deck,deck_height)
                    not_deck_1 = np.where((az < min_az))
                    not_deck_2 = np.where((az > max_az))
                    not_deck_height = np.where((alt <= max_alt) & (alt >= else_min_alt))
                    second = np.intersect1d(not_deck_1,not_deck_height)
                    third = np.intersect1d(not_deck_2,not_deck_height)
                    good = np.concatenate((first,second,third))
                    if len(good > 0):
                        e_i.append(np.round(((t[good[0]]+TimeDelta(s_i[i]*60,format='sec'))-t[0]).jd*24*60,1))
                        if t[good[-1]].jd < stop.jd:
                            l_i.append(np.round((t[good[-1]]-t[0]).jd*24*60,1))
                        elif t[good[-1]].jd >= stop.jd:
                            l_i.append(np.round((stop-t[0]).jd*24*60,1))
                    else:
                        print('Error, target {} does not meet observability requirements'.format(all_targets_frame.loc[ind_to_id[i],'Starname']))
            
            def to_wrap_frame(angle):
                angle += 90
                if angle >= 360:
                    angle -= 360
                return angle
            
            start_slots = Time(np.linspace(start.jd,stop.jd,T,endpoint=False),format='jd')
            #Move to middle of slot for coordinates
            shift = (end.jd-start.jd)/(2*T)
            start_slots += shift
            logger.info('Calculating Distance Tensor for qn {}'.format(qn))
            
            coordinate_matrix = []
            coordinate_matrix.append(np.zeros(len(start_slots)))
            for i in range(len(nightly_targets)):
                id_num = nightly_targets[i]
                ra,dec = get_ra_dec(all_targets_frame,id_num)
                coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
                target = apl.FixedTarget(name=id_num, coord=coords)
                coordinate_matrix.append(keck.altaz(start_slots,target).az.deg)
            coordinate_matrix.append(np.zeros(len(start_slots)))
            coordinate_matrix = np.array(coordinate_matrix)

            def distance(targ1,targ2,slot):
                if targ1 == 0 and targ2 == R-1:
                    return 0
                if targ1 == R-1 and targ2 == 0:
                    return 0
                if targ1 == 0 or targ1 == R-1:
                    #return to_wrap_frame(coordinate_matrix[targ2,slot])
                    return 0
                if targ2 == 0 or targ2 == R-1:
                    #return to_wrap_frame(coordinate_matrix[targ1,slot])
                    return 0
                t1 = to_wrap_frame(coordinate_matrix[targ1,slot])
                t2 = to_wrap_frame(coordinate_matrix[targ2,slot])

                #diff = (t1[0]-t2[0], t1[1]-t2[1])
                #return math.sqrt(diff[0]*diff[0]+diff[1]*diff[1])
                return np.abs(t1-t2)

            dist = defaultdict(float)
            for i in range(len(start_slots)):
                for targ1,targ2 in permutations(range(R),2):
                    dist[(targ1,targ2,i)] = np.round(distance(targ1,targ2,i)/60,1)

            #Time windows
            w = []
            for m in range(T):
                w.append(np.round((start_slots[m] - start_slots[0]).jd*24*60,1))
            w.append(np.round((stop - start_slots[0]).jd*24*60,1))

            ####Model and Variables####
            Mod = Model('TDTSP')
            Mod.Params.OutputFlag = output_flag

            yi = Mod.addVars(range(R),vtype=GRB.BINARY,name='yi')
            xijm = Mod.addVars(range(R),range(R),range(T),vtype=GRB.BINARY,name='xijm')
            tijm = Mod.addVars(range(R),range(R),range(T),vtype=GRB.CONTINUOUS,name='tijm')
            ti = Mod.addVars(range(R),vtype=GRB.CONTINUOUS,lb=0,name='ti')

            tijmdef = Mod.addConstrs((ti[i] == gp.quicksum(tijm[i,j,m] for j in range(R)[1:] for m in range(T)) 
                                for i in range(R)[:-1]),'tijm_def')
            start_origin = Mod.addConstr(gp.quicksum(xijm[0,j,m] for j in range(R) for m in range(T)) == 1,
                               'start_origin')
            end_origin = Mod.addConstr(gp.quicksum(xijm[i,R-1,m] for i in range(R) for m in range(T)) == 1,
                               'end_origin')
            visit_once = Mod.addConstrs((gp.quicksum(xijm[i,j,m] for i in range(R)[:-1] for m in range(T)) == yi[j]
                             for j in range(R)[1:]), 'visit_once')
            flow_constr = Mod.addConstrs(((gp.quicksum(xijm[i,k,m] for i in range(R)[:-1] for m in range(T))
                               - gp.quicksum(xijm[k,j,m] for j in range(R)[1:] for m in range(T)) == 0)
                              for k in range(R)[:-1][1:]), 'flow_constr')
            exp_constr = Mod.addConstrs((ti[j] >= tijm[i,j,m] + (dist[(i,j,m)] + s_i[j])*xijm[i,j,m] 
                             for i in range(R)[:-1] for j in range(R)[1:] for m in range(T))
                          , 'exp_constr')
            t_min = Mod.addConstrs(((tijm[i,j,m] >= w[m]*xijm[i,j,m]) for j in range(R) for m in range(T)
                             for i in range(R)),'t_min')
            t_max = Mod.addConstrs((tijm[i,j,m] <= w[m+1]*xijm[i,j,m] for j in range(R) for m in range(T) 
                             for i in range(R)),'t_max')
            rise_constr = Mod.addConstrs((ti[i] >= e_i[i]*yi[i] for i in range(R)),'rise_constr')
            set_constr = Mod.addConstrs((ti[i] <= l_i[i]*yi[i] for i in range(R)),'set_constr')
            priority_param = 50
            slew_param = 1/100

            Mod.setObjective(priority_param*gp.quicksum(yi[j] for j in range(R)) 
                             -slew_param *gp.quicksum(dist[(i,j,m)]*xijm[i,j,m] for i in range(R)[:-1] 
                                           for j in range(R)[1:] for m in range(T))
                                ,GRB.MAXIMIZE)
            logger.info('Solving TDTSP for qn {}'.format(qn))
            Mod.params.TimeLimit = 300
            Mod.params.MIPGap = .0005
            Mod.update()
            Mod.optimize()

            scheduled_targets = []
            for i in range(R)[1:][:-1]:
                if np.round(yi[i].X,0) == 1:
                    v = ti[i]
                    scheduled_targets.append((int(v.VarName[3:-1]),int(np.round(v.X))))

            start_times = []
            unordered_times = []

            for i in range(len(scheduled_targets)):
                unordered_times.append(int(scheduled_targets[i][1]))
            order = np.argsort(unordered_times)
            scheduled_targets = [scheduled_targets[i] for i in order]
            

            if plot_results:
                logger.info('Plotting qn {}'.format(qn))

                obs_time=[]
                az_path=[]
                alt_path=[]
                names=[]
                targ_list=[]
                for pair in scheduled_targets:
                    targ=pair[0]
                    time=t[pair[1]]
                    exp=s_i[targ]
                    targ_list.append(targ)
                    names.append(all_targets_frame.loc[ind_to_id[targ],'Starname'])
                    time1=time-TimeDelta(60*exp,format='sec')
                    obs_time.append(time1.jd)
                    obs_time.append(time.jd)
                    az_path.append(get_alt_az(all_targets_frame,ind_to_id[targ],time1,keck)[1])
                    az_path.append(get_alt_az(all_targets_frame,ind_to_id[targ],time,keck)[1])
                    alt_path.append(get_alt_az(all_targets_frame,ind_to_id[targ],time1,keck)[0])
                    alt_path.append(get_alt_az(all_targets_frame,ind_to_id[targ],time,keck)[0])
                    
                time_array = []
                for i in range(len(obs_time)):
                    if i % 2 == 0:
                        time_array.append((obs_time[i] - obs_time[0])*1440)

                time_counter = Time(obs_time[0],format='jd')
                time_strings = []
                tel_az = []
                tel_zen = []
                observed_at_time = []
                while time_counter.jd < obs_time[-1]:
                    time_strings.append(time_counter.isot)
                    ind = np.flatnonzero(time_array <= (time_counter.jd - obs_time[0])*1440)[-1]
                    req = ind_to_id[targ_list[ind]]
                    observed_at_time.append(ind)
                    tel_az.append(get_alt_az(all_targets_frame,req,time_counter,keck)[1] * (np.pi/180))
                    tel_zen.append(90-get_alt_az(all_targets_frame,req,time_counter,keck)[0])
                    time_counter += TimeDelta(30,format='sec')

                time_counter = Time(obs_time[0],format='jd')
                total_azimuth_list = []
                total_zenith_list = []
                while time_counter.jd < obs_time[-1]:
                    azimuth_list = []
                    zenith_list = []
                    for targ in targ_list:
                        azimuth_list.append(get_alt_az(all_targets_frame,ind_to_id[targ],time_counter,keck)[1] * (np.pi/180))
                        zenith_list.append(90-get_alt_az(all_targets_frame,ind_to_id[targ],time_counter,keck)[0])
                    total_azimuth_list.append(azimuth_list)
                    total_zenith_list.append(zenith_list)
                    time_counter += TimeDelta(30,format='sec')
                
                plotting.plot_path_2D(obs_time,az_path,alt_path,names,targ_list,
                                    plotpath,current_day,int(plan.loc[qn,'stop']*4))
                plotting.animate_telescope(time_strings,total_azimuth_list,total_zenith_list,
                                        tel_az,tel_zen,observed_at_time,plotpath,int(plan.loc[qn,'stop']*4))
            
            for pair in scheduled_targets:
                ordered_requests.append(ind_to_id[pair[0]])

        #Turn all the nights targets into starlists   
        write_starlist(all_targets_frame,ordered_requests,condition,current_day)

def semester_schedule(observers_sheet,twilight_times,allocated_nights,marked_scripts,schedule_dates,output_flag,
                                                                                    equalize_programs,plot_results):

    if len(schedule_dates) == 1:
        high_production_mode = False
    if len(schedule_dates) > 1:
        high_production_mode = True
    
    current_day = schedule_dates[0]
    
    ############Load Observer information and data files, Pandas is an easy way to manage this############
    keck = apl.Observer.at_site('W. M. Keck Observatory')

    #Retrieve the night allocations as csv from jump-config, drop the RM observation nights not part of Community Cadence
    obs_plan = pd.read_csv(allocated_nights)
    obs_plan = obs_plan[obs_plan['Date'] != '2022-02-09']
    #obs_plan = obs_plan[obs_plan['Date'] != '2022-06-24']
    #obs_plan = obs_plan[obs_plan['Date'] != '2022-07-06']
    obs_plan.reset_index(inplace=True,drop=True)

    #Retrieve the generated twilight times information and the HIRES observers sheet as csv's
    twilight_frame = pd.read_csv(twilight_times, parse_dates=True, index_col=0)
    all_targets_frame = pd.read_csv(observers_sheet)
    
    #Turn the observers sheet into a dataframe that is simplified and useful for caluculations. This assumes the 
    #sheet is left downloaded directly
    all_targets_frame = all_targets_frame.drop([0,1,2,3,4])


    ############Compute coordinates in decimal RA hours and decimal DEC deg############
    all_targets_frame['ra'] = all_targets_frame['RAH'] + (1/60)*all_targets_frame['RAM'] + (1/3600)*all_targets_frame['RAS']

    #Reading in dec can be troublesome
    all_targets_frame['dec'] = '' 
    for index,row in all_targets_frame.iterrows():
        if row['DECD'][0] == '-' or row['DECD'][0] == '−':
            dec = (int(row['DECD'][1:]) + (1/60)*row['DECM'] + (1/3600)*row['DECS'])
            all_targets_frame.at[index,'dec'] = -np.abs(dec)
        else:
            dec = (int(row['DECD']) + (1/60)*row['DECM'] + (1/3600)*row['DECS'])
            all_targets_frame.at[index,'dec'] = np.abs(dec)
            
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
    all_targets_frame = all_targets_frame[all_targets_frame['Starname'] != 'TIC69997672']
    all_targets_frame = all_targets_frame[all_targets_frame['Starname'] != 'TIC230075120']
    all_targets_frame = all_targets_frame[all_targets_frame['Starname'] != 'TIC139474683']

    #Remove the duplicated targets from our list that are calibration shots not relevant to community cadence
    dup = all_targets_frame[all_targets_frame.duplicated(subset='Starname',keep=False)]
    all_targets_frame = all_targets_frame.drop(dup[dup['N_obs(full_semester)'] == 1].index.tolist())
    all_targets_frame = all_targets_frame[all_targets_frame['Done'] != 'done']
    all_targets_frame = all_targets_frame[all_targets_frame['N_obs(full_semester)'] > 0]

    #Finally assign a unique identifier to each request to be used in lookup tables
    all_targets_frame.reset_index(inplace=True,drop=True)
    all_targets_frame['request_number'] = all_targets_frame.index.tolist()
    
    '''#Targets forcefully scheduled will have different rules
    force_sched = all_targets_frame[(all_targets_frame['Include'] == 'Y') | 
                                        (all_targets_frame['Include'] == 'y')].index.tolist()'''

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
    reservations,min_separations = generate_reservation_list(all_targets_frame,plan,twilight_frame)

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
    for targ in range(len(reservations)):
        for slot in range(len(reservations[targ])):
            if reservations[targ][slot] == 1:
                reservation_dict[targ].append(slot)

    ############Process completed observations############
    observed_dict = process_scripts(all_targets_frame,plan,marked_scripts,current_day)

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
    day_buffer = 4
    if high_production_mode == False:
        ind = obs_plan[obs_plan['Date'] == current_day].index.values[0] + day_buffer
    if high_production_mode == True:
        ind = obs_plan[obs_plan['Date'] == schedule_dates[-1]].index.values[0] + day_buffer
    

    if ind in obs_plan.index.values:
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

    #Generate three lists for various conditions for observers to bounce between
    for condition_type in ['nominal']:

        logger.info('Scheduling semester for condition: {}'.format(condition_type))

        #Initialize our Gurobi Model
        m = Model('Semester_Plan')
        m.Params.OutputFlag = output_flag

        #Determine the size of our qn buckets
        interval_dict = defaultdict(int)
        for index, row in plan.iterrows():
            interval_dict[index] = (int(np.round((Time(row['qn_stop'],format='jd') - 
                                            Time(row['qn_start'],format='jd')).jd * 24 * 60,0)))

        #Group together our community cadence targets for the optimizer
        cadenced_obs = all_targets_frame[(all_targets_frame['N_obs(full_semester)'] > 1) 
                                    & (all_targets_frame['Program code'] != 'Ex')].index.tolist()
        
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
        fill_slots = plan[plan['Date'].isin(schedule_dates)].index.values

        #Force the next night to contain an amount of stars necessary for different conditions
        if condition_type == 'nominal':
            #These bound the size of the upcoming nights 'bin'
            lb = 0.95
            ub = 1.0
            mag_lim = np.inf
            outpfile = '2023A_nominal.csv'
        if condition_type == 'weathered':
            lb = 0.6
            ub = 0.8
            mag_lim = 12
            outpfile = '2023A_weathered.csv'
        if condition_type == 'poor':
            lb = 0.3
            ub = 0.5
            mag_lim = 11
            outpfile = '2023A_poor.csv'
        fill_above = m.addConstrs((gp.quicksum(yrt[r,t] * all_targets_frame.loc[r,'discretized_duration'] 
                            for r in target_ids) >= lb * interval_dict[t] for t in fill_slots)
                                ,'constr_fill_above')
        fill_below = m.addConstrs((gp.quicksum(yrt[r,t] * all_targets_frame.loc[r,'discretized_duration'] 
                            for r in target_ids) <= ub * interval_dict[t] for t in fill_slots)
                                ,'constr_fill_below')
        
        #Magnitude limitation can be turned on using this code
        #for r in all_targets_frame[all_targets_frame['Vmag'] > mag_lim].index.values:
        #    m.addConstr(gp.quicksum(yrt[r,t] for t in next_slots) == 0,
        #                    'constr_lim_magnitude_{}'.format(r))


        #Brtt = 1 if target r is scheduled to both slot t and t2. This is how we build a cadence framework
        Brtt = m.addVars(cadenced_obs,slots,slots,vtype=GRB.BINARY)
        for r in cadenced_obs:
            for t in reservation_dict[r]:
                for t2 in reservation_dict[r]:
                    if ((t2 < t) and (t2 >= starting_slot or t >= starting_slot)):
                        m.addGenConstrAnd(Brtt[r,t,t2],[yrt[r,t],yrt[r,t2]],
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

        #Drdt = 1 if target r is scheduled to two slots with interval dt between them           
        Drdt = m.addVars(cadenced_obs,dtdict.keys(), vtype=GRB.BINARY)
        for r in cadenced_obs:
            for dt in dtdict.keys():

                #This is activated for ANY two slots that produce this difference dt
                m.addGenConstrOr(Drdt[r,dt], 
                [Brtt[r,t,t2] for (t,t2) in dtdict[dt]],
                    "slot_dt_indicator_{}_{}".format(r,dt))

        #These constraints are built on interface from the human schedulers
        #They may enforce that a target must or must not be scheduled in the next night
        dont_sched = all_targets_frame[(all_targets_frame['Include'] == 'N') | 
                                        (all_targets_frame['Include'] == 'n')].index.tolist()
        dont_include = m.addConstrs((yrt[r,t] == 0 for r in dont_sched for t in next_slots)
                                    ,'dont_sched_{}_{}'.format(r,t))
        
        '''for r in force_sched:
            #Require that the request can go through only if it is also observable
            if can_force(next_slots,reservation_dict[r]):
                m.addConstr((gp.quicksum(yrt[r,t] for t in next_slots) == 1),
                            'forced_sched_{}_{}'.format(r,next_slots[0]))'''
                
        #This activates the minimum slot separation constraint
        relax_coeff = 1
        constr_min_slotsep = m.addConstrs((Drdt[r,dt] == 0 for r in cadenced_obs for dt in dtdict.keys() 
                                       if dt < relaxed_minimum(min_separations[r],relax_coeff)), 'constr_min_slot_sep')


        #We can also choose to enforce program equality. Research suggests this is often harmful to the end result,
        #for now we've set it to false
        if equalize_programs:
            equity_list = [item for item in program_dict.keys()]
            yp = m.addVars(program_dict.keys(),vtype=GRB.INTEGER)
            ytot = m.addVar(vtype=GRB.INTEGER)
            m.addConstrs((yp[program] == gp.quicksum(yrt[r,t] for r in program_dict[program] for t in reservation_dict[r])
                        for program in program_dict.keys()),'N_per_program')
            m.addConstr((ytot == gp.quicksum(yrt[r,t] for r in target_ids for t in slots)),'N_total')
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
            m.setObjective(gp.quicksum(yrt[r,t] * priority_dict[r] for r in target_ids for t in slots)
                    - cadence_scalar * gp.quicksum(Drdt[r,dt] * dt * 1/(Nobs[r]) for r in cadenced_obs for dt in dtdict.keys())
                    - program_scalar * gp.quicksum(abs_dev[program] for program in equity_list)
                    , GRB.MAXIMIZE)

        else:
            m.setObjective(gp.quicksum(yrt[r,t] for r in target_ids for t in slots)
                    - cadence_scalar * gp.quicksum(Drdt[r,dt] * dt * 1/(Nobs[r]) for r in cadenced_obs for dt in dtdict.keys())
                    , GRB.MAXIMIZE)

        #Optimization almost always complete or plateaued within 6 minutes
        m.Params.TimeLimit = 300
        m.update()
        m.optimize()

        #Often, for good conditions, with our hard cadence cap it may not be easy to schedule a full night
        while m.Status == GRB.INFEASIBLE:
            #Remove the constraint so it can be redone
            m.remove(constr_min_slotsep)
            relax_coeff -= 0.1
            constr_min_slotsep = m.addConstrs((Drdt[r,dt] == 0 for r in cadenced_obs for dt in dtdict.keys() 
                                        if dt < relaxed_minimum(min_separations[r],relax_coeff)), 'constr_min_slot_sep')

            #It's possible the solve issue isn't with our fill constraints
            #If so, this allows us to view our limiting constraints that break the model
            if relax_coeff <= 0:
                logger.error('Issue with solve for condition {}. Proceeding to IIS computation'.format(condition_type))
                break
            m.update()
            m.optimize()
        
        #Search out model issues
        if m.Status == GRB.INFEASIBLE:
            logger.info('Model remains infeasible. Searching for invalid constraints')
            conflicting_constraints(m,300)
        
        if m.Status != GRB.INFEASIBLE:
            #Create lists of target ids that correspond to each quarter night slot
            scheduled_targets = []
            for v in yrt.values():
                #If decision variable = 1, append the id to that slot
                if np.round(v.X,0) == 1:
                    #Perhaps theres a better way than parsing the names, I haven't found it!
                    scheduled_targets.append(v.VarName[4:][:-1].split(','))

            unordered_times = []

            #Reorder the quarter nights
            for i in range(len(scheduled_targets)):
                unordered_times.append(int(scheduled_targets[i][1]))
            order = np.argsort(unordered_times)
            scheduled_targets = [scheduled_targets[i] for i in order]
            
            targets_observed = defaultdict(list)
            for i in range(len(scheduled_targets)):
                slot = int(scheduled_targets[i][1])
                targets_observed[slot].append(int(scheduled_targets[i][0]))
            
            #Create nightly lists of star requests to be assigned there
            starlists = []
            for i in range(len(plan)):
                time_log = plan.loc[i,'Date']
                s = [time_log]
                for j in range(len(scheduled_targets)):
                    if int(scheduled_targets[j][1]) == i:
                        obs = int(scheduled_targets[j][0])
                        s.append(all_targets_frame.loc[obs,'request_number'])
                starlists.append(s)

            #TODO update accounting logs for high production mode
            #accounting.completion_logs(all_targets_frame,observed_dict,starlists,current_day)

            if plot_results:
                #Plot program CDF's
                logger.info('Plotting Program CDFs')
                if high_production_mode == True:
                    for day in schedule_dates[1:]:
                        plotpath = '{}_plots'.format(day)
                        if not os.path.isdir(plotpath):
                            os.mkdir(plotpath)
                plotpath = '{}_plots'.format(current_day)
                if not os.path.isdir(plotpath):
                    os.mkdir(plotpath)
                plotting.plot_program_cdf(plan,program_dict,targets_observed,Nobs,plotpath,current_day)
                #plotting.plot_one_cdf(plan,program_dict,targets_observed,Nobs,plotpath,current_day,'DG')
                
                #Plot target cadence by program
                logger.info('Plotting Program Cadences')
                plotting.plot_program_cadence(plan,all_targets_frame,twilight_frame,starlists,
                                    min_separations,plotpath,current_day)
                #plotting.plot_cadence_night_resolution(plan,all_targets_frame,twilight_frame,starlists,
                #                     min_separations,plotpath)


            #Write three individual files. Print them so it's easy to inspect what the top level is producing
            logger.info('Writing to {}'.format(outpfile))
            with open(outpfile, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerows(starlists)
    
    if high_production_mode == True:
        for date_to_schedule in schedule_dates:
            salesman_scheduler(all_targets_frame,plan,date_to_schedule,output_flag,plot_results)
    if high_production_mode == False:
        salesman_scheduler(all_targets_frame,plan,current_day,output_flag,plot_results)
