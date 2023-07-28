import argparse

parser = argparse.ArgumentParser(description='Runtime tests for the Traveling Telescope Problem through Gurobi')
parser.add_argument('-t','--num_targs', help='Number of targets to sample')
parser.add_argument('-s','--num_slots', help='Number of slots used to discretize time interval')
parser.add_argument('-d','--sim_type', help='Duration of time interval: Quarter=1, Half=2, or Full=3')
parser.add_argument('-g','--mip_gap', help='Absolute optimality gap in minutes, default is 0',default=0)
parser.add_argument('-l','--time_limit', help='Time limit for Gurobi optimizer, default is none')
parser.add_argument('-o','--gurobi_output',action='store_true', help='Activate Gurobi console outputs',
                default=False)

args = parser.parse_args()
num_targs = int(args.num_targs)
num_slots = int(args.num_slots)
gap = int(args.mip_gap)
sim_num = int(args.sim_type)

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
from itertools import permutations



if sim_num == 1:
    sim_type = 'Quarter'
if sim_num == 2:
    sim_type = 'Half'
if sim_num == 3:
    sim_type = 'Full'

def to_wrap_frame(angle):
    angle += 90
    if angle >= 360:
        angle -= 360
    return angle

def get_ra_dec(all_targets_frame,id_num):
    ra = all_targets_frame.loc[id_num,'ra']
    dec = all_targets_frame.loc[id_num,'dec']
    
    return ra,dec

def discretize(seconds):
    return np.round(seconds/60,1)

def get_alt_az(frame,id_num,time,observatory):
    ra,dec = get_ra_dec(frame,id_num)
    
    coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
    target = apl.FixedTarget(name=id_num, coord=coords)
    AZ = observatory.altaz(time, target)
    
    return AZ.alt.deg,AZ.az.deg

def generate_reservation_list(all_targets_frame,start,stop):

    keck = apl.Observer.at_site('W. M. Keck Observatory')

    step = TimeDelta(60,format='sec')
    t = np.arange(start.jd,stop.jd,step.jd)
    t = Time(t,format='jd')

    #These coordinates designate observability at Keck
    min_az = 5.3
    max_az = 146.2
    min_alt = 33.3
    else_min_alt = 25
    max_alt = 90

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
    
    reservation_matrix = np.zeros(len(all_targets_frame))
    AZ = keck.altaz(t, target_list, grid_times_targets=True)
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
        if (len(a) - np.bincount(a)[0]) >= len(a)/2:
            reservation_matrix[i] = 1
    
    return reservation_matrix

############Load Observer information and data files, Pandas is an easy way to manage this############
keck = apl.Observer.at_site('W. M. Keck Observatory')
cls_frame = pd.read_csv('query-cls-legacy-targets.csv')

output = {'Duration': [], 'NumberTargets': [], 'NumberSlots': [], 'Runtime(sec)': [], 'RealSlewTime(min)': [], 
          'RealMedianSlew(sec)': [], 'TauSlewTime(min)': [], 'TauMedianSlew(sec)': [], 'LogFile': []}
outputframe = pd.DataFrame.from_dict(output)

if sim_type == 'Quarter':
    start = Time('2023-07-06T05:59:00.000',format='isot')
    end = Time('2023-07-06T08:12:00.000',format='isot')
if sim_type == 'Half':
    start = Time('2023-07-06T05:59:00.000',format='isot')
    end = Time('2023-07-06T10:26:00.000',format='isot')
if sim_type == 'Full':
    start = Time('2023-07-06T05:59:00.000',format='isot')
    end = Time('2023-07-06T14:55:00.000',format='isot')
    
duration = int(np.round((end.jd - start.jd) * 24 * 60,0))

cls_frame['discretized_duration'] = [(duration-3*num_targs)/num_targs]*len(cls_frame)
reservations = generate_reservation_list(cls_frame,start,end)
visible_targets = []
for i in range(len(reservations)):
    if reservations[i] == 1:
        visible_targets.append(i)

np.random.seed(7062023)
success = False
while success == False:

    random_sample = np.random.randint(len(visible_targets),size=num_targs)
    while len(random_sample) != len(set(random_sample)):
            random_sample = np.random.randint(len(visible_targets),size=num_targs)
    nightly_targets = [visible_targets[index] for index in random_sample]

    ind_to_id = defaultdict(int)
    id_to_ind = defaultdict(int)
    i = 1
    for targ in nightly_targets:
        ind_to_id[i] = targ
        id_to_ind[targ] = i
        i += 1

    dur = np.round((end.jd-start.jd)*24*60,0)*60
    stop = start + TimeDelta(dur,format='sec')
    step = TimeDelta(60,format='sec')
    t = np.arange(start.jd,stop.jd,step.jd)
    t = Time(t,format='jd')

    R = num_targs + 2
    constrained_targets = [i for i in range(R)[1:-1]]

    T = num_slots
    slots = range(T)

    min_az = 5.3
    max_az = 146.2
    min_alt = 33.3
    else_min_alt = 25
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
            exp = cls_frame.loc[targ,'discretized_duration']
            s_i.append(exp)
            ra,dec = get_ra_dec(cls_frame,targ)
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
            e_i.append(np.round(((t[good[0]]+TimeDelta(s_i[i]*60,format='sec'))-t[0]).jd*24*60,1))
            if t[good[-1]].jd < stop.jd:
                l_i.append(np.round((t[good[-1]]-t[0]).jd*24*60,1))
            elif t[good[-1]].jd >= stop.jd:
                l_i.append(np.round((stop-t[0]).jd*24*60,1))

    start_slots = Time(np.linspace(start.jd,stop.jd,T,endpoint=False),format='jd')
    #Move to middle of slot for coordinates
    shift = (end.jd-start.jd)/(2*T)
    slot_times = start_slots + shift
    coordinate_matrix = []
    coordinate_matrix.append([(0,0) for slot in slot_times])
    for i in range(len(nightly_targets)):
        id_num = nightly_targets[i]
        ra,dec = get_ra_dec(cls_frame,id_num)
        coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
        target = apl.FixedTarget(name=id_num, coord=coords)
        altaz = keck.altaz(slot_times,target)
        coordinate_matrix.append([(item.alt.deg,item.az.deg) for item in altaz])
    coordinate_matrix.append([(0,0) for slot in slot_times])
    coordinate_matrix = np.array(coordinate_matrix)

    def distance(targ1,targ2,slot):
        if targ1 == 0 or targ1 == R-1:
            return 0
        if targ2 == 0 or targ2 == R-1:
            return 0

        alt1 = coordinate_matrix[targ1,slot][0]
        alt2 = coordinate_matrix[targ2,slot][0]                             
        az1 = to_wrap_frame(coordinate_matrix[targ1,slot][1])
        az2 = to_wrap_frame(coordinate_matrix[targ2,slot][1])

        return max(np.abs(alt1-alt2),np.abs(az1-az2))

    dist = defaultdict(float)
    for i in range(len(slot_times)):
        for targ1,targ2 in permutations(range(R),2):
            dist[(targ1,targ2,i)] = np.round(distance(targ1,targ2,i)/60,2)

    w = []
    for m in range(T):
        w.append(np.round((start_slots[m] - start_slots[0]).jd*24*60,1))
    w.append(np.round((stop - start_slots[0]).jd*24*60,1))

    Mod = Model('TDTSP')
    Mod.Params.OutputFlag = args.gurobi_output

    yi = Mod.addVars(range(R),vtype=GRB.BINARY,name='yi')
    force_sched = Mod.addConstrs((yi[i] == 1 for i in range(R)[1:-1]),'force_sched')
    xijm = Mod.addVars(range(R),range(R),range(T),vtype=GRB.BINARY,name='xijm')
    tijm = Mod.addVars(range(R),range(R),range(T),vtype=GRB.CONTINUOUS,name='tijm')
    ti = Mod.addVars(range(R),vtype=GRB.CONTINUOUS,lb=0,name='ti')
    start_origin = Mod.addConstr(gp.quicksum(xijm[0,j,m] for j in range(R)[1:-1] for m in range(T)) == 1,
                                'start_origin')
    end_origin = Mod.addConstr(gp.quicksum(xijm[i,R-1,m] for i in range(R)[1:-1] for m in range(T)) == 1,
                                'end_origin')
    visit_once = Mod.addConstrs((gp.quicksum(xijm[i,j,m] for i in range(R)[:-1] for m in range(T)) == yi[j]
                                for j in range(R)[1:]), 'visit_once')
    flow_constr = Mod.addConstrs(((gp.quicksum(xijm[i,k,m] for i in range(R)[:-1] for m in range(T))
                                - gp.quicksum(xijm[k,j,m] for j in range(R)[1:] for m in range(T)) == 0)
                                for k in range(R)[1:-1]), 'flow_constr')
    tijmdef = Mod.addConstrs((ti[i] == gp.quicksum(tijm[i,j,m] for j in range(R)[1:] for m in range(T)) 
                    for i in range(R)[:-1]),'tijm_def')
    exp_constr = Mod.addConstrs((ti[j] >= tijm[i,j,m] + (dist[(i,j,m)] + s_i[j])*xijm[i,j,m] 
                                    for i in range(R)[:-1] for j in range(R)[1:] for m in range(T))
                                , 'exp_constr')
    t_min = Mod.addConstrs(((tijm[i,j,m] >= w[m]*xijm[i,j,m]) for j in range(R) for m in range(T)
                                    for i in range(R)),'t_min')
    t_max = Mod.addConstrs((tijm[i,j,m] <= w[m+1]*xijm[i,j,m] for j in range(R) for m in range(T) 
                                    for i in range(R)),'t_max')
    rise_constr = Mod.addConstrs((ti[i] >= e_i[i]*yi[i] for i in range(R)),'rise_constr')
    set_constr = Mod.addConstrs((ti[i] <= l_i[i]*yi[i] for i in range(R)),'set_constr')
    slew_tracker = Mod.addVar(vtype=GRB.CONTINUOUS,name='slew_tracker')
    slew_tracker_def = Mod.addConstr((slew_tracker == gp.quicksum(dist[(i,j,m)]*xijm[i,j,m] for i in range(R)[1:-1] 
                                    for j in range(R)[1:-1] for m in range(T))),'slew_def')

    if args.time_limit is not None:
        Mod.params.TimeLimit = int(args.time_limit)
    Mod.params.MIPGap = gap

    Mod.setObjective(slew_tracker,GRB.MINIMIZE)
    Mod.update()
    logfile = '{}_{}t_{}s_{}'.format(sim_type,num_targs,num_slots,Mod.fingerprint)
    Mod.Params.LogFile = logfile
    Mod.optimize()

    if Mod.SolCount > 0:
        tau_slews = []
        for i in range(R)[1:-1]:
            for j in range(R)[1:-1]:
                for m in range(T):
                    if np.round(xijm[i,j,m].X,0) ==1:
                        tau_slews.append(dist[i,j,m]*60)

        real_slews = []
        for i in range(R)[1:-1]:
            for j in range(R)[1:-1]:
                for m in range(T):
                    if np.round(tijm[i,j,m].X,1) != 0:
                        minutes = np.round(tijm[i,j,m].X,1)
                        time_of_slew = start + TimeDelta(minutes*60,format='sec')
                        altaz1 = get_alt_az(cls_frame,ind_to_id[i],time_of_slew,keck)
                        altaz2 = get_alt_az(cls_frame,ind_to_id[j],time_of_slew,keck)
                        real_slews.append(np.round(max(np.abs(altaz1[0]-altaz2[0]),np.abs(altaz1[1]-altaz2[1])),2))

        total_tau_slew = sum(tau_slews)/60
        med_tau_slew = np.median(tau_slews)
        total_real_slew = sum(real_slews)/60
        med_real_slew = np.median(real_slews)
        runtime = Mod.Runtime

        outputframe.loc[len(outputframe.index)] = [sim_type,num_targs,num_slots,runtime,total_real_slew,med_real_slew,total_tau_slew,med_tau_slew,logfile]
        outputframe.to_csv('statistics_{}_{}t_{}s.csv'.format(sim_type,num_targs,num_slots))

        success = True
        
    Mod.dispose()