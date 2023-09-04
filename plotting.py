from collections import defaultdict
import matplotlib.pyplot as plt
import imageio
import os
import numpy as np
import astropy as apy
import astroplan as apl
from astropy.time import Time
from astropy.time import TimeDelta
import astropy.units as u
import pandas as pd
import math
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

def plot_path_2D(obs_time,az_path,alt_path,names,targ_list,outputdir,current_day):

    fig, axs = plt.subplots(2,sharex=True,sharey=False,figsize = (14,8))
    fig.patch.set_alpha(1)
    axs[0].plot(obs_time,az_path,color = 'indigo')
    axs[0].vlines(obs_time,0,360,linestyle = 'dashed', alpha = .5, color = 'gray')
    axs[0].set_yticks([0,120,240,360],[0,120,240,360])
    ax2 = axs[0].twiny()
    ax2.set_xlim(axs[0].get_xlim())

    topticks = []
    index = 0
    while index < len(obs_time):
        val = (obs_time[index+1]+obs_time[index])/2
        topticks.append(val)
        index+=2

    ax2.set_xticks(topticks)
    ax2.set_xticklabels(names,rotation=45)
    axs[1].plot(obs_time,alt_path,color = 'seagreen')
    axs[1].vlines(obs_time,0,90,linestyle = 'dashed', alpha = .5, color = 'gray')
    #axs[2].plot(obs_time,airmass_path,color = 'firebrick')
    #axs[2].vlines(obs_time,1,3.5,linestyle = 'dashed', alpha = .5, color = 'gray')

    bottomticks = []
    bottomticks.append(obs_time[0])
    for i in range(1,4):
        val = obs_time[0] + i*(obs_time[-1]-obs_time[0])/4
        bottomticks.append(val)
    bottomticks.append(obs_time[-1])

    i = 0
    while i < len(az_path):
        axs[0].fill_betweenx([0,360],obs_time[i],obs_time[i+1],color = 'orange',alpha = .25)
        axs[1].fill_betweenx([0,90],obs_time[i],obs_time[i+1],color = 'orange',alpha = .25)
        #axs[2].fill_betweenx([1,3.5],obs_time[i],obs_time[i+1],color = 'orange',alpha = .25)
        i += 2

    axs[0].set(ylabel='Azimuth Angle (Deg)')
    axs[1].set(ylabel='Elevation Angle (Deg)')
    #axs[2].set(ylabel='Airmass')
    axs[1].set(xlabel='Observation Time (JD)')
    plt.title('Telescope Path {}'.format(current_day))
    plt.savefig(os.path.join(outputdir,'Telescope_Path'))
    plt.close()

def animate_telescope(time_strings,total_azimuth_list,total_zenith_list,tel_az,tel_zen,
                      observed_at_time,plotpath):
           
    theta = np.arange(5.3/180, 146.2/180, 1./180)*np.pi
    total_azimuth_list = np.array(total_azimuth_list)
    total_zenith_list = np.array(total_zenith_list)
    tel_ims_dir = os.path.join(plotpath,'tel_ims')
    if not os.path.isdir(tel_ims_dir):
        os.mkdir(tel_ims_dir)
    
    filenames = []
    for i in range(len(time_strings)):
        if i % 60 == 0:
            fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
            ax.set_ylim(0,70)
            ax.set_title(time_strings[i])
            ax.set_yticklabels([])
            ax.fill_between(theta,56.7,70,color = 'red',alpha=.7)
            #ax.fill_between(np.arange(0,2,1/180)*np.pi,50,70,color= 'red',alpha=.4)
            ax.set_theta_zero_location('N')
            observed_list = observed_at_time[:i]
            for j in set(observed_list):
                ax.scatter(total_azimuth_list[i][j],total_zenith_list[i][j],color='orange',marker='*')
            for j in set(observed_at_time):
                if j not in set(observed_list):
                    ax.scatter(total_azimuth_list[i][j],total_zenith_list[i][j],color='white',marker='*')
            ax.plot(tel_az[:i],tel_zen[:i],color='orange')
            ax.set_facecolor('black')

            # create file name and append it to a list
            filename = f'{i}.png'
            filenames.append(filename)

            # save frame
            plt.savefig(os.path.join(tel_ims_dir,filename),dpi=100)
            plt.close()

    # build gif
    with imageio.get_writer(os.path.join(plotpath,'Observing_Animation.gif'), mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(os.path.join(tel_ims_dir,filename))
            writer.append_data(image)

    # Remove files
    for filename in set(filenames):
        os.remove(os.path.join(tel_ims_dir,filename))
    try:
        os.remove(tel_ims_dir)
    except:
        logger.debug('Cannot remove redundant tel_ims directory due to file permissions')

def plot_program_cdf(plan,program_dict,targets_observed,Nobs,outputdir,current_day):

    fig,axs = plt.subplots(1,figsize = (7,5))
    fig.patch.set_alpha(1)
    colors = plt.cm.rainbow(np.linspace(0, 1, len(program_dict.keys())))
    i = 0
    for program in program_dict.keys():
        tot_obs = sum(Nobs[r] for r in program_dict[program])
        starting_date = plan.loc[0,'Date']
        x = []
        y = []
        targets = []
        for date in plan.Date.unique():
            indeces = plan[plan['Date'] == date].index.values.tolist()
            for index in indeces:
                for targ in targets_observed[index]:
                    if targ in program_dict[program]:
                        targets.append(targ)
            x.append((Time(date,format='iso') - Time(starting_date,format='iso')).jd)
            y.append((len(targets)/tot_obs)*100)
        plt.plot(x,y,label=program,color=colors[i])
        i += 1
    axs.vlines((Time(current_day,format='iso') - Time(starting_date,format='iso')).jd,0,100,color='black',linestyles='dashed')
    axs.set_xticks([0,50,100,150,185])
    axs.set_xlim(0,225)
    axs.set_ylim(-5,105)
    axs.set_ylabel('Completion %')
    axs.set_xlabel('Days into Semester')
    plt.legend()
    plt.savefig(os.path.join(outputdir,'Program_CDFs.png'),dpi=200)
    plt.close()

def plot_one_cdf(plan,program_dict,targets_observed,Nobs,outputdir,current_day,program):

    fig,axs = plt.subplots(1)
    fig.patch.set_alpha(1)
    tot_obs = sum(Nobs[r] for r in program_dict[program])
    starting_date = plan.loc[0,'Date']
    x = []
    y = []
    targets = []
    for date in plan.Date.unique():
        indeces = plan[plan['Date'] == date].index.values.tolist()
        for index in indeces:
            for targ in targets_observed[index]:
                if targ in program_dict[program]:
                    targets.append(targ)
        x.append((Time(date,format='iso') - Time(starting_date,format='iso')).jd)
        y.append((len(targets)/tot_obs)*100)
    plt.plot(x,y,label=program)
    axs.vlines((Time(current_day,format='iso') - Time(starting_date,format='iso')).jd,0,100,color='black',linestyles='dashed')
    axs.set_xticks([0,50,100,150,185])
    axs.set_ylim(-5,105)
    axs.set_ylabel('Completion %')
    plt.legend()
    plt.savefig(os.path.join(outputdir,'{}_CDF.png'.format(program)),dpi=200)
    plt.close()

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

def plot_program_cadence(instrument,plan,all_targets_frame,twilight_frame,starlists,min_separations,outputdir,
                         current_day,accessibility_constant,cadence_mode=False):

    dates = plan.Date.to_list()
    
    sta = []
    sto = []
    dat = []
    for day in pd.date_range(dates[0],dates[-1]).strftime(date_format='%Y-%m-%d').tolist():
        for i in range(4):
            sta.append(i*.25)
            sto.append((i+1)*.25)
            dat.append(day)
    all_qn_plan_dict = {'Date' : dat,'start' : sta, 'stop' : sto}
    all_qn_plan = pd.DataFrame.from_dict(all_qn_plan_dict)

    #Call the interval calculator function
    all_qn_plan = calculate_intervals(all_qn_plan,twilight_frame)
    
    start = Time(dates[0],format='iso',scale='utc')
    stop = start + TimeDelta(1,format='jd')
    step = TimeDelta(60,format='sec')
    t = np.arange(start.jd,stop.jd,step.jd)
    t = Time(t,format='jd')
    min_az = 5.3
    max_az = 146.2
    min_alt = 33.3
    else_min_alt = 25
    #min_alt = 40
    #else_min_alt = 40
    max_alt = 85
    keck = apl.Observer.at_site('W. M. Keck Observatory')
    cadence_folder = os.path.join(outputdir,'Cadence_Plots')
    if not os.path.isdir(cadence_folder):
        os.mkdir(cadence_folder)

    for program in all_targets_frame['Program code'].unique().tolist():
        fig,axs = plt.subplots(1,figsize=(12,12))
        fig.patch.set_alpha(1)
        axs.set_xlim(-5,260)
        height = 0
        if cadence_mode:
            program_frame = all_targets_frame[(all_targets_frame['Program code'] == program) & (all_targets_frame['N_obs(full_semester)'] > 1)]
        if not cadence_mode:
            program_frame = all_targets_frame[(all_targets_frame['Program code'] == program)]
        if len(program_frame) > 0:
            y_positions = []
            y_labels = []
            nobs_list = []
            nvisits_list = []
            target_list = []
            request_list = []
            targ_names = []
            for request_id in program_frame['request_number'].tolist():
                request_list.append(request_id)
                name = program_frame.loc[request_id,'Starname']
                ra = program_frame.loc[request_id,'ra']
                dec = program_frame.loc[request_id,'dec']
                coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
                target = apl.FixedTarget(name=name, coord=coords)
                target_list.append(target)
                targ_names.append(name)
                nobs_list.append(int(program_frame.loc[request_id,'N_obs(full_semester)']))
                if instrument == 'KPF':
                    nvisits_list.append(int(program_frame.loc[request_id,'Nvisits']))

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
            
            b_dict = defaultdict(list)
            c_dict = defaultdict(list)
            for qn in [0,1,2,3]:
                for day in pd.date_range(dates[0],dates[-1]).strftime(date_format='%Y-%m-%d').tolist():
                    day_start = Time(day,format='iso',scale='utc')
                    day_diff = (day_start-start).value
                    offset = math.floor(day_diff/15)
                    shift = day_diff * 4 - offset
                    night_frame = all_qn_plan[all_qn_plan['Date'] == day].reset_index()                        
                    qn_starts = Time(night_frame.loc[qn,'qn_start'],format='jd')
                    qn_stops = Time(night_frame.loc[qn,'qn_stop'],format='jd')
                    b_dict[qn].append(int(np.round((qn_starts - day_start).jd * 24*60,0) + shift))
                    c_dict[qn].append(int(np.round((qn_stops - day_start).jd * 24*60,0) + 1 + shift))

            for i in range(len(observability_matrix)):
                height = i*5
                name = targ_names[i]
                nobs = nobs_list[i]
                if instrument == 'KPF':
                    nvisits = nvisits_list[i]
                total_observed = 0
                a = observability_matrix[i]
                targ = request_list[i]
                min_cadence = int(min_separations[targ])
                total_observed = 0
                for qn in [0,1,2,3]:
                    observability = []
                    if qn == 3:
                        y_positions.append(height+.5)
                        y_labels.append(name)
                    for j in range(len(b_dict[qn])):
                        b = b_dict[qn][j]
                        c = c_dict[qn][j]
                        if c > len(a):
                            one = a[b:]
                            remainder = c % len(a)
                            two = a[:remainder]
                            d = np.concatenate((one,two))
                        else:
                            d = a[b:c]
                        if (len(d) - np.bincount(d)[0]) >= len(d)/accessibility_constant:
                            observability.append(1)
                        else:
                            observability.append(0)
            
                    if len(np.bincount(observability)) > 1:
                        condition = True
                        points = []
                        for j in range(len(observability)):
                            if observability[j] == condition:
                                points.append(j)
                                condition = not condition
                        if len(points) % 2 != 0:
                            points.append(len(observability)-1)
                        axs.hlines([height,height+1],0,190,color = 'black',alpha = .2,linewidth=0.5)
                        j = 0
                        while j < (len(points) - 1):
                            axs.fill_between(np.linspace(points[j],points[j+1]),height,height+1,color = 'lime',alpha = .5)
                            j += 2
                        observed = []
                        for j in range(len(starlists)):
                            if targ in starlists[j]:
                                if plan.loc[j,'start'] == qn/4:
                                    observed.append((Time(plan.loc[j,'Date'],format='iso')-start).jd)
                        axs.vlines(observed,height,height+1,color='black')
                        total_observed += len(observed)
                        allocated = []
                        for index,row in plan.iterrows():
                            if row['start'] == qn/4:
                                allocated.append((Time(row['Date'],format='iso')-start).jd)
                        axs.vlines(allocated,height,height+1,color='blue',alpha=.2)
                        if qn == 3:
                            if instrument == 'KPF':
                                plt.text(190,height-.5,str(min_cadence)
                                                        + ' Days Minimum Cadence, '+ str(total_observed*nvisits) + 
                                                        '/{} Observations Achieved ({} visits each tick)'.format(nobs,nvisits),fontsize=4)
                            else:
                                plt.text(190,height-.5,str(min_cadence)
                                                        + ' Days Minimum Cadence, '+ str(total_observed) + 
                                                        '/%s Observations Achieved' % nobs,fontsize=4)
                    else:
                        observed = []
                        axs.hlines([height,height+1],0,190,color = 'black',alpha = .2,linewidth=0.5)
                        for j in range(len(starlists)):
                            if targ in starlists[j]:
                                if plan.loc[j,'start'] == qn/4:
                                    day = Time(plan.loc[j,'Date'],format='iso')
                                    observed.append((day-start).jd)
                        axs.vlines(observed,height,height+1,color='black')
                        total_observed += len(observed)
                        allocated = []
                        for index,row in plan.iterrows():
                            if row['start'] == qn/4:
                                allocated.append((Time(row['Date'],format='iso')-start).jd)
                        axs.vlines(allocated,height,height+1,color='blue',alpha=.2)
                        if qn == 3:
                            if instrument == 'KPF':
                                plt.text(190,height-.5,str(min_cadence)
                                                        + ' Days Minimum Cadence, '+ str(total_observed*nvisits) + 
                                                        '/{} Observations Achieved ({} visits each tick)'.format(nobs,nvisits),fontsize=4)
                            else:
                                plt.text(190,height-.5,str(min_cadence)
                                                        + ' Days Minimum Cadence, '+ str(total_observed) + 
                                                        '/%s Observations Achieved' % nobs,fontsize=4)
                    height += 1
            axs.vlines((Time(current_day,format='iso')-start).jd,0,height,color='black',alpha=0.15)
            axs.set_yticks(y_positions,y_labels,fontsize=4)
            axs.set_xticks([0,85,171])
            axs.set_title('Cadence for Program %s' % program)
            plt.savefig(os.path.join(cadence_folder,'Program_{}.pdf'.format(program)),dpi=200)
            plt.close()

def plot_cadence_night_resolution(plan,all_targets_frame,twilight_frame,starlists,min_separations,outputdir):

    dates = plan.Date.to_list()
    
    sta = []
    sto = []
    dat = []
    
    start = Time(dates[0],format='iso',scale='utc')
    stop = start + TimeDelta(1,format='jd')
    step = TimeDelta(60,format='sec')
    t = np.arange(start.jd,stop.jd,step.jd)
    t = Time(t,format='jd')
    min_az = 5.3
    max_az = 146.2
    min_alt = 33.3
    else_min_alt = 25
    #min_alt = 40
    #else_min_alt = 40
    max_alt = 85
    keck = apl.Observer.at_site('W. M. Keck Observatory')
    cadence_folder = os.path.join(outputdir,'Cadence_Plots')
    if not os.path.isdir(cadence_folder):
        os.mkdir(cadence_folder)

    for program in all_targets_frame['Program code'].unique().tolist():
        fig,axs = plt.subplots(1,figsize=(12,12))
        fig.patch.set_alpha(1)
        axs.set_xlim(-5,195)
        program_frame = all_targets_frame[all_targets_frame['Program code'] == program]
        y_positions = []
        y_labels = []
        height = 0
        for targ in program_frame['request_number'].tolist():
            name = program_frame.loc[targ,'Starname']
            ra = program_frame.loc[targ,'ra']
            dec = program_frame.loc[targ,'dec']
            dur = program_frame.loc[targ,'discretized_duration']
            nobs = program_frame.loc[targ,'N_obs(full_semester)']
            if nobs > 1:
                coords = apy.coordinates.SkyCoord(ra * u.hourangle, dec * u.deg, frame='icrs')
                target = apl.FixedTarget(name=name, coord=coords)
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
                a = np.zeros(len(t),dtype=int)
                a[good] = 1
                min_cadence = min_separations[targ]
                total_observable = 0
                total_observed = 0
                y_positions.append(height+.5)
                y_labels.append(name)
                observability = []
                for day in pd.date_range(dates[0],dates[-1]).strftime(date_format='%Y-%m-%d').tolist():
                    day_start = Time(day,format='iso',scale='utc')
                    day_diff = (day_start-start).value
                    offset = math.floor(day_diff/15)
                    shift = day_diff * 4 - offset
                    qn_starts = Time(twilight_frame.loc[day,'12_evening'],format='jd')
                    qn_stops = Time(twilight_frame.loc[day,'12_morning'],format='jd')
                    b = int(np.round((qn_starts - day_start).jd * 24*60,0) + shift)
                    c = int(np.round((qn_stops - day_start).jd * 24*60,0) + 1 + shift)
                    if c > len(a):
                        one = a[b:]
                        remainder = c % len(a)
                        two = a[:remainder]
                        d = np.concatenate((one,two))
                    else:
                        d = a[b:c]
                    if (len(d) - np.bincount(d)[0]) >= 75:
                        observability.append(1)
                    else:
                        observability.append(0)
                
                if len(np.bincount(observability)) > 1:
                    condition = True
                    points = []
                    for i in range(len(observability)):
                        if observability[i] == condition:
                            points.append(i)
                            condition = not condition
                    if len(points) % 2 != 0:
                        points.append(len(observability)-1)
                    axs.hlines([height,height+1],0,171,color = 'black',alpha = .2,linewidth=0.5)
                    i = 0
                    while i < (len(points) - 1):
                        axs.fill_between(np.linspace(points[i],points[i+1]),height,height+1,color = 'lime',alpha = .5)
                        i += 2
                    observed = []
                    for i in range(len(starlists)):
                        if targ in starlists[i]:
                            observed.append((Time(plan.loc[i,'Date'],format='iso')-start).jd)
                    axs.vlines(observed,height,height+1,color='black')
                    total_observed += len(observed)
                    allocated = []
                    for index,row in plan.iterrows():
                        allocated.append((Time(row['Date'],format='iso')-start).jd)
                    axs.vlines(list(set(allocated)),height,height+1,color='blue',alpha=.1)
                    plt.text(182,height,str(total_observed) + '/%s' % nobs)
                    height += 1
                else:
                    observed = []
                    axs.hlines([height,height+1],0,171,color = 'black',alpha = .2,linewidth=0.5)
                    for i in range(len(starlists)):
                        if targ in starlists[i]:
                            day = Time(plan.loc[i,'Date'],format='iso')
                            observed.append((day-start).jd)
                    axs.vlines(observed,height,height+1,color='black')
                    total_observed += len(observed)
                    for index,row in plan.iterrows():
                        allocated.append((Time(row['Date'],format='iso')-start).jd)
                    axs.vlines(list(set(allocated)),height,height+1,color='blue',alpha=.1)
                    plt.text(182,height,str(total_observed) + '/%s' % nobs)
                    height += 1
                height += 1
        axs.set_yticks(y_positions,y_labels)
        axs.set_xticks([0,50,100,150,185])
        axs.set_xlabel('Days into Semester')
        axs.set_title('Cadence for Program %s' % program)
        plt.savefig(os.path.join(cadence_folder,'Nightly_Program_{}.png'.format(program)),dpi=200)
        plt.close()