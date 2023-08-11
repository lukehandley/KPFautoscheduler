import numpy as np
import pandas as pd
import logging
import os

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def discretize(seconds):
    return np.round(seconds/60,1)

def format_observers_sheet(observers_sheet,instrument):
    pd.set_option('display.precision', 30)
    if instrument == 'HIRES':
        #Turn the observers sheet into a dataframe that is simplified and useful for caluculations. This assumes the 
        #sheet is left downloaded directly
        all_targets_frame = pd.read_csv(observers_sheet)
        all_targets_frame = all_targets_frame.drop([0,1,2,3])
        all_targets_frame['Starname'] = all_targets_frame['Starname'].str.slice(0, 16)


        ############Compute coordinates in decimal RA hours and decimal DEC deg############
        all_targets_frame['ra'] = all_targets_frame['RAH'] + (1/60)*all_targets_frame['RAM'] + (1/3600)*all_targets_frame['RAS']

        #Reading in dec can be troublesome
        all_targets_frame['dec'] = '' 
        for index,row in all_targets_frame.iterrows():
            #Remove Spaces
            row['Starname'] = row['Starname'].replace(" ", "")
            if row['DECD'] < 0:
                dec = (np.abs((int(row['DECD']))) + (1/60)*row['DECM'] + (1/3600)*row['DECS'])
                all_targets_frame.at[index,'dec'] = -np.abs(dec)
            else:
                dec = (int(row['DECD']) + (1/60)*row['DECM'] + (1/3600)*row['DECS'])
                all_targets_frame.at[index,'dec'] = np.abs(dec)
        
        all_targets_frame['N_obs(full_semester)'] = pd.to_numeric(all_targets_frame['N_obs(full_semester)'])
        all_targets_frame['N_obs'] = pd.to_numeric(all_targets_frame['N_obs'])
        all_targets_frame['T_exp(sec)'] = pd.to_numeric(all_targets_frame['T_exp(sec)'])
        all_targets_frame['T_max(sec)'] = pd.to_numeric(all_targets_frame['T_max(sec)'])
        all_targets_frame['Cadence'] = pd.to_numeric(all_targets_frame['Cadence'])

        all_targets_frame['discretized_duration'] = all_targets_frame['T_exp(sec)'].apply(discretize)
        #For triple shots, we assume the total exposure time is triple the listed nominal for our calculations
        for index,row in all_targets_frame.iterrows():
            if row['N_obs'] > 1:
                all_targets_frame.at[index,'discretized_duration'] = row['discretized_duration'] * row['N_obs']

        #Remove the duplicated targets from our list that are calibration shots not relevant to community cadence
        dup = all_targets_frame[all_targets_frame.duplicated(subset='Starname',keep=False)]
        all_targets_frame = all_targets_frame.drop(dup[dup['N_obs(full_semester)'] == 1].index.tolist())
        all_targets_frame = all_targets_frame[(all_targets_frame['Done'] != 'done') & (all_targets_frame['Done'] != 'Y')
                                               & (all_targets_frame['Done'] != 'y')]
        all_targets_frame = all_targets_frame[all_targets_frame['N_obs(full_semester)'] > 0]

        #Finally assign a unique identifier to each request to be used in lookup tables
        all_targets_frame.reset_index(inplace=True,drop=True)
        all_targets_frame['request_number'] = all_targets_frame.index.tolist()

        return all_targets_frame

    if instrument == 'KPF':
        #Turn the observers sheet into a dataframe that is simplified and useful for caluculations. This assumes the 
        #sheet is left downloaded directly
        all_targets_frame = pd.read_csv(observers_sheet,dtype={'Gaia_DR3_id': object})
        all_targets_frame = all_targets_frame.drop([0,1,2])
        all_targets_frame['Starname'] = all_targets_frame['Starname'].str.slice(0, 16)


        ############Compute coordinates in decimal RA hours and decimal DEC deg############
        all_targets_frame['ra'] = all_targets_frame['RAH'] + (1/60)*all_targets_frame['RAM'] + (1/3600)*all_targets_frame['RAS']

        #Reading in dec can be troublesome
        all_targets_frame['dec'] = '' 
        for index,row in all_targets_frame.iterrows():
            #Remove Spaces
            row['Starname'] = row['Starname'].replace(" ", "")
            if row['DECD'] < 0:
                dec = (np.abs((int(row['DECD']))) + (1/60)*row['DECM'] + (1/3600)*row['DECS'])
                all_targets_frame.at[index,'dec'] = -np.abs(dec)
            else:
                dec = (int(row['DECD']) + (1/60)*row['DECM'] + (1/3600)*row['DECS'])
                all_targets_frame.at[index,'dec'] = np.abs(dec)
        
        all_targets_frame['N_obs(full_semester)'] = pd.to_numeric(all_targets_frame['N_obs(full_semester)'])
        all_targets_frame['Nobs per'] = pd.to_numeric(all_targets_frame['Nobs per'])
        all_targets_frame['T_exp(sec)'] = pd.to_numeric(all_targets_frame['T_exp(sec)'])
        all_targets_frame['T_max(sec)'] = pd.to_numeric(all_targets_frame['T_max(sec)'])
        all_targets_frame['Cadence'] = pd.to_numeric(all_targets_frame['Cadence'])
        all_targets_frame['Jmag'] = pd.to_numeric(all_targets_frame['Jmag'])
        all_targets_frame['Gmag'] = pd.to_numeric(all_targets_frame['Gmag'])
        all_targets_frame['Nvisits'] = all_targets_frame['Nvisits'].fillna(1)

        all_targets_frame['discretized_duration'] = all_targets_frame['T_exp(sec)'].apply(discretize)
        #For triple shots or multi visit targets, we assume the total exposure time is multiplied by the number
        for index,row in all_targets_frame.iterrows():
            if row['Nobs per'] > 1:
                all_targets_frame.at[index,'discretized_duration'] = row['discretized_duration'] * row['Nobs per']
                #Iron out if Nobs semester should be divided by per obs
                #all_targets_frame.at[index,'N_obs(full_semester)'] = row['N_obs(full_semester)']/row['Nobs per']
            try:
                if row['Nvisits'] > 1:
                    all_targets_frame.at[index,'discretized_duration'] = row['discretized_duration'] * row['Nvisits']
            except:
                if int(row['Nvisits'][-1]) > 1:
                    all_targets_frame.at[index,'discretized_duration'] = row['discretized_duration'] * int(row['Nvisits'][-1])

        #Remove the RM targets corresponding to the RM nights
        all_targets_frame = all_targets_frame[(all_targets_frame['RM or'] != 'Y') & (all_targets_frame['RM or'] != 'y')]

        #Remove done targets
        all_targets_frame = all_targets_frame[(all_targets_frame['Done'] != 'done') & (all_targets_frame['Done'] != 'Y')
                                               & (all_targets_frame['Done'] != 'y')]
        all_targets_frame = all_targets_frame[all_targets_frame['N_obs(full_semester)'] > 0]

        #Finally assign a unique identifier to each request to be used in lookup tables
        all_targets_frame.reset_index(inplace=True,drop=True)
        all_targets_frame['request_number'] = all_targets_frame.index.tolist()
        
        return all_targets_frame

def format_allocated_nights(allocated_nights,instrument):
    obs_plan = pd.read_csv(allocated_nights,dtype={'start': float, 'stop': float})
    if instrument == 'HIRES':
        obs_plan = obs_plan[obs_plan['Date'] != '2023-11-22']
        obs_plan = obs_plan[obs_plan['Date'] != '2022-12-02']
        obs_plan.reset_index(inplace=True,drop=True)

        return obs_plan

    if instrument == 'KPF':
        obs_plan = obs_plan[obs_plan['Queue'] == 'y']
        obs_plan.reset_index(inplace=True,drop=True)

        return obs_plan

def write_starlist(instrument,frame,requests,condition,current_day,outputdir):
    script_dir = os.path.join(outputdir,'{}_scripts'.format(instrument))
    if not os.path.isdir(script_dir):
        os.mkdir(script_dir)
    script_file = os.path.join(script_dir,'{}_{}_{}.txt'.format(instrument,current_day,condition))
    logger.info('Writing starlist to ' + script_file)
    if instrument == 'HIRES':
        #Retrieve columns from observers sheet that are relevant to the starlist
        columns = ['Starname','RAH','RAM','RAS','DECD','DECM','DECS','epoch','vmag=','Vmag',
            'T_exp(sec)','T_max(sec)','Exp_meter_counts(k)','decker','N_obs','Iodine',
            'Priority','Program code','Telescope Propsoal Code','Comment']

        formatted_frame = frame.loc[requests][columns]

        lines = []
        for index,row in formatted_frame.iterrows():
            #Just a bunch of string formatting. This prints standard starlists as ordered by the salesman optimization

            namestring = ' '*(16-len(row['Starname'][:16])) + row['Starname'][:16]

            rastring = ('0'*(2-len(str(int(row['RAH'])))) + str(int(row['RAH'])) + ' '
                            + '0'*(2-len(str(int(row['RAM'])))) + str(int(row['RAM'])) + ' '
                                + '0'*(4-len(str(np.round(row['RAS'],1)))) + str(np.round(row['RAS'],1)))

            starter = '+'
            if row['DECD'] < 0:
                starter = '-'
                decstring = (starter + '0'*(2-len(str(abs(int(row['DECD']))))) + str(abs(int(row['DECD']))) + ' '
                                + '0'*(2-len(str(int(row['DECM'])))) + str(int(row['DECM'])) + ' '
                                    + '0'*(2-len(str(int(row['DECS'])))) + str(int(row['DECS'])))
            else:
                decstring = (starter + '0'*(2-len(str(abs(int(row['DECD']))))) + str(abs(int(row['DECD']))) + ' '
                                + '0'*(2-len(str(int(row['DECM'])))) + str(int(row['DECM'])) + ' '
                                    + '0'*(2-len(str(int(row['DECS'])))) + str(int(row['DECS'])))

            magstring = (row['vmag='] + str(np.round(float(row['Vmag']),1)) + ' '*(4-len(str(np.round(row['Vmag'],1)))))

            exposurestring = (' '*(4-len(str(int(row['T_exp(sec)'])))) + str(int(row['T_exp(sec)'])) + '/' 
                                + str(int(row['T_max(sec)'])) + ' '*(4-len(str(int(row['T_max(sec)'])))))

            countstring = (' '*(3-len(str(int(row['Exp_meter_counts(k)']))))
                                + str(int(row['Exp_meter_counts(k)'])) + 'k')

            line = (namestring + ' ' + rastring + ' ' + decstring + ' ' + str(int(row['epoch'])) + ' '
                        + magstring + ' ' + exposurestring + ' ' + countstring + ' ' + row['decker'] + ' ' + 
                        str(int(row['N_obs'])) + 'x' + ' ' + ' '*(3-len(row['Iodine'])) + row['Iodine']
                        + ' ' + row['Priority'] + ' ' + 
                        row['Program code'] + ' '*(3-len(row['Program code'])) + ' ' + row['Telescope Propsoal Code'])

            if not pd.isnull(row['Comment']):
                line += (' ' + str(row['Comment']))
            
            lines.append(line)

            #Formatted starlists appear as text files in the directory
            with open(script_file, 'w') as f:
                f.write('\n'.join(lines))

    if instrument == 'KPF':
        columns = ['Starname','RAH','RAM','RAS','DECD','DECM','DECS','epoch','jmag=','Jmag',
            'T_exp(sec)','T_max(sec)','Nobs per','Format1','Nvisits','Format2','Simulcal',
            'gmag=','Gmag','Teff', 'Teff_val','Gaia_DR3','Gaia_DR3_id','Program code','Format3',
            'priority','Telescope','Comment']
        
        formatted_frame = frame.loc[requests][columns]

        lines = []
        for index,row in formatted_frame.iterrows():
            #Just a bunch of string formatting. This prints standard starlists as ordered by the salesman optimization

            namestring = ' '*(16-len(row['Starname'][:16])) + row['Starname'][:16]

            rastring = ('0'*(2-len(str(int(row['RAH'])))) + str(int(row['RAH'])) + ' '
                            + '0'*(2-len(str(int(row['RAM'])))) + str(int(row['RAM'])) + ' '
                                + '0'*(4-len(str(np.round(row['RAS'],1)))) + str(np.round(row['RAS'],1))) 

            starter = '+'
            if row['DECD'] < 0:
                starter = '-'
                decstring = (starter + '0'*(2-len(str(abs(int(row['DECD']))))) + str(abs(int(row['DECD']))) + ' '
                                + '0'*(2-len(str(int(row['DECM'])))) + str(int(row['DECM'])) + ' '
                                    + '0'*(4-len(str(row['DECS']))) + str(int(row['DECS'])))
            else:
                decstring = (starter + '0'*(2-len(str(abs(int(row['DECD']))))) + str(abs(int(row['DECD']))) + ' '
                                + '0'*(2-len(str(int(row['DECM'])))) + str(int(row['DECM'])) + ' '
                                    + '0'*(4-len(str(row['DECS']))) + str(int(row['DECS'])))

            jmagstring = (row['jmag='] + str(np.round(float(row['Jmag']),1)) + ' '*(4-len(str(np.round(row['Jmag'],1)))))

            exposurestring = (' '*(4-len(str(int(row['T_exp(sec)'])))) + str(int(row['T_exp(sec)'])) + '/' 
                                + str(int(row['T_max(sec)'])) + ' '*(4-len(str(int(row['T_max(sec)'])))))
            
            try:
                ofstring = ('1of' + str(int(row['Nvisits'])))
            except:
                ofstring = row['Nvisits']
                
            scstring = 'sc=' + row['Simulcal']
            
            numstring = (str(row['Nobs per']) + row['Format1'])
                
            gmagstring = (row['gmag='] + str(np.round(float(row['Gmag']),1)) + ' '*(4-len(str(np.round(row['Gmag'],1)))))
            
            teffstr = row['Teff'] + str(int(row['Teff_val'])) + ' '*(4-len(str(int(row['Teff_val']))))
            
            gaiastring = 'Gaia DR3 ' + str(int(row['Gaia_DR3_id'])) + ' '*(19-len(str(int(row['Gaia_DR3_id']))))
            
            programstring = row['Program code']
            
            priostring = row['Format3'] + str(int(row['priority'])) + ' ' + row['Telescope']

            line = (namestring + ' ' + rastring + ' ' + decstring + ' ' + str(int(row['epoch'])) + ' '
                        + jmagstring + ' ' + exposurestring + ' ' + ofstring + ' ' + scstring +  ' '
                        + numstring + ' '+ gmagstring + ' ' + teffstr + ' ' + gaiastring + ' ' 
                                + programstring + ' ' + priostring)
            
            if not pd.isnull(row['Comment']):
                line += (' ' + str(row['Comment']))

            lines.append(line)

            #Formatted starlists appear as text files in the directory
            with open(script_file, 'w') as f:
                f.write('\n'.join(lines))
