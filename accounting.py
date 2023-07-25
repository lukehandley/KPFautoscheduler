import pandas as pd

def completion_logs(all_targets_frame,observed_dict,starlists,current_day):
    targets = []
    requested = []
    observed_past = []
    forecasted = []
    programs = []
    flattened_starlist = [star for sublist in starlists for star in sublist]
    for index,row in all_targets_frame.iterrows():
        targets.append(row['Starname'])
        requested.append(row['N_obs(full_semester)'])
        observed_past.append(len(observed_dict[row['request_number']]))
        forecasted.append(flattened_starlist.count(row['request_number']))
        programs.append(row['Program code'])

    data = {'Name':targets, 'Program':programs, 'Requested':requested, 'Observed':observed_past, 'Forecasted':forecasted}
    df = pd.DataFrame.from_dict(data)
    df.to_csv('{}_completion_logs.csv'.format(current_day))
