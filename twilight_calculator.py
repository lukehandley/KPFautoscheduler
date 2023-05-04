import astropy as apy
import astroplan as apl
import pandas as pd
from astropy.time import Time

#Initialize the astroplan observer object
keck = apl.Observer.at_site('W. M. Keck Observatory')

#Create a range of dates to calculate (some before and some after semester in this case)
twilight_frame = pd.date_range('2022-06-01','2024-06-01').to_frame()
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
twilight_frame.to_csv('2023_data/twilight_times.csv')