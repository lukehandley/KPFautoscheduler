import argparse

parser = argparse.ArgumentParser(description='Generate Schedules for Upcoming Night')

parser.add_argument('-o','--observers_sheet', help='Path to observers csv')
parser.add_argument('-t','--twilight_times', help='Path to twilight times csv',default='HIRES_2023B_data/twilight_times.csv')
parser.add_argument('-n','--allocated_nights', help='Path to night allocations')
parser.add_argument('-m','--marked_scripts', help='Path to script directory')
parser.add_argument('-i','--instrument', help='Instrument to schedule, currently HIRES or KPF')
parser.add_argument('-d','--schedule_dates',action='append',help='Date(s) to be scheduled as YYYY-MM-DD. Must be in allocated_nights')
parser.add_argument('-g','--gurobi_output',action='store_true', help='Activate Gurobi console outputs',
                default=False)
parser.add_argument('-e','--equalize',action='store_true',help='EXPERIMENTAL: Alter Objective to prioritize program equity',
                default=False)
parser.add_argument('-p','--plot_results',action='store_true', help='Turn on plotting for semester/upcoming quarter nights',
                default=False)
parser.add_argument('-f','--folder', help='Folder in which to save generated scripts and plots')
parser.add_argument('-l','--time_limit', help='Maximum time spent in each optimization instance in seconds',type=int,
                default=300)
parser.add_argument('-a','--accessibility_constant',help='Fraction of night targets must be visible',
                default=6)

args = parser.parse_args()

if args.folder is None:
    if args.instrument == 'KPF':
        args.folder = 'KPF_2023B_data'
    if args.instrument == 'HIRES':
        args.folder = 'HIRES_2023B_data'
if args.observers_sheet is None:
    if args.instrument == 'KPF':
        args.observers_sheet = 'KPF_2023B_data/KPF_target_Requests - 2023B_requests.csv'
    if args.instrument == 'HIRES':
        args.observers_sheet = 'HIRES_2023B_data/HIRES_KPF_observers - HIRES_Requests2023B.csv'
if args.allocated_nights is None:
    if args.instrument == 'KPF':
        args.allocated_nights = 'KPF_2023B_data/kpf_schedule_2023B.csv'
    if args.instrument == 'HIRES':
        args.allocated_nights = 'HIRES_2023B_data/hires_schedule_2023B.csv'
if args.marked_scripts is None:
    if args.instrument == 'KPF':
        args.marked_scripts = 'KPF_2023B_data/MarkedScripts'
    if args.instrument == 'HIRES':
        args.marked_scripts = 'HIRES_2023B_data/MarkedScripts'


import optimize

#Call the scheduling function
optimize.semester_schedule(args.instrument,args.observers_sheet,args.twilight_times,args.allocated_nights,
                        args.marked_scripts,args.schedule_dates,args.gurobi_output,args.equalize,args.plot_results,
                        args.folder,args.time_limit,args.accessibility_constant)


