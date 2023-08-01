import argparse
import optimize


parser = argparse.ArgumentParser(description='Generate Schedules for Upcoming Night')

parser.add_argument('-o','--observers_sheet', help='Path to observers csv',
                default='2023_data/HIRES_observers - HIRES_Requests2023B.csv')
parser.add_argument('-t','--twilight_times', help='Path to twilight times csv',
                default='2023_data/twilight_times.csv')
parser.add_argument('-a','--allocated_nights', help='Path to night allocations',
                default='2023_data/hires_schedule_2023B.csv')
parser.add_argument('-m','--marked_scripts', help='Path to script directory',
                default='2023_data/MarkedScripts')
parser.add_argument('-d','--schedule_dates',action='append',help='Date(s) to be scheduled as YYYY-MM-DD. Must be in allocated_nights')
parser.add_argument('-g','--gurobi_output',action='store_true', help='Activate Gurobi console outputs',
                default=False)
parser.add_argument('-e','--equalize',action='store_true',help='EXPERIMENTAL: Alter Objective to prioritize program equity',
                default=False)
parser.add_argument('-p','--plot_results',action='store_true', help='Turn on plotting for semester/upcoming quarter nights',
                default=True)

args = parser.parse_args()


#Call the scheduling function
optimize.semester_schedule(args.observers_sheet,args.twilight_times,args.allocated_nights,
                        args.marked_scripts,args.schedule_dates,args.gurobi_output,args.equalize,args.plot_results)
