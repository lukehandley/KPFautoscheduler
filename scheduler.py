import argparse
import optimize


parser = argparse.ArgumentParser(description='Generate Schedules for Upcoming Night')

parser.add_argument('-o','--observers_sheet', help='Path to observers csv',
                default='example_data/HIRES_observers - PI-requests2022B.csv')
parser.add_argument('-t','--twilight_times', help='Path to twilight times csv',
                default='example_data/twilight_times.csv')
parser.add_argument('-a','--allocated_nights', help='Path to night allocations',
                default='example_data/hires_schedule_2022B.csv')
parser.add_argument('-m','--marked_scripts', help='Path to script directory',
                default='example_data/MarkedScripts')
parser.add_argument('-d','--current_day', help='Date to be scheduled YYYY-MM-DD. Must be in allocated_nights',
                default='2022-10-04')
parser.add_argument('-g','--gurobi_output',action='store_true', help='Activate Gurobi console outputs',
                default=False)
parser.add_argument('-p','--equalize',action='store_true',help='EXPERIMENTAL: Alter Objective to prioritize program equity',
                default=False)

args = parser.parse_args()


#Call the scheduling function
optimize.semester_schedule(args.observers_sheet,args.twilight_times,args.allocated_nights,
                        args.marked_scripts,args.current_day,args.gurobi_output,args.equalize)


