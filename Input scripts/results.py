#!/usr/bin/env python
import os
import sys
import glob

cwd = os.getcwd()

dG_solvent=0
error_solvent=0
dG_complex=0
error_complex=0
found_solvent_result=False
found_complex_result=False

os.chdir("output")
output_files = glob.glob('*.out')
for output_file in output_files:
	output_data = open(output_file)
	for line in output_data:
		line.strip()
		data=line.split()
		if(len(data) != 0  and data[0] == "Relative"):
			if(data[5] == "solvent"):
				dG_solvent=data[7]
				error_solvent=data[10]
				found_solvent_result=True
			elif(data[5] == "complex"):
				dG_complex=data[7]
				error_complex=data[10]
				found_complex_result=True
	output_data.close()

if(not found_solvent_result):
	if(not found_complex_result):
		print("FATAL ERRROR: Found neither solvent nor complex results. Exiting")
		exit(0)
	else:
		print("FATAL ERRROR: Found no solvent results. Exiting")
		exit(0)
elif(not found_complex_result):
	print("FATAL ERRROR: Found no complex results. Exiting")
	exit(0)


os.chdir(cwd)

dG_solvent=float(dG_solvent)
dG_complex=float(dG_complex)
error_solvent=float(error_solvent)
error_complex=float(error_complex)
ddG=dG_complex-dG_solvent
results_filename="ddG.res"
results_data=open(results_filename,'w')
results_data.write("ddG total: {0:.2f} kcal/mol - \n".format(ddG))
results_data.write("dG complex: {0:.2f} kcal/mol +/- {1:.2f} \n".format(dG_complex, error_complex))
results_data.write("dG solvent: {0:.2f} kcal/mol +/- {1:.2f} \n".format(-dG_solvent, error_solvent))
results_data.close()
