import sys
import os

files = os.listdir()
# print(files)
tables = []
for f in files:
    if 'dynacEzTab' in f: 
        tables.append(f)
# print(tables)

plot_command = 'gnuplot -e \"'
plot_command +=   "set xlabel '[m]'"
plot_command += "; set ylabel '[V/m]'"
plot_command += "; plot "
for f in tables:
    plot_command += f" '{f}' skip 1 with lines,"
plot_command +='; pause mouse keypress"'
# process gnuploy
os.system(plot_command)


