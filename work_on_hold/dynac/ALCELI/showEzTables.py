import sys
import os

files = os.listdir()
# print(files)
tables = []
for f in files:
    if 'dynacEzTab' in f: 
        tables.append(f)
# print(tables)

plot_command = 'gnuplot -persist -e \"'
plot_command +=   "set xlabel '[m]'"
plot_command += "; set ylabel '[V/m]'"
plot_command += "; plot "
for cnt,f in enumerate(tables):
    plot_command += f" '{f}' skip 1 title '{f}' noenhanced with lines,"
    # if cnt == 0: plot_command += f" '{f}' skip 1 title '{f}' noenhanced with lines,"
plot_command +=';"'
# process gnuploy
os.system(plot_command)
