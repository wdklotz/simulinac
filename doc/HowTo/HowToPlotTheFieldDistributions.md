## Plot the field distributions with gnuplot like that
* _change directory_
* **cd .../generated**
* _start gnuplot interactively_
* $ **gnuplot**
* _enter these commands_
* gnuplot> **set xlabel "[m]"**
* gnuplot> **set ylabel "[V/m]"**
* gnuplot> **plot "dynacEzTab_1" skip 1 with lines, "dynacEzTab_481" skip 1 with lines**
* _or run python script as alternative:_
* $ **python showEzTables.py**