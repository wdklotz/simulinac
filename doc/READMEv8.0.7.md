# SIMULINAC
Is a python simulator/tracker for linac structures.

I worked on a user friendly version of my code. Here it is!!!

## Getting started:
* You have to decompress the zip-archive `v8.0.7.zip` then change directory:
* $**cd v8.0.7** then decompress the zip-archive `SIMULINAC.zip`:
* then change directory **cd SIMULINAC** and run the demo:
* $**python simu.py yml/simuIN.yml**
* You must have python 3 (I tested with 3.4).
* The program depends on the modules: _matplolib, numpy, scipy, PyYaml and tkinter._
* The demo input file is **yml/simuIN.yml**. Copy it and modify at your will.
* You can track particle bunches:
* **$python tracker.py yml/trackIN.yml**
* The demo input file is **yml/trackIN.yml**. Copy it and modify at your will.
* You can get help from me by mail to wdklotz@gmail.com
## Structure of the input file:
* A **LATTICE** is an array of N **LINES**:
  <pre>- [125,*line]</pre>
* A **LINE** is an array of **CELLS**:
  <pre>- [*cell]</pre>
*  A **CELL** is an array of **SEGMENTS**:
   <pre>- [*sqf1,*srfc,*sqd,*srfc,*sqf2]</pre>
* A **SEGMENT** is an array of **NODES**: 
  <pre>- SQD:   &sqd  [*d1,*qd,*d1]  # QD + LR drifts</pre>
* A **NODE** is an **ELEMENT**: 
  <pre>- D1:   &d1   D3</pre>
* An **ELEMENT** is an array of key:value pairs:     
  <pre>- D3:                         # ID: &alias
       - type:     D            # D class
       - length:   0.03         # [m]
       - sec:      *HE          # is part of section</pre>
* **PARAMETERS** and **FLAGS** are similar to ELEMENTS, i.e. array of (key:value) pairs.

### WARNING: _You have to follow strictly the YAML syntax in the input files!_
* The top level blocks titled: **FLAGS:**, **SECTIONS:**, **PARAMETERS:**, **ELEMENTS:**, **NODES:**, **SEGMENTS:**, **CELL:**, **LINE:** and **LATTICE:** are mandatory. The lattice-parser will not parse correctly when you replace or rename them.

-- have fun --
