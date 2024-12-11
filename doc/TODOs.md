# Roadmap & TODOs
# done items
* SIMULINAC with SIGMA matrices?
* pretty lattice print? 1/2
* low energy lattice? 5Mev - 50Mev?
* nonlinear long. dynamics? non-linear (see Shislo, Holmes paper!), not for lin. SIMULINAC!
* SIMULINAC tracking with nonlin. cavities? better with pyORBIT!
* elements with  Position member
* pyORBIT vergleich (MAC,Windows,Linux), my candidate for PIC!
* neuen python driver tracker.py (python for simple traking) schreiben
* common python lattice parser for simulinac & pyOrbit (with MAD syntax?) * maybe later!
* Ripken-Schmidt DRIFT & QUADS DESY_95_063 anstelle von T3D drifts? wichtig, denn thin-lens matrices make most of the length of the machine!!
* refactor: position
* refactor map invocation
* Tanke/Valero gap-model
# dropped items
*  test with electrons why?
*  ELEGANT vergleich? have it on Windows - is it running?
*  DYNAC vergleich (Linux)
*  Monte Carlo für verluste (mit DYNAC?)
*  Cavitychain: phase advance from (i) to (f); kann man vergessen siehe helper5.py
# postponed items
* ELEGANT vergleich & custom cavity?
* generic Display module
* Display-Function class
* acceptance ellipse? not needed?
# maybe items
* doublet & triplet lattices? pre-tests with TRACE3D (Windows & wine on Mac)
* with magnetic- & alignment-errors!
# todo items
* end energy variabel 80Mev-200MeV!
* injection at ~15MeV!!!
* rf-phase rippel!!!
* 180 deg bend!!
* TEAPOT spacing of thin-qudrupoles!!!
* check dwf with new element classes
* Auswertung Bunch tracking
* BMAD vergleich!
* DYNAC with Python from ESS
* track whole bunch (all particles) through element then go to next element instead of reverse. this will avoid recalculation of map for each particle. i think i already do that.
* mehr 'namedtupel' benutzen
* [Covariance Ellipse](https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html)
