## Roadmap & TODOs
* DONE SIMULINAC with SIGMA matrices?
* DONE pretty lattice print? 1/2
* DONE low energy lattice? 5Mev - 50Mev?
* DONE nonlinear long. dynamics? non-linear (see Shislo, Holmes paper!), not for lin. SIMULINAC!
* DONE SIMULINAC tracking with nonlin. cavities? better with pyORBIT!
* DONE elements with  Position member
* DONE pyORBIT vergleich (MAC,Windows,Linux), my candidate for PIC!
* DONE neuen python driver tracker.py (python for simple traking) schreiben
* DONE common python lattice parser for simulinac & pyOrbit (with MAD syntax?) * maybe later!
* DONE Ripken-Schmidt DRIFT & QUADS DESY_95_063 anstelle von T3D drifts? wichtig, denn thin-lens matrices make most of the length of the machine!!
* DONE refactor: position
* DONE refactor map invocation
* DONE Tanke/Valero gap-model!

* DROPPED test with electrons why?
* DROPPED ELEGANT vergleich? have it on Windows - is it running?
* DROPPED DYNAC vergleich (Linux)
* DROPPED Monte Carlo für verluste (mit DYNAC?)
* DROPPED Cavitychain: phase advance from (i) to (f); kann man vergessen siehe helper5.py

* POSTPONED ELEGANT vergleich & custom cavity?
* POSTPONED generic Display module
* POSTPONED Display-Function class
* POSTPONED acceptance ellipse? not needed?

* VIELLEICHT doublet & triplet lattices? pre-tests with TRACE3D (Windows & wine on Mac)
* VIELLEICHT with magnetic- & alignment-errors!

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
* #TODO: mehr 'namedtupel' benutzen
* #TODO: [Covariance Ellipse](https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html)
