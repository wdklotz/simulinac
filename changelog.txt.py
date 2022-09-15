"""
Version v7.1.3a4
*) input_files for Windows now static now
*) Ez0.py refactored. Neue Intervalteilung
*) improved speed of calculation in OXAL.py
*) added TmStamp class in setutil.py to check program flow 
*) removed depricated time.clock() in tracker.py
*) removed emitw_i as initial parameter. Vorgabe is now relative kinetic energy spread DT/T @ entrance
"""

"""
Version v7.1.3a4
*) updated documents in doc-dir
*) longitudinal phase space now initialized by DT/T kinetic energy spread
*) new SUMMARY print-out
*) helper: bmParams.py to calculate longitudinal parameters for comparisons
"""
"""
Version v7.1.3
*) corection of oxal-map for wrong formulas in S&H paper

Version v7.1.2
*) added the openXAL linear gap-model from Shishlo & Holmes
*) the implemetation with sympy does not produce correct results
*) various checks and bug fixes
*) code cleaning and commenting

Version v7.1.1
*) Functions class for easier plotting
*) TTFG: map completed with inclusion of S(k) ans S'(k) coefficients
*) code cleaning: class _Node attributes are no longer overriden by inheritance

Version v7.1.0a1
*) Lattice is now double linked list
*) All models now implemented as DKD nodes
*) dyn mapping working
*) EzAvg-bug fixed
*) aperture-bug fixed

Version v7.0.7a2
*) totally new DYNAC rf-gap model replaces all ealier which were incorrect
*) sign error for z correction/step in _DYN_G fixed 
*) mapping is global now

Version v7.0.6a2
*) added loss plots to traker
*) unified the coordinate handling: two vectors: twiss-vector (Ktw) and track point (Ktp)
*) solved the mystery between twiss-envelopes and sigma-envelopes differences. Replaced the old
twiss-envelope-calculation by a new, simpler one. Both show now similar and seamingly
correct results.

Version v7.0.5
*) code cleaning
*) Marker-actions reduced from list to single variable
*) cleaned scatter plotting
*) nbsigma PARAMETER activated
*) useaper FLAG activated
*) aperture checks with nbsigma in simu and tracker
*) input of aperture improved

Version v7.0.4
*) Particle is not a dictionary anymore instead it has track as attribute now.
*) fixed bugs when FLAG['dWf'] is set to no acceleration
*) launching without m4 preprocessing like 'python (simu.py | tracker.py) <input>' possible now
*) ERROR corrected: RFG objects did not receive the correct wavelength. 
*) PoincareAction class introduced (in module marker_actions.py). Generation of poincare-sections
   now possible at any marker position. Frames go to directory ./frames.

Version v7.0.3
*) Bug fixed: some elements lost their dictionary values.
*) Bug fixed: the rf-gaps did not have the T3D linear matrix They were treated as unit matrices and 
showed emittance growths with energy when using the non-linear maps. I never understood
this behaviour when comparing with the t3d mapping. This artefact has been eliminated.
*) simu.py can be started with 3 file arguments: 
    'python simu.py <tmplate-file>  <macro-definition-file>  <input_file>'
*) renamed worktmpl.yml to tmpl.yml per default
*) Track class is now a simple list of Tpoint objects
*) lattice version is now printed in simu and tracker
*) code clean up

Version v7.0.2
*) elements D,QF,QD,RFG,RCAV,GAP have all object attibutes in their private dictionary to
    make the construction of lattice-converter scripts easier
*) preliminary User's Guide as tex-document
*) new parameter-print with KVout=True. Includes now traceX and traceY for stability estimation.

Version v7.0.1
removed bug with show flag in tracker.py
"""