"""
Version v7.0.6alpha
*) added loss plots to traker
*) unified the coordinate handling: two vectors: twiss-vector (Ktw) and track point (Ktp)
*) verified correctness of twiss-envelopes and sigma-envelopes

Version v7.0.5
*) code cleaning
*) Marker-actions reduced from list to single variable
*) cleaned scatter plotting
*) n_sigma PARAMETER activated
*) useaper FLAG activated
*) aperture checks with n_sigma in simu and tracker
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