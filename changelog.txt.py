Version v7.0.3
*) Bug fixed: some elements lost their dictionary values.
*) Bug fixed: the rf-gaps did not have the T3D linear matrix They were treated as unit matrices and 
showed emittance growths with energy when using the non-linear maps. I never understood
this behaviour when comparing with the t3d mapping. This artefact has been eliminated.
*) simu.py and tracker.py can now be called with 3 file arguments.
*) renamed worktmpl.yml to tmpl.yml per default
*) Track class is now a simple list of Tpoint objects
*) lattice version is now printed in simu and tracker
*) code clean up

Version v7.0.2
*) elements D,QF,QD,RFG,RCAV,GAP have all object attibutes in their private dictionary to
    make the construction of lattice-converter scripts easier
*) simu.py can be started as: 
    'python simu.py <tmplate-file>  <macro-definition-file> <input_file>'
*) preliminary User's Guide as tex-document
*) new parameter-print with KVout=True. Includes now traceX and traceY for stability estimation.

Version v7.0.1
removed bug with show flag in tracker.py