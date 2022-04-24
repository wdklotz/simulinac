    wdklotz@xps8500:/mnt/c/Users/wdklotz/SIMULINAC$ docker start -a tracker
    tracker.py v8.0.6a2 on python 3.9.2 on linux
    run Version 20.02.2019_nlat
    macros=macros_20.02.2019_nlat
    template=tmpl_20.02.2019_nlat
    input=yml/trackIN.yml
-----------------------track_bunch with 1750 particles---

          ==================Tracker Log==================
           DT/T-kin................ :    0.006
           Dp2p.................[%] : 0.303944
           T-kin..............[MeV] :       25
           acceptance Dp2p......[%] : 0.559983
           accpetance z........[mm] :  20.7895

           betaw_i............[rad] :  2071.87
           betax_i..............[m] :    3.617
           betay_i..............[m] :    0.709
           betaz_i..........[m/rad] :   1.4423
           emitw_i,wmx........[rad] : 5.29526e-05 0.00029454
           emitx_i..............[m] :    4e-06
           emity_i..............[m] :    4e-06
           emitz_i..............[m] : 1.33243e-05
           lattice version......... : 29.01.19_pyO
           mapping................. : base
           sigma(Dphi,w)i..([rad,]) : 0.331226 0.000159868
           sigma(x,x')i...([m,rad]) : 0.00380368 0.0014148
           sigma(y,y')i...([m,rad]) : 0.00168404 0.00319555
           sigma(z,Dp2p)i....([m,]) : 0.00438379 0.00303944

    (track design) (track bunch)   100% done 0% lost
    TRACKING DONE (live particles 1750, lost particles 0)


|                   | WSL2 (Docker) | W10 (native) | VBox (LxMint+Docker) | W10 (PYZO) | VBox (FreeBSD) |
| ----------------- | ------------- | ------------ | -------------------- | ---------- | -------------- |
|   tot. time [sec] |       622.028 |     656.516  |         751.945      | 573.672    | 730.508        |
| track bunch [sec] |       575.461 |     652.719  |         730.879      | 500.734    | 722.697        |
| relative          | 1             | 1.13         | 1.27                 | 0.87       | 1.26           |