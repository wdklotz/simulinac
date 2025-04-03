class RFC_TODO(Node):   
    """ Rf cavity as product DKD*RFG*DKD """
    pass
class RFG_OLD(Node):
    """ KAPUTT-BROKEN-HORS SERVICE  KAPUTT-BROKEN-HORS SERVICE  KAPUTT-BROKEN-HORS SERVICE  KAPUTT-BROKEN-HORS SERVICE """
    """  RF-gap of zero length with different kick gap-models """
    # def __init__(self, label, EzPeak, phisoll, gap, cavlen,freq, SFdata=0, particle=UTIL.Proton(UTIL.PARAMS['injection_energy']), position=(0.,0.,0.), aperture=None, dWf=UTIL.FLAGS['dWf'], mapping='t3d'):
    def __init__(self, label, **kwargs):
        super().__init__()
        self.label        = label
        self.length       = 0. # 0. because it's a kick
        self.viseo        = 0.25
        self.accelerating = True
        self.dWf          = UTIL.FLAGS['dWf']                 # dWf=1 with acceleration =0 else

        self.EzPeak    = kwargs.get('EzPeak',0)                 # [MV/m] peak gap field
        self.phisoll   = kwargs.get('phisoll',None)             # [radians] soll phase
        self.gap       = kwargs.get('gap',None)                 # [m] rf-gap
        self.cavlen    = kwargs.get('cavlen',None)              # [m] cavity length
        self.freq      = kwargs.get('freq',None)                # [Hz]  RF frequenz
        self.SFdata    = kwargs.get('SFdata',None)              # SuperFish data
        self.particle  = kwargs.get('particle',None)
        self.position  = kwargs.get('position',None)
        self.aperture  = kwargs.get('aperture',None)
        self.mapping   = kwargs.get('mapping',None)             # map model

        self.EzPeak    = self.EzPeak*self.dWf         # [MV/m] peak gap field
        self.omega     = twopi*self.freq
        self.lamb      = UTIL.PARAMS['clight']/self.freq
        self.ttf       = None
        self.E0L       = None
        self.qE0LT     = None
        self.deltaW    = None
        self.particlef = None
        self.rf_gap    = None
    def dispatch_model_matrix(self):
        """ dispatching to different gap models """
        if self.mapping   == 't3d' :   #NOTE: t3d mapping is matrix multiplication
            self.matrix    = self.gap_object.T3D_matrix()
            # DEBUG_OFF('matrix',self.matrix)
            pass
        elif self.mapping == 'simple':
            self.matrix    = self.gap_object.simple_matrix()  #NOTE: simple mapping is matrix multiplication
            # DEBUG_OFF('matrix',self.matrix)
            pass
        elif self.mapping == 'oxal':
            (self.matrix,self.ttf,self.deltaW)  = self.gap_object.OXAL_matrix(self.particle.tkin) #NOTE: OXAL mapping is matrix multiplication
            DEBUG_OFF(f'{self.gap_object},matrix',self.matrix)
            pass
        elif self.mapping == 'base':
            self.matrix    = self.T3D_matrix(self.ttf,self.particle,self.particlef,self.E0L,self.phisoll,self.lamb,self.deltaW,self.length)
            self.particlef = None
            self.map = self.base_map_1
        elif self.mapping == 'ttf':
            self.matrix    = self.T3D_matrix(self.ttf,self.particle,self.particlef,self.E0L,self.phisoll,self.lamb,self.deltaW,self.length)
            self.particlef = None
            # self.map  =>  # TTF_G has its own mapping method
        elif self.mapping == 'dyn':
            self.matrix    = self.T3D_matrix(self.ttf,self.particle,self.particlef,self.E0L,self.phisoll,self.lamb,self.deltaW,self.length)
            self.particlef = None
            # self.map  =>  # DYN_G has its own mapping method
        # TODO other mappings not tested
        else:
            raise(UserWarning(f"INFO: RFG is a kick-model and does not work with {self.mapping} mapping! Use one of [t3d,simple,base,ttf,oxal]."))
            sys.exit()
    @property
    def gap_object(self):
        return self.rf_gap
    @gap_object.setter
    def gap_object(self,rf_gap):
        self.rf_gap = rf_gap
    def T3D_matrix(self):
        def ttf(lamb, gap, beta):
            """ Panofsky transit-time-factor (see Lapostolle CERN-97-09 pp.65, T.Wangler pp.39) """
            x = gap/(beta*lamb)
            res =NP.sinc(x)
            return res
        """ RF gap-matrix nach Trace3D pp.17 (LA-UR-97-886) """
        m              = NP.eye(MDIM,MDIM)
        self.E0L       = self.EzPeak*self.gap
        self.ttf       = ttf(self.lamb,self.gap,self.particle.beta)
        self.qE0LT     = self.E0L*self.ttf
        self.deltaW    = self.E0L*self.ttf*M.cos(self.phisoll)
        self.particlef = UTIL.Proton(self.particle.tkin+self.deltaW)
        Wavg    = self.particle.tkin+self.deltaW/2.   # average tkin
        pavg    = UTIL.Proton(Wavg)
        bavg    = pavg.beta
        gavg    = pavg.gamma
        m0c2    = pavg.e0
        kz      = twopi*self.E0L*self.ttf*M.sin(self.phisoll)/(m0c2*bavg*bavg*self.lamb)
        ky      = kx = -0.5*kz/(gavg*gavg)
        bgi     = self.particle.gamma_beta
        bgf     = self.particlef.gamma_beta
        bgi2bgf = bgi/bgf
        m       = NP.eye(MDIM,MDIM)
        m[XPKOO, XKOO] = kx/bgf;    m[XPKOO, XPKOO] = bgi2bgf
        m[YPKOO, YKOO] = ky/bgf;    m[YPKOO, YPKOO] = bgi2bgf
        m[ZPKOO, ZKOO] = kz/bgf;    m[ZPKOO, ZPKOO] = bgi2bgf   # koppelt z,z'
        m[EKOO, DEKOO] = self.deltaW
        m[SKOO, DSKOO]  = 0.
        return m
    def simple_matrix(self):
        """ Simplified Matrix Model. (A.Shislo 4.1) """
        def ttf(lamb, gap, beta):
            """ Panofsky transit-time-factor (see Lapostolle CERN-97-09 pp.65, T.Wangler pp.39) """
            x = gap/(beta*lamb)
            res =NP.sinc(x)
            return res
        self.E0L       = self.EzPeak*self.gap
        self.ttf       = ttf(self.lamb,self.gap,self.particle.beta)
        self.qE0LT     = self.E0L*self.ttf
        self.deltaW    = self.E0L*self.ttf*M.cos(self.phisoll)

        # xi   = i_track[XKOO]       # 0
        # xpi  = i_track[XPKOO]      # 1
        # yi   = i_track[YKOO]       # 2
        # ypi  = i_track[YPKOO]      # 3
        # zi   = i_track[ZKOO]       # 4 z
        # zpi  = i_track[ZPKOO]      # 5 Dp/p
        # T    = i_track[EKOO]       # 6 kinetic energy REF
        # dT   = i_track[DEKOO]      # 7 delta energy REF
        # S    = i_track[SKOO]       # 8 position REF
        # dS   = i_track[DSKOO]      # 9 delta position REF

        C         = UTIL.PARAMS['clight']
        particle  = self.particle
        m0c2      = particle.e0
        WIN       = particle.tkin
        betai     = particle.beta
        gammai    = particle.gamma
        gbi       = particle.gamma_beta

        deltaW    = self.deltaW
        lamb      = self.lamb
        qE0LT     = self.qE0LT
        phisoll   = self.phisoll

        WOUT       = WIN + deltaW
        particlef  = copy(particle)(tkin = WOUT)
        betaf      = particlef.beta
        gammaf     = particlef.gamma
        gbf        = particlef.gamma_beta

        condTdP = 1./(m0c2*betaf**2*gammaf)  # conversion DW --> Dp/p
        DDW = self.omega/C/betai*self.qE0LT*M.sin(self.phisoll) # A.Shishlo/J.Holmes (4.1.8)

        m = NP.eye(MDIM,MDIM)
        m[1,0] = m[3,2] = -(UTIL.pi*qE0LT/(m0c2*lamb*gbi*gbi*gbf))*M.sin(phisoll)   # x',x = y',y   (4.1.11)
        # m[1,1] = m[3,3] = m[5,5] = gbi/gbf # x',x' = y',y'  m[5,5] like T3D
        m[1,1] = m[3,3] = gbi/gbf            # x',x' = y',y': m[5,5]=1 like in Shishlo's paper
        m[4,4] = betaf/betai                 # z,z
        m[5,4] = DDW*condTdP                 # z,z'
        m[6,7] = deltaW                      # T,dT
        return m
    def base_map_0(self, i_track):
        """alte map version bis 02.02.2022"""
        # def DEBUG_TRACK(inout,track):
        #     print('{} {} {}'.format('base_map',inout,track))
        # function body ================= function body ================= function body ================= 
        """ Mapping (i) to (f) in BASE RF-Gap Model. (A.Shislo 4.2) """
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z
        zp       = i_track[ZPKOO]      # [5] dp/p
        T        = i_track[EKOO]       # [6] kinetic energy SOLL
        S        = i_track[SKOO]       # [8] position SOLL

        particle = self.particle
        m0c2     = particle.e0
        betai    = particle.beta
        gammai   = particle.gamma
        gbi      = particle.gamma_beta
        tki      = particle.tkin
        freq     = self.freq
        lamb     = self.lamb
        phisoll  = self.phisoll
        qE0LT    = self.qE0LT
        deltaW   = self.deltaW
        
        # if 0: 
        #     DEBUG_ON()
        #     DEBUG_TRACK('tr_i',i_track)
        max_r  = 0.05              # max radial excursion [m]
        r      = M.sqrt(x**2+y**2)  # radial coordinate
        if r > max_r:
            raise UTIL.OutOfRadialBoundEx(S)
        Kr     = (twopi*r)/(lamb*gbi)
        i0     = UTIL.I0(Kr)                               # bessel function I0
        i1     = UTIL.I1(Kr)                               # bessel function I1
        # if 0: print('Kr=',Kr,'r=',r,'gbi=',gbi,'i0=',i0,'i1=',i1)
        # SOLL
        WIN       = tki                               # energy (i)
        DELTAW    = deltaW                       # energy kick
        WOUT      = WIN + DELTAW                      # energy (f) (4.1.6) A.Shishlo/J.Holmes
        # PARTICLE
        converter = UTIL.WConverter(WIN,freq)
        # phin      = -z * twopi/(betai*lamb) + phis     # phase (i)  alte methode
        phin      = converter.zToDphi(z) + phisoll          # phase (i)
        deltaW    = qE0LT*i0*M.cos(phin)                   # energy kick
        # win     = (zp * (gammai+1.)/gammai +1.) * WIN  # energy (i) dp/p --> dT alte methode
        win       =  converter.Dp2pToDW(zp) + WIN         # energy (i) dp/p --> dT
        wout      = win + deltaW                         # energy (f)   (4.2.3) A.Shishlo/J.Holmes
        dw        = wout - WOUT                          # d(deltaW)

        # DEBUG_OFF('base_map: (deltaW,qE0LT,i0,phis) = ({},{},{},{})'.format(deltaW,qE0LT,i0,phis))

        particlef = UTIL.Proton(WOUT)       # !!!IMPORTANT!!! SOLL particle (f)
        betaf     = particlef.beta
        gammaf    = particlef.gamma
        gbf       = particlef.gamma_beta

        # converter = WConverter(WOUT,frq)
        z         = betaf/betai*z                     # z (f) (4.2.5) A.Shishlo/J.Holmes
        # zpf     = gammaf/(gammaf+1.) * dw/WOUT      # dW --> dp/p (f)  alte methode
        zpf       = converter.DWToDp2p(dw)            # dW --> dp/p (f)
        # if 0: print('z ',z,'zpf ',zpf)

        commonf = qE0LT/(m0c2*gbi*gbf)*i1             # common factor
        if r > 0.:
            xp  = gbi/gbf*xp - x/r*commonf*M.sin(phin)  # Formel 4.2.6 A.Shishlo/J.Holmes
            yp  = gbi/gbf*yp - y/r*commonf*M.sin(phin)  # should be phi-middle
        elif r == 0.:
            xp  = gbi/gbf*xp
            yp  = gbi/gbf*yp

        f_track = NP.array([x, xp, y, yp, z, zpf, T+deltaW, 1., S, 1.])

        # for DEBUGGING
        # if 0: DEBUG_TRACK('tr_f',f_track)
        if 0:
            UTIL.arrprnt([x*1.e3 for x in i_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')
            UTIL.arrprnt([x*1.e3 for x in f_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')
        # the parent reads these attributes below
        self.particlef = particlef
        return f_track
    def base_map_1(self, i_track):
        """Neue map Version ab 03.02.2022 ist ein Remake um Korrecktheit der Rechnung zu testen. 
           Produziert dasselbe Verhalten wie base_map_0 """
        # def DEBUG_TRACK(inout,track):
        #     print('{} {} {}'.format('base_map',inout,track))
        # function body ================= function body ================= function body ================= 
        """ Mapping (i) to (O) in BASE RF-Gap Model. (A.Shislo 4.2) """
        x        = i_track[XKOO]       # [0]
        xp       = i_track[XPKOO]      # [1]
        y        = i_track[YKOO]       # [2]
        yp       = i_track[YPKOO]      # [3]
        z        = i_track[ZKOO]       # [4] z
        zp       = i_track[ZPKOO]      # [5] dp/p
        T        = i_track[EKOO]       # [6] kinetic energy ref Teilchen
        S        = i_track[SKOO]       # [8] position gap

        particleRi = self.particle   # ref Teilchen (I)
        m0c2       = particleRi.e0
        betai      = particleRi.beta
        gammai     = particleRi.gamma
        gbi        = particleRi.gamma_beta
        wRi        = particleRi.tkin
        freq       = self.freq
        lamb       = self.lamb
        phisoll    = self.phisoll
        deg_phisoll= M.degrees(phisoll)
        qE0LT      = self.qE0LT
        deltaW     = self.deltaW
        
        # if 0: 
        #     DEBUG_ON()
        #     DEBUG_TRACK('tr_i',i_track)

        max_r  = 0.05              # max radial excursion [m]
        r      = M.sqrt(x**2+y**2)   # radial coordinate
        if r > max_r:
            raise UTIL.OutOfRadialBoundEx(S)
        Kr     = (twopi*r)/(lamb*gbi)
        i0     = UTIL.I0(Kr)            # bessel function I0
        i1     = UTIL.I1(Kr)            # bessel function I1

        # if 0: print('Kr=',Kr,'r=',r,'gbi=',gbi,'i0=',i0,'i1=',i1)

        # ref Teilchen
        wRo = wRi + deltaW                           # ref Teilchen energy (O)
 
        # Teilchen
        converter   = UTIL.WConverter(wRi,freq)
        deg_converter = M.degrees(converter.zToDphi(z)) 
        phiin       = converter.zToDphi(z) + phisoll 
        deg_phiin   = M.degrees(phiin)        # Teilchen phase (I)
        wo_wi       = qE0LT*i0*M.cos(phiin)                 # energy kick (Shislo 4.2.3)
        wi          =  converter.Dp2pToDW(zp) + wRi        # Teilchen energy (I) dp/p --> dT
        wo          = wi + wo_wi                          # Teilchen energy (O)   
        dw          = wo - wRo                            # Differenz der energy kicks von Teilchen und ref Teilchen (entspricht delta**2)

        # DEBUG_OFF('base_map: (deltaW,qE0LT,i0,phis) = ({},{},{},{})'.format(deltaW,qE0LT,i0,phis))

        particleRo = UTIL.Proton(wRo)
        betao      = particleRo.beta
        gammao     = particleRo.gamma
        gbo        = particleRo.gamma_beta

        zo         = betao/betai*z                     # z (O) (4.2.5) A.Shishlo/J.Holmes
        zpo        = converter.DWToDp2p(dw)            # dW --> dp/p (O)

        # if 0: print('z ',z,'zpf ',zpf)

        factor = qE0LT/(m0c2*gbi*gbo)*i1               # common factor
        if r > 0.:
            xp  = gbi/gbo*xp - x/r*factor*M.sin(phiin)   # Formel 4.2.6 A.Shishlo/J.Holmes
            yp  = gbi/gbo*yp - y/r*factor*M.sin(phiin)
        elif r == 0.:
            xp  = gbi/gbo*xp
            yp  = gbi/gbo*yp

        f_track = NP.array([x, xp, y, yp, zo, zpo, T+deltaW, 1., S, 1.])

        # for DEBUGGING
        # if 0: DEBUG_TRACK('tr_f',f_track)
        # if 0:
        #     arrprnt([x*1.e3 for x in i_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')
        #     arrprnt([x*1.e3 for x in f_track[:-4]], fmt = '{:7.4g},', txt = 'base_map:i_track(x[mm],xp,y[mm],yp,z[mm],zp)=')

        # """ the parent reads these attributes below """
        self.particlef = particleRo
        return f_track
    def adjust_energy(self, tkin):
        if self.gap_object == self:
            adjusted = RFG(self.label,EzPeak=self.EzPeak,phisoll=self.phisoll,gap=self.gap,cavlen=self.cavlen,freq=self.freq,SFdata=self.SFdata,particle=UTIL.Proton(tkin),position=self.position,aperture=self.aperture,dWf=self.dWf,mapping=self.mapping)
            adjusted.gap_object = adjusted
            adjusted.dispatch_model_matrix()
        elif self.mapping == 'oxal':
            self.particle = UTIL.Proton(tkin)
            self.gap_object.OXAL_matrix(tkin)
            adjusted = self
            # adjusted.dispatch_model_matrix()
        return adjusted
    def waccept(self):
        """ 
        Calculate longitudinal acceptance, i.e. phase space ellipse parameters: T.Wangler (6.47-48) pp.185
        (w/w0)**2 + (Dphi/Dphi0)**2 = 1
        emitw = w0*Dphi0 = ellipse_area/pi
        """
        rf_gap    = self.gap_object      # this RF gap to use: can be self or others like OXAL_G or TTF_G

        Ez0       = rf_gap.EzPeak
        ttf       = rf_gap.ttf
        phisoll   = rf_gap.phisoll         # [rad]
        lamb      = rf_gap.lamb            # [m]
        freq      = rf_gap.freq            # [Hz]
        particle  = rf_gap.particle

        E0T       = Ez0*ttf              # [MV/m]
        m0c2      = particle.e0          # [MeV]
        gb        = particle.gamma_beta
        beta      = particle.beta
        gamma     = particle.gamma
        tkin      = particle.tkin
        # DEBUG_OFF("waccept",dict(E0T=E0T,phisoll=degrees(phisoll),lamb=lamb,freq=freq,m0c2=m0c2,gb=gb,beta=beta,gamma=gamma,tkin=tkin))

        # converter for this node
        conv = UTIL.WConverter(tkin,freq)

        try:
            # LARGE amplitude oscillations (T.Wangler pp. 175 6.28). w = Dgamma = DW/m0c2 normalized energy spread """
            # DEBUG_OFF(f'w2phi {(1,m0c2,Ez0,ttf,gamma,beta,lamb,phisoll,phisoll)}')                                                                                                                                                              
            w0large = M.sqrt(w2phi(1,m0c2,Ez0,ttf,gamma,beta,lamb,phisoll,phisoll))
            # DEBUG_OFF(f'w0large {w0large}')                                                                                                                                                              
        except ValueError as ex:
            exception = ex
            w0large = -1
        try:
            # SMALL amplitude oscillations separatrix (T.Wangler pp.185) """
            w0small = M.sqrt(2.*E0T*gb**3*lamb*phisoll**2*M.sin(-phisoll)/(M.pi*m0c2))
            # DEBUG_OFF(f'w0small {w0small}')                                                                                                                                                              
        except ValueError as ex:
            exception = ex
            w0small = -1
        if w0large != -1: 
            wmax = w0large
        elif w0large == -1 and w0small != -1:
            wmax = w0small
        else:
            DEBUG_OFF(f'{exception} reason: ttf={rf_gap.ttf}, E0T={E0T}')
            sys.exit(1)

        # this node Dp/p max on separatrix
        Dp2pmax = conv.wToDp2p(wmax) 

        #  convert T.Wangler units {Dphi,w} to {z,dp/p} units with 1st cavity parameters
        betaw_i,alfaw_i,gammaw,emitw_i = UTIL.PARAMS['twiss_w_i']()
        Dphi0_i = UTIL.PARAMS['Dphi0_i']
        w0_i = (gamma-1.)*UTIL.PARAMS['DT2T_i']
        z0_i,Dp2p0_i,emitz_i,betaz_i = conv.wtoz((Dphi0_i,w0_i,emitw_i,betaw_i))
        alfaz_i = 0.

        # omega sync for this node
        omgl0zuomg = M.sqrt(E0T*lamb*M.sin(-phisoll)/(2*M.pi*m0c2*gamma**3*beta))
        omgl_0     = omgl0zuomg*twopi*freq   # [Hz]

        # longitudinal acceptance check (always done)     #TODO  is this needed?
        # if wmax <= w0_i:
        #     si,sm,sf = self.position
        #     warnings.showwarning(
        #         colors.RED+'out of energy acceptance @ s={:.1f} [m]'.format(si)+colors.ENDC,
        #         UserWarning,'elements.py',
        #         'waccept')

        # phase acceptance (REMARK: phase limits are not dependent on Dp/p aka w)
        phi_2=2.*phisoll
        phi_1=-phisoll

        res =  dict (
                emitw_i         = emitw_i,      # 1st cavity emittance {Dphi,w} units [rad,1]
                z0_i            = z0_i,         # 1st cavity ellipse z-axe crossing (1/2 axis) [m]
                Dp2p0_i         = Dp2p0_i,      # 1st cavity ellipse dp/p-axe crossing (1/2 axis)
                twiss_z_i       = UTIL.Twiss(betaz_i, alfaz_i, emitz_i), # 1st cavity twis parameters
                DWmax           = wmax*m0c2,    # this node max delta-W on separatrix [MeV]
                Dp2pmax         = Dp2pmax,      # this node Dp/p max on separatrix [1]
                phaseacc        = (conv,phi_2,phisoll,phi_1), # this node phase acceptance [rad]
                omgl_0          = omgl_0,       # this node synchrotron oscillation [Hz]
                wmax            = wmax,         # this node w max on separatrix [1] (large amp. oscillations)
                zmax            = conv.DphiToz(-phisoll) # this node z max on separatrix [m] (large amp. oscillations -- Wrangler's approximation (pp.178) is good up to -58deg)
                )
        return res
    def aper_check(self,new_tp,s,**kwargs):
        new_point=new_tp()
        fifo_z = kwargs['fifo_z']
        sfifo_z= kwargs['sfifo_z']
        fifo_xy= kwargs['fifo_xy']
        sfifo_xy=kwargs['sfifo_xy']

        rf_gap = self.gap_object      # this RF gap to use: can be self or others like OXAL_G or TTF_G

        lost=False
        tkin=rf_gap.particle.tkin

        res = self.waccept()
        dummy,phi_2,phisoll,phi_1 = res['phaseacc']   # dummy kept for old track_node
        conv = UTIL.WConverter(tkin,rf_gap.freq)
        Dphi=conv.zToDphi(new_point[Ktp.z])
        phi = phisoll+Dphi

        # longitudinal acceptance
        if not (phi_2 < phi and phi < phi_1):  # Wrangler's approximation (pp.178) is good up to -58deg
            fifo_z.append(f'loss (z) {new_point[Ktp.z]:.3e} at {s:.4e} m')
            sfifo_z.append(s)
            lost = True
        elif abs(new_point[Ktp.zp]) > res['Dp2pmax']:
            fifo_z.append(f'loss (zp) {new_point[Ktp.zp]:.3e} at {s:.4e} m')
            sfifo_z.append(s)
            lost = True
        # transverse apertures
        elif rf_gap.aperture != None and not (abs(new_point[Ktp.x]) < rf_gap.aperture or abs(new_point[Ktp.y]) < rf_gap.aperture):
            fifo_xy.append(f'loss (x|y) ({new_point[Ktp.x]:.3e},{new_point[Ktp.y]:.3e}) at {s:.4e} m')
            sfifo_xy.append(s)
            lost = True
        return lost
