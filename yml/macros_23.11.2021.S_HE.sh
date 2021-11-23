#============== Injektion energy
T=17.5                     # kinetic energy in [Mev]

#============== Injektion energy spread
DT2T=1.0e-3             # delta-T/T kinetic

#============== X,Y emittances
# EMITX=1.e-6             # x emittance in [m*rad]
# EMITY=1.e-6             # y emittance in [m*rad]
EMITX=4.e-6             # x emittance in [m*rad]
EMITY=4.e-6             # y emittance in [m*rad]
BETAX=3.617             # twiss beta x in [m]  (T,B')=(25,23)
BETAY=0.709             # twiss beta x in [m]  (T,B')=(25,23)

#============== RF
PHISY=-30.              # synchronous phase in [deg]
FREQ=816.e6             # common rf-frequency [Hz] (T>=25)

#============== Quad gradients
BGRAD=23.000            # quad gradient [T/m] (T=25)

#============== cavity mapping
MAP=t3d
# MAP=simple
# MAP=oxal
# MAP=base
# MAP=ttf
# MAP=dyn

ARGS="-D _TKIN=$T"
ARGS="$ARGS -D _DT2T=$DT2T"
ARGS="$ARGS -D _EMITX=$EMITX"
ARGS="$ARGS -D _EMITY=$EMITY"
ARGS="$ARGS -D _BETAX=$BETAX"
ARGS="$ARGS -D _BETAY=$BETAY"
ARGS="$ARGS -D _PHISY=$PHISY"
ARGS="$ARGS -D _FREQ=$FREQ"
ARGS="$ARGS -D _BGRAD=$BGRAD"
ARGS="$ARGS -D _MAPPING=$MAP"

# invoke m4
m4 $ARGS $1 > $2
