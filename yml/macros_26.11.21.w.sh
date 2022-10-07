#============== Injektion energy
# T=17.5                     # kinetic energy in [Mev]
# T=6.

#============== cavity mapping
#MAP=t3d
# MAP=simple
# MAP=oxal
MAP=base
# MAP=ttf
# MAP=dyn

#============== Injektion energy spread
# DT2T=1.0e-3             # delta-T/T kinetic

#============== X,Y emittances
EMITX=1.e-6             # x emittance in [m*rad]
EMITY=1.e-6             # y emittance in [m*rad]
# EMITX=4.e-6             # x emittance in [m*rad]
# EMITY=4.e-6             # y emittance in [m*rad]


#============== twiss parameters
# FODO matched
BETAX=0.692             # twiss beta x in [m]
BETAY=1.600             # twiss beta y in [m]
ALFAX=0.000             # twiss alfa x []
ALFAY=0.000             # twiss alfa y []
# 12 cav matched
# BETAX=9.729             # twiss beta x in [m]
# BETAY=2.685             # twiss beta y in [m]
# ALFAX=0.184             # twiss alfa x []
# ALFAY=0.149             # twiss alfa y []

# BETAX=2.205
# BETAY=2.0
# ALFAX=0.057
# ALFAY=0.0

# BETAY=1.137
# BETAX=0.9
# ALFAY=0.1310
# ALFAY=0.05

# NC=300
#============== RF
PHISY=-30.              # synchronous phase in [deg]
FREQ=816.e6             # common rf-frequency [Hz] (T>=25)

#============== Quad gradients
BGRAD=9.000            # quad gradient - FODO matched
BGRAD=16.5

ARGS="-D _TKIN=$T"
ARGS="$ARGS -D _DT2T=$DT2T"
ARGS="$ARGS -D _EMITX=$EMITX"
ARGS="$ARGS -D _EMITY=$EMITY"
ARGS="$ARGS -D _BETAX=$BETAX"
ARGS="$ARGS -D _BETAY=$BETAY"
ARGS="$ARGS -D _ALFAX=$ALFAX"
ARGS="$ARGS -D _ALFAY=$ALFAY"
ARGS="$ARGS -D _PHISY=$PHISY"
ARGS="$ARGS -D _FREQ=$FREQ"
ARGS="$ARGS -D _BGRAD=$BGRAD"
ARGS="$ARGS -D _MAPPING=$MAP"
ARGS="$ARGS -D _NC=$NC"

# invoke m4
m4 $ARGS $1 > $2
