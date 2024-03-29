#============== Injektion energy
TIN=20.
TIN=5.

#============== cavity mapping
# MAP=t3d
# MAP=simple
# MAP=oxal
MAP=base
# MAP=ttf
# MAP=dyn

#============== cavity frquency
FREQ=408.e6
# FREQ=816.e6

#============== Injektion energy spread
# DT2T=1.0e-3             # delta-T/T kinetic

#============== X,Y emittances
# EMITX=1.e-6             # x emittance in [m*rad]
# EMITY=1.e-6             # y emittance in [m*rad]
# EMITX=4.e-6             # x emittance in [m*rad]
# EMITY=4.e-6             # y emittance in [m*rad]

#============== twiss parameters
# FODO matched
# BETAX=0.692             # twiss beta x in [m]
# BETAY=1.600             # twiss beta y in [m]
# ALFAX=0.000             # twiss alfa x []
# ALFAY=0.000             # twiss alfa y []

#============== RF
# PHISY=-30.              # synchronous phase in [deg]
# FREQ=816.e6             # common rf-frequency [Hz] (T>=25)

#============== Quad gradients
BGRAD=25.0

      ARGS="-D _TIN=$TIN"
ARGS="$ARGS -D _MAPPING=$MAP"
ARGS="$ARGS -D _BGRAD=$BGRAD"
ARGS="$ARGS -D _FREQ=$FREQ"
# ARGS="$ARGS -D _DT2T=$DT2T"
# ARGS="$ARGS -D _EMITX=$EMITX"
# ARGS="$ARGS -D _EMITY=$EMITY"
# ARGS="$ARGS -D _BETAX=$BETAX"
# ARGS="$ARGS -D _BETAY=$BETAY"
# ARGS="$ARGS -D _ALFAX=$ALFAX"
# ARGS="$ARGS -D _ALFAY=$ALFAY"
# ARGS="$ARGS -D _PHISY=$PHISY"

# invoke m4
m4 $ARGS $1 > $2
