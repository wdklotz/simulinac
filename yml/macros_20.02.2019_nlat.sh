T=50.                   # kinetic energy in [Mev]
T=80.                   # kinetic energy in [Mev]
DT2T=6.0e-3             # delta-T/T kinetic
DT2T=1.0e-3             # delta-T/T kinetic

EMITX=1.e-6             # x emittance in [m*rad]
EMITY=1.e-6             # y emittance in [m*rad]

BETAX=3.65              # twiss beta x in [m]
BETAY=0.66              # twiss beta x in [m]

PHISY=-30.              # synchronous phase in [deg]

BGRAD=43.150            # quad gradient [T/m]

NL=54                   # nboff lines a.k.a. cells

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
ARGS="$ARGS -D _BGRAD=$BGRAD"
ARGS="$ARGS -D _NL=$NL"
ARGS="$ARGS -D _MAPPING=$MAP"

# invoke m4
m4 $ARGS $1
