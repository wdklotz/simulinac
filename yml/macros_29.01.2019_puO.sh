TKIN=25.                # kinetic energy in [Mev]

EMITX=1.e-6             # x emittance in [m*rad]
EMITY=1.e-6             # y emittance in [m*rad]
EMITW=1.e-5             # w emittance in [rad]

BETAX=3.3               # twiss beta x in [m]
BETAY=0.55              # twiss beta x in [m]

PHISY=-25.              # synchronous phase in [deg]

NCELL=100

MAP=base

ARGS="-D _TKIN=$TKIN"
ARGS="$ARGS -D _EMITW=$EMITW"
ARGS="$ARGS -D _EMITX=$EMITX"
ARGS="$ARGS -D _EMITY=$EMITY"
ARGS="$ARGS -D _BETAX=$BETAX"
ARGS="$ARGS -D _BETAY=$BETAY"
ARGS="$ARGS -D _PHISY=$PHISY"
ARGS="$ARGS -D _NCELL=$NCELL"
ARGS="$ARGS -D _MAPPING=$MAP"

# invoke m4
m4 $ARGS $1