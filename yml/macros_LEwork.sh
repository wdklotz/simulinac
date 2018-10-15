TKIN=10.                # kinetic energy in [Mev]
BGRAD=21.575

EMITW=1.0e-5            # w emittance in [rad]
EMITX=1.e-6             # x emittance in [m*rad]
EMITY=1.e-6             # y emittance in [m*rad]

BETAX=2.030             # best w/o rf th=10
BETAY=0.600             # best w/o rf tk=10

ALFAX=-5
ALFAY=0

PHISY=-30.              # synchronous phase in [deg]
NCELL=285
NCELL=1

ARGS="-D _TKIN=$TKIN"
ARGS="$ARGS -D _EMITW=$EMITW"
ARGS="$ARGS -D _EMITX=$EMITX"
ARGS="$ARGS -D _EMITY=$EMITY"
ARGS="$ARGS -D _BETAX=$BETAX"
ARGS="$ARGS -D _BETAY=$BETAY"
ARGS="$ARGS -D _PHISY=$PHISY"
ARGS="$ARGS -D _NCELL=$NCELL"
ARGS="$ARGS -D _ALFAX=$ALFAX"
ARGS="$ARGS -D _ALFAY=$ALFAY"
ARGS="$ARGS -D _BGRAD=$BGRAD"

# invoke m4
m4 $ARGS $1
