TKIN=50.                # kinetic energy in [Mev]
##EMITW=12.6e-6         # w emittance in [rad]
##EMITW=12.6e-5         # w emittance in [rad]
EMITW=9.6e-5            # w emittance in [rad]
##EMITX=1.e-6           # x emittance in [m*rad]
##EMITY=1.e-6           # y emittance in [m*rad]
EMITX=2.e-6             # x emittance in [m*rad]
EMITY=2.e-6             # y emittance in [m*rad]
BETAX=2.285             # twiss beta x in [m]
BETAY=0.328             # twiss beta x in [m]
PHISY=-20.              # synchronous phase in [deg]

ARGS="-D _TKIN=$TKIN"
ARGS="$ARGS -D _EMITW=$EMITW"
ARGS="$ARGS -D _EMITX=$EMITX"
ARGS="$ARGS -D _EMITY=$EMITY"
ARGS="$ARGS -D _BETAX=$BETAX"
ARGS="$ARGS -D _BETAY=$BETAY"
ARGS="$ARGS -D _PHISY=$PHISY"

# invoke m4
m4 $ARGS $1
