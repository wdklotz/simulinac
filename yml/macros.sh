TKIN=5.                 # kinetic energy in [Mev]
TKIN=70.                # kinetic energy in [Mev]
BGRAD=43.150

MAPPING=t3d          # Trace 3D linear map model
# MAPPING=simple       # Shishlo/Holmes linear map model
# MAPPING=base         # Shishlo/Holmes base map model
# MAPPING=ttf          # Shishlo/Holmes three point TTF RF gap-model
# MAPPING=dyn          # Tanke/Valero DYNAC RF gap-model

EMITW=1.0e-5            # w emittance in [rad]
EMITW=1.0e-5            # w emittance in [rad]
EMITX=1.0e-6            # x emittance in [m*rad]
EMITY=1.0e-6            # y emittance in [m*rad]

BETAX=3.3               # best w/o rf
# BETAX=0.65              # good with tkin=5
BETAY=2.5               # best w/o rf
# BETAY=0.49              # good with tkin=5

ALFAX=0.
ALFAY=0.

PHISY=-30.              # synchronous phase in [deg]
NCELL=1
NCELL=6
# NCELL=41
NCELL=166
# NCELL=240

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
ARGS="$ARGS -D _MAPPING=$MAPPING"

# invoke m4
m4 $ARGS $1
