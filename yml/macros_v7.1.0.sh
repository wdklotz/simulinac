TKIN=5.                 # kinetic energy in [Mev]
TKIN=50.                # kinetic energy in [Mev]
TKIN=80.                # kinetic energy in [Mev]
BGRAD=43.150

MAPPING=t3d          # Trace 3D linear map model
# MAPPING=simple       # Shishlo/Holmes linear map model
# MAPPING=oxal           # openXAL linear model
# MAPPING=base         # Shishlo/Holmes base map model
# MAPPING=ttf          # Shishlo/Holmes three point TTF RF gap-model
# MAPPING=dyn          # Tanke/Valero DYNAC RF gap-model

EMITW=1.0e-5            # w emittance in [rad]
# EMITW=1.0e-7            # w emittance in [rad]
EMITX=1.0e-6            # x emittance in [m*rad]
EMITY=1.0e-6            # y emittance in [m*rad]

BETAX=3.3               # tao
BETAY=0.75              # tao
BETAX=3.713    # closed
BETAY=2.501    # closed

ALFAX=0.
ALFAY=0.

PHISY=-30.              # synchronous phase in [deg]
NCELL=5
NCELL=20
# NCELL=23    # dynac limit  80->113 Mev
# NCELL=54
# NCELL=85
# NCELL=40
# NCELL=60
# NCELL=80
# NCELL=160
# NCELL=200
# NCELL=300

####################################################
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
