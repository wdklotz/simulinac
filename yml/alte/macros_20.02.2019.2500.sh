# T=5.                    # kinetic energy in [Mev]
# T=20.                   # kinetic energy in [Mev]
T=25.                   # kinetic energy in [Mev]
# T=50.                   # kinetic energy in [Mev]
# T=80.                   # kinetic energy in [Mev]

DT2T=6.0e-3             # delta-T/T kinetic
# DT2T=1.0e-3             # delta-T/T kinetic

EMITX=1.e-6             # x emittance in [m*rad]
EMITY=1.e-6             # y emittance in [m*rad]
EMITX=4.e-6             # x emittance in [m*rad]
EMITY=4.e-6             # y emittance in [m*rad]

BETAX=3.617             # twiss beta x in [m]  (T,B')=(25,23)
BETAY=0.709             # twiss beta x in [m]  (T,B')=(25,23)

FREQ=408.e6             # common rf-frequency [Hz]
FREQ=816.e6             # common rf-frequency [Hz] (T>=25)

PHISY=-30.              # synchronous phase in [deg]

BGRAD=43.150            # quad gradient [T/m] (T=25)
BGRAD=23.000            # quad gradient [T/m] (T=25)
BGRAD=21.000            # quad gradient [T/m] (T=25)
# BGRAD=13.000            # quad gradient [T/m]
# BGRAD=6.7               # quad gradient [T/m]

NL=54                   # nboff lines a.k.a. cells
NL=125                  # (Ti,Tf)=(25,200)

MAP=t3d
# MAP=simple
MAP=oxal
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
ARGS="$ARGS -D _NL=$NL"
ARGS="$ARGS -D _MAPPING=$MAP"

# invoke m4
m4 $ARGS $1 > $2
