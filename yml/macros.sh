#TKIN  = kinetic energy in
#EMITW = w emittance
#EMITX = x emittance
#EMITY = y emittance
#BETAX = beta x
#BETAY = beta y
#PHISY = synchronous phase

TKIN=50.
##EMITW=12.6e-6
##EMITW=12.6e-5
EMITW=9.6e-5
##EMITX=1.e-6
##EMITY=1.e-6
EMITX=2.e-6
EMITY=2.e-6
BETAX=2.285
BETAY=0.328
PHISY=-20.

# invoke m4
m4 -D_TKIN=$TKIN -D_EMITW=$EMITW -D_EMITY=$EMITY -D_EMITX=$EMITX -D_BETAX=$BETAX -D_BETAY=$BETAY -D_PHISY=$PHISY $1
