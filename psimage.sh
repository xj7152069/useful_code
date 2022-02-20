filename="hessmat.simple.bin"
epsname="${filename}.eps"
#figure="pswigb"
figure="psimage"
:<<!
$figure< $filename \
n1=2048 wbox=40 hbox=50 axeswidth=12 labelfont=Times-Roman \
legend=0 f1=-2.0 d1=0.004 f2=-2000 d2=40.0 \
x1beg=0 x1end=2 \
x2beg=0 x2end=4000 \
label1='t / s' label2='x / m' \
labelsize=140 \
perc=99.9>$epsname
!
:<<!
$figure< $filename \
n1=2048 wbox=40 hbox=50 axeswidth=12 labelfont=Times-Roman \
legend=0 f1=-2 d1=0.004 f2=-0.00125 d2=0.000005 \
x1beg=-0 x1end=2 \
x2beg=-0.00125 x2end=0.00125 \
label1='t / s' label2='p / (s/m)' \
labelsize=140 \
perc=99.998>$epsname
!

$figure< $filename \
n1=401 wbox=20 hbox=20 axeswidth=9 labelfont=Times-Roman \
legend=1 lwidth=1 lheight=20 lstyle=vertright lx=22.0 \
f1=-0.001999999 d1=0.00000999 f2=-0.001999999 d2=0.00000999 \
label1='p / (s/m)' label2='p / (s/m)' \
labelsize=75 \
bclip=200 wclip=-50.0>$epsname

#x2beg=0 x2end=1000\
#x2beg= x2end=\
#bclip=1.24932 wclip=-0.606472\
#lstyle=vertright lwidth=0.2 lheight=3 \
#bclip=0.001 wclip=-0.001 x1beg=500 x1end=3000 \
#threecolor=1 bps=24 \


#red white blue: wrgb=0,0,1.0 grgb=1.0,1.0,1.0 brgb=1.0,0,0
