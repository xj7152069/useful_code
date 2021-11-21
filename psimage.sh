filename="real.rebuild.bin"
epsname="${filename}.eps"
#figure="pswigb"
figure="psimage"

$figure< $filename \
n1=3000 wbox=40 hbox=60 \
legend=0 f1=-0 d1=0.001 f2=-00 d2=10 \
x1beg=0 x1end=1 \
x2beg=0 x2end=1000 \
label1='Time(s)' label2='X(m)' \
labelsize=140 \
perc=99.9>$epsname

:<<!
$figure< $filename \
n1=1024 wbox=20 hbox=30 \
legend=0 f1=-0 d1=0.001 f2=-0.0000995 d2=0.000001 \
x1beg=0 x1end=1 \
x2beg=-0.000051 x2end=0.000051 \
label1='Time(s)' label2='p(s/m)' \
labelsize=70 \
bclip=0.00268943 wclip=-0.001421>$epsname
!
#x2beg=0 x2end=1000\
#x2beg= x2end=\
#bclip=1.24932 wclip=-0.606472\
#lstyle=vertright lwidth=0.2 lheight=3 \
#bclip=0.001 wclip=-0.001 x1beg=500 x1end=3000 \
#threecolor=1 bps=24 \


#red white blue: wrgb=0,0,1.0 grgb=1.0,1.0,1.0 brgb=1.0,0,0
