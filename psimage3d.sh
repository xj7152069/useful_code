filename="data3d.tp.bin"
vplname="${filename}.vpl"
epsname="${filename}.eps"
:<<!
echo n1=1024 d1=0.001 o1=0\
 n2=200 o2=0 d2=10\
 n3=200 o3=0 d3=10\
 in=$filename data_format=native_float>rsftp

sfbyte < rsftp gainpanel=a pclip=99|\
sfgrey3 frame1=0 frame2=200 frame3=0 color= flat=n\
 label1="t/s" label2="x/m" label3="y/m" labelfat=2 gridfat=2\
 point1=0.65 point2=0.7 screenratio=1 \
 >$vplname 
!

echo n1=1024 d1=0.001 o1=0\
 n2=80 o2=-0.000200 d2=0.000005\
 n3=50 o3=-0.000125 d3=0.000005\
 in=$filename data_format=native_float>rsftp

sfbyte < rsftp gainpanel=a pclip=99.95 |\
 sfgrey3 frame1=0 frame2=60 frame3=25 color= flat=n\
 label1="t/s" label2="px/(s/m)" label3="py/(s/m)" labelfat=2\
 point1=0.65 point2=0.7 screenratio=1 \
 >$vplname 


vpconvert $vplname bgcolor=w color=y $epsname
