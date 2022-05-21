filename="model.3d.nz101.nx301.ny301.bin"
vplname="${filename}.vpl"
epsname="${filename}.eps"

echo n1=101 o1=0 d1=5\
 n2=301 o2=0 d2=5\
 n3=301 o3=0 d3=5\
 in=$filename data_format=native_float>rsftp

sfbyte < rsftp gainpanel=a pclip=100 |\
sfgrey3 frame1=0 frame2=100 frame3=100 color=g flat=n dframe=0\
 label1="depth (m)" label2="X (m)" label3="Y (m)" labelfat=2 gridfat=2\
 point1=0.65 point2=0.7 screenratio=1 \
 framelabel1=n	framelabel2=n framelabel3=n \
 framelabelcol=32767 font=3\
 >$vplname 

:<<!
echo n1=3000 d1=0.001 o1=0\
 n2=81 o2=-0.000825 d2=0.000015\
 n3=41 o3=-0.0003 d3=0.000015\
 in=$filename data_format=native_float>rsftp

sfbyte < rsftp gainpanel=a pclip=99. |\
 sfgrey3  wantaxis1=y axis1fat=20 frame1=0 frame2=57 frame3=20 color= flat=n\
 label1="Time (s)" label2="Ray parameter px (s/m)" label3="Ray parameter py (s/m)" labelfat=2\
 axiscol=7 axisfat=5000 wantaxis=y \
 point1=0.65 point2=0.7 screenratio=1 \
 framelabel1=n	framelabel2=n framelabel3=n \
 framelabelcol=32767 font=3\
 >$vplname 
!
vpconvert $vplname bgcolor=w color=y $epsname
