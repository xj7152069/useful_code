#!/bin/bash

file_in=$1
file_rsf=${file_in}.rsf
title=` basename $1 `
#vmin=2000
#vmax=5500

nvz=2901
nvx=1300
nvy=470
zmax=5.0

frame1=1000
frame2=700
frame3=250

echo "n1=$nvz n2=$nvx n3=$nvy \
	d1=0.002 d2=1 d3=1\
	o1=0 o2=0 o3=0\
	in=$file_in \
	data_format=float_native"\
	>$file_rsf

	#minval=$vmin maxval=$vmax	\
<$file_rsf sfwindow max1=$zmax min1=0 | sfbyte gainpanel=all bar=gbar.rsf	\
	|\
	sfgrey3 flat=n frame1=$frame1 frame2=$frame2 frame3=$frame3 \
	point2=0.75 point1=0.6 \
	title=$title \
	label1='Time' unit1='s'  label2='CDP' label3='line' unit2='cdp-no' unit3='line-no' \
	labelfat=6 scalebar=y color=e  \
	|sfpen bgcolor=w &
	#labelfat=6 scalebar=y color=j  \
