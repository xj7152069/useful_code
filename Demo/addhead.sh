set -x
suaddhead < shot_cnooc_fault_fd_all.pc.dat ns=800 | sushw key=dt a=1000 | \
sushw key=fldr,tracf a=1,1 b=0,1 c=1,0 j=101,101|\
sushw key=sx,gx a=500,0 b=0,10 c=10,10 j=101,101|\
suchw key1=offset key2=gx key3=sx b=1 c=-1|\
sushw key=cdp a=1 b=1 c=2 j=101> yang_model_shot.su 

#a:initial value in ensemble
#b:increment in ensemble
#c:increment between ensembles
#j:traces in a ensemble

#suchw: calculate some unknown header with some known headers
#       key1=a+b*key2+c*key3

