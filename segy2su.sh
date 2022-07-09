#find ./ -name "*.su*" -type f -print -exec rm -rf {} \;

paths=$(ls ./data)
for spath in ${paths}
do
folder1="./data/${spath}/*.segy"
folder="./data/${spath}/"
softfiles=$(ls $folder1)
for sfile in ${softfiles}
do 
    echo "${sfile}"
#    segyread tape=${folder}${sfile} verbose=0 endian=0>"${folder}${sfile}.su"
#    susort <"${folder}${sfile}.su">"${folder}${sfile}.su.sort" tracf
#   sustrip <"${sfile}.su.sort">"${sfile}.bin"
    filename="${sfile}.bin"
    epsname="${filename}.eps"
    figure="psimage"
    $figure< $filename \
    n1=4096 wbox=100 hbox=50 axeswidth=9 labelfont=Times-Roman \
    legend=0 lwidth=1 lheight=20 lstyle=vertright lx=22.0 \
    f1=0 d1=0.002 f2=0 d2=1 x2beg=9690 x2end=10982\
    label1='Time (s)' label2='Trace' \
    labelsize=75 \
    perc=98>$epsname
done
done