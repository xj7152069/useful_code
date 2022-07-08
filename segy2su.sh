#find ./ -name "*.su*" -type f -print -exec rm -rf {} \;

paths=$(ls ./data)
for spath in ${paths}
do
folder1="./data/${spath}/*.segy"
folder="./data/${spath}/"
echo $folder
softfiles=$(ls $folder)
for sfile in ${softfiles}
do 
    echo "${folder}${sfile}"
    segyread tape=${folder}${sfile} verbose=0 endian=0>"${folder}${sfile}.su"
    susort <"${folder}${sfile}.su">"${folder}${sfile}.su.sort" tracf
    #segyhdrs <${sfile}
    #segyname="./segy/${sfile}.segy"
    #segywrite <${sfile} tape=$segyname endian=0
done
done