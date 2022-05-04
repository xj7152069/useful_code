folder="./"
softfiles=$(ls $folder)
for sfile in ${softfiles}
do 
    echo "${sfile}"
    segyhdrs <${sfile}
    segyname="./segy/${sfile}.segy"
    segywrite <${sfile} tape=$segyname endian=0
done
