name1='seg-'
name2='-v1-a1.ts'
id1=1
id2=123

migration='fad897a3296e30c1e9817de351cab044'
key='C21F969B5F03D33D43E04F8F136E7682'

echo $migration
echo $key
wait

k=10000

for((j=$[id1];j<=$[id2];j++))
do

echo $[j]
openssl aes-128-cbc -d -in $name1$[j]$name2 -out ./output/$[j+k].ts -nosalt -iv $migration -K $key
wait

done
wait


