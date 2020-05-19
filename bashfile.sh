<<'COMMENT'
批量注释：
    参数含义：
    运行时，使用命令： nohup bash XXX.sh &  ；可以回车挂起运行。

#编译
g++ filename.cpp -o filename.cx -larmadillo

#输入参数运行程序
echo par1 par2 \
    par3 par4 \
    str_filename \ | ./filename.cx

#用循环并行程序运行：共运行5*20次
for((j=0;j<=5;j++))
do

#内循环一次运行20个程序，要用20个线程
for ((i=1;i<=20;i++))
do
c++ filename.cpp -o filename_$[i].cx
    echo \
    par1 \
    par2 \
    $[$[i]*5+$[j]*20*5+95]_filename | nohup ./filename_$[i].cx & #后台运行程序，将输出存入'nohup'
done
wait #等待内循环调用的程序运行完成

done
COMMENT

g++ wave.poynting.hilbert.cpp -o test.parallel

echo 150 150 10 | nohup ./test.parallel &
echo 200 200 10 | nohup ./test.parallel &
echo 250 250 10 | nohup ./test.parallel &

ps -e | grep 'test.parallel' | awk '{print $1}'
for pid in `ps -e|grep 'test.parallel'|awk '{print $1}'`;
do
taskset -pc 5,6,7 $pid;
done

<<'COMMENT'
echo 50 90 10 | nohup ./test.parallel &
echo 100 140 10 | nohup ./test.parallel &
echo 150 190 10 | nohup ./test.parallel &
echo 200 240 10 | nohup ./test.parallel &
echo 250 290 10 | nohup ./test.parallel &
echo 300 350 10 | nohup ./test.parallel &
COMMENT
