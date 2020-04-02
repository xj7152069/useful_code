<<'COMMENT'
批量注释：
    参数含义：
    运行时，使用命令： nohup bash XXX.sh &  ；可以回车挂起运行。

COMMENT

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
