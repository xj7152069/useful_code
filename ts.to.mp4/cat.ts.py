
import os

#ts文件绝对路径
ts_path = '/home/dell/桌面/children/output'  
#读取ts文件夹下所有的ts文件
path_list = os.listdir(ts_path)
#对文件进行排序
path_list.sort()
#将排序后的ts的绝对路径放入列表中
li = [os.path.join(ts_path,filename) for filename in path_list]
#类似于[001.ts|00.2ts|003.ts]
input_file = '|'.join(li)
#指定输出文件名称
output_file = ts_path + '/combine.mp4'
#使用ffmpeg将ts合并为mp4
command = 'ffmpeg -i "concat:%s" -acodec copy -vcodec copy -absf aac_adtstoasc %s'%	(input_file,output_file)
#指行命令
os.system(command)

