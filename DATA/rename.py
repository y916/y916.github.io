import os
path = r"H:\ABAQUS\DATA\JZT22"
filelist = os.listdir(path) #该文件夹下所有的文件（包括文件夹）

for file in filelist:
    print(file)
for file in filelist:   #遍历所有文件
    Olddir=os.path.join(path,file)   #原来的文件路径
    if os.path.isdir(Olddir):   #如果是文件夹则跳过
        continue
    filename=os.path.splitext(file)[0]   #文件名
    filetype=os.path.splitext(file)[1]   #文件扩展名
    Newdir=os.path.join(path,filename+'_加钢板'+filetype)
    os.rename(Olddir,Newdir)#重命名
