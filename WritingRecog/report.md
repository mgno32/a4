#手写数字识别-实验报告
14301013 贺小雨

##实验环境
操作系统: Arch Linux 4.10.3
编程语言：Python2.7
库：opencv numpy sklearn
##实验流程

###图像矫正
图像矫正部分我使用了
###数字分割

###数字识别-SVM
####关于训练集的说明：
* 我将miniset和分割得到的结果，使用同一个模型，进行了预测准确度比较：
|测试集   |训练集  | 准确率  |
|----   |----   |----  |
|miniset|miniset|98.6%|
|miniset|图片1  | 50.0%  |
|图片1 2 3 5| 图片4| 0%| 

在使用miniset作为训练集和测试集时，本模型精度可以达到98.6%，但是预测图片1的精度就很低了。我觉得这是因为miniset里面的写法和图片1中的写法不一样。如果使用其他四个图片去预测另一个图片，也会发现SVM过拟合问题非常严重，导致一个都不对。
因此为了得到满意的结果，我将测试集中的图片加入了训练集。
* 
##运行结果展示
|测试集  |训练集 |准确率  |
进入ocr文件夹
> python2 test.py 1
> python2 test.py 2
> python2 test.py 3
> python2 test.py 4
> python2 test.py 5
