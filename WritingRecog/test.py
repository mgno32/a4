#coding=utf-8
import numpy as np
import cv2
#import net
import svm
import sys
from align import *

FSIZE = 100
LINE_H = 100
imid = 5
if len(sys.argv) > 1:
    imid = int(sys.argv[1])
filename = "./%d.jpg" % imid

def is_pic(num):
    return ("%d.jpg" % num) in filename

im, source = align(filename) 
if is_pic(4):
    edges = cv2.Canny(im, 50,500, apertureSize = 3)
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (5,5))
    edges = cv2.morphologyEx(edges, cv2.MORPH_CLOSE, kernel)
elif is_pic(5):
    edges = cv2.Canny(im, 50,300, apertureSize = 3)
else:
    edges = cv2.Canny(im, 100,200, apertureSize = 3)
kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (5,5))
edges = cv2.morphologyEx(edges, cv2.MORPH_CLOSE, kernel)
#kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (5,5))
#edges = cv2.morphologyEx(edges, cv2.MORPH_OPEN, kernel) # 去噪声

b = edges > 0
b = b.astype(np.uint8) * 255
valve = 170
bf = ((im[:,:,0] < valve) & (im[:,:,1] < valve) & (im[:,:,2] < valve)).astype(np.uint8) * 255

'''
kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (3,3))
b = cv2.morphologyEx(b, cv2.MORPH_OPEN, kernel)
'''
kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (5,5))
b = cv2.morphologyEx(b, cv2.MORPH_CLOSE, kernel)
#寻找联通分量
b, contours, hierarchy = cv2.findContours(b, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

bc = im.copy()#np.dstack([bf,bf,bf])
delete_ratio = []
cs = []
areas = []
ws = []
hs = []
for c in contours:
    minc = np.min(c, 0).reshape(2) #(min(x), min(y))
    maxc = np.max(c, 0).reshape(2)
    h = abs(maxc[1] - minc[1])
    w = abs(maxc[0] - minc[0])
    ratio = 1.0 * w / h
    #
    multi = 20
    if (w > 5 and h > 10  and w < 600 and h <600 and  ratio > 1.0/multi and ratio < multi ): 
        print (ratio)
        area = abs(maxc[0] - maxc[1]) * abs(minc[1] - minc[0])
        areas.append(area)
        delete_ratio.append(c)
        ws.append(w)
        hs.append(h)

avg_area = sum(areas) / len(areas)
def is_number_rect(avg_area,area,minc,maxc):
    return True
    multi = 20
    if (area < (1.0 / multi) * avg_area or area > multi * avg_area):
        return False
    return True

if is_pic(4):
    for c in contours:
        minc = np.min(c, 0).reshape(2)
        maxc = np.max(c, 0).reshape(2)
        cs.append([minc, maxc])
else:
    for c in delete_ratio:
        minc = np.min(c, 0).reshape(2) #(min(x), min(y))
        maxc = np.max(c, 0).reshape(2)
        area = abs(maxc[0] - maxc[1]) * abs(minc[1] - minc[0])
        if is_number_rect(avg_area,area,minc,maxc):
            print ("ok")
            cs.append([minc,maxc])
        print ("——————————————————————————")

def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K(object):
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0  
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K

def cmpxy(a,b):
    dy = a[0][1] - b[0][1]
    if abs(dy) < LINE_H:
        if a[0][0] < b[0][0]:
            return -1
        if a[0][0] > b[0][0]:
            return 1
    else:
        if a[0][1] < b[0][1]:
            return -1
        if a[0][1] > b[0][1]:
            return 1
    return 0

gray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY) 
cs.sort(key = cmp_to_key(cmpxy)) 
pres = [[]]
layer = pres[0] 
lastY = None
print ("num?: ", len(cs))
for minc, maxc in cs:
    pt1 = (minc[0], minc[1])
    pt2 = (maxc[0], maxc[1])
    nz = bf[minc[1]:maxc[1], minc[0]:maxc[0]]
    w = pt2[0] - pt1[0]
    h = pt2[1] - pt1[1]
    '''
    if w > 100 or h > 100:
        continue
    if h < 12:
        continue
    '''
    if is_pic(4):
        ratio = w * 1.0 / h
        if ratio > 1.5:
            continue
        if 1.0 / ratio > 15:
            continue
    #it is digit
    bc = cv2.rectangle(bc, pt1, pt2, (255,0,255), 10) 
    if w > h:
        nw = FSIZE
        nh = FSIZE * h // w
        x = 0
        y = (FSIZE - nh) // 2
    else:
        nh = FSIZE
        nw = FSIZE * w // h
        y = 0
        x = (FSIZE - nw) // 2
    if nh < 1 or nw < 1:
        continue
    nz = cv2.resize(nz, (nw, nh))
    pz = np.zeros((FSIZE, FSIZE))
    #nz = np.pad(nz, ((1,1),(2,2)), "constant")
    pz[y:y+nh, x:x+nw] = nz

    pre = svm.predict(pz) 
    py = minc[1]
    if lastY is not None and py - lastY > LINE_H:
        pres.append([])
        layer = pres[-1]
    lastY = py 
    layer.append(pre)
    '''
    import os
    path = "/home/hal/Downloads/Final+Project/dataset/"
    name = path + "%s.png" % 0
    k = 0
    while os.path.exists(name): 
        k += 1
        name = path + "%s.png" % k
    '''
    #cv2.imwrite(name, pz)
print ("\n\n\n识别结果：")
for i in range(len(pres)):
    print ("Line {}: {}".format(i, pres[i]))

cv2.imshow("canny", R(edges))
cv2.imshow("im", R(bc))
cv2.imshow("source", R(source)) 
cv2.waitKey(0)

