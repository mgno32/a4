#coding=utf-8
import numpy as np
import cv2
import net
from align import *

FSIZE = 100
LINE_H = 30
filename = "../4.jpg"
im = align(filename) 
edges = cv2.Canny(im, 50,500, apertureSize = 3)
kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (5,5))
edges = cv2.morphologyEx(edges, cv2.MORPH_CLOSE, kernel)
#kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (5,5))
#edges = cv2.morphologyEx(edges, cv2.MORPH_OPEN, kernel) # 去噪声
cv2.imshow("canny", R(edges))
cv2.waitKey(0)

b = edges > 0
b = b.astype(np.uint8) * 255
valve = 170
bf = ((im[:,:,0] < valve) & (im[:,:,1] < valve) & (im[:,:,2] < valve)).astype(np.uint8) * 255


#kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (2,2))
#b = cv2.morphologyEx(b, cv2.MORPH_OPEN, kernel) # 去噪声
kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (5,5))
b = cv2.morphologyEx(b, cv2.MORPH_CLOSE, kernel)

b, contours, hierarchy = cv2.findContours(b, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
bc = im.copy()#np.dstack([bf,bf,bf])
cs = []
for c in contours:
    minc = np.min(c, 0).reshape(2)
    maxc = np.max(c, 0).reshape(2)
    cs.append([minc, maxc])

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
#cv2.imshow("gray", gray)
#cv2.waitKey(0)

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
    ratio = w * 1.0 / h
    if ratio > 1.5:
        continue
    #it is digit
    bc = cv2.rectangle(bc, pt1, pt2, (255,0,255), 1) 
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

    import os
    path = "/home/hal/zk/HelloShirley/ocr/dataset/"
    name = path + "%s.png" % 0
    k = 0
    while os.path.exists(name): 
        k += 1
        name = path + "%s.png" % k
    cv2.imwrite(name, pz)
    
    '''
    pre = net.predict(pz)
    py = minc[1]
    if lastY is not None and py - lastY > LINE_H:
        pres.append([])
        layer = pres[-1]
    lastY = py 
    layer.append(pre)

    '''

print ("\n\n\n识别结果：")
for i in range(len(pres)):
    print ("Line {}: {}".format(i, pres[i]))
    '''
    cv2.imshow("pz",pz)
    cv2.waitKey(500)
    cv2.waitKey(0)
    '''
cv2.imshow("im", R(bc))
cv2.waitKey(0)
