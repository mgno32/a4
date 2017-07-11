#coding=utf-8
import cv2
import numpy as np
from sklearn.cluster import KMeans

def R(im):
    shp = im.shape
    rows, cols = shp[0:2]
    ratio = 512.0 / rows
    r = cv2.resize(im, (int(rows * ratio), int(cols * ratio))) 
    return r
PAPER_WIDTH = 1190.0
PAPER_HEIGHT = 1684.0
'''
PAPER_HEIGHT = 1190.0
PAPER_WIDTH = 1684.0
'''

def align(filename):
    im = cv2.imread(filename) 
    gray = cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
    gray = cv2.GaussianBlur(gray, (5,5), 0)
    cv2.imshow("gray",R(gray))
    edges = cv2.Canny(gray, 10, 100, apertureSize = 3)
    cv2.imshow("edges", R(edges))
    cv2.waitKey(0)
    lines = cv2.HoughLines(edges, 1, np.pi / 360, 650)
    print ("lines: ", len(lines))

    def merge_lines(lines):
        nlines = [line[0] for line in lines]
        km = KMeans(4)
        km.fit(nlines)
        return km.cluster_centers_

    lines = merge_lines(lines)

    print (lines)
    result = im.copy()
    h,w,c = result.shape
    for line in lines:
        # line: xcos(t) + ysin(y) - p = 0
        rho, theta = line
        _cos = np.cos(theta)
        _sin = np.sin(theta)
        if _cos != 0:
            pt1 = (int(rho / _cos), 0)
            pt2 = (int((rho - h * _sin) / _cos), h)
        else:
            #horizontal
            pt1 = (0, rho)
            pt2 = (w, rho)
        cv2.line(result, pt1, pt2, (0,0,255), 3)


    # 求交点
    def get_intersection(a,b):
        #a: (rho1, theta1)
        #b: (rho2, theta2)
        '''
        cs = np.array([[np.cos(a[1]), np.sin(a[1])], [np.cos(b[1]), np.sin(b[1])]])
        p = np.array([[a[0]],[b[0]]])
        # todo: x = 0
        w = np.dot(np.asmatrix(cs).I, p).T
        return np.array([w[0, 0], w[0, 1]])
        '''
        _cos0 = np.cos(a[1])
        _sin0 = np.sin(a[1])
        _cos1 = np.cos(b[1])
        _sin1 = np.sin(b[1])
        if _sin1 != 0:
            t = _sin0 / _sin1
            x = (a[0] - t * b[0]) / (_cos0 - _cos1 * t)
            y = (b[0] - x * _cos1) / _sin1
            return np.array([x,y])
        x = b[0] / _cos1
        y = (a[0] - x * _cos0) / _sin0
        return np.array([x,y])

    tts = np.vstack([get_intersection(lines[3], lines[j]) for j in range(3)])
    rec = np.sum(np.multiply(tts, tts), 1)
    ai = np.argmax(rec)
    #lines[3] 和 lines[ai] 几乎平行
    pts = []
    for j in range(3):
        if j != ai:
            pts.append(tts[j])
    pts.extend([get_intersection(lines[ai], lines[j]) for j in range(3) if j != ai])
    for pt in pts:
        xy =  (int(pt[0]), int(pt[1]))
        cv2.circle(result, xy, 32, (255,0,0), 32)

    cnts = [0] * 4
    for i in range(4):
        for j in range(4):
            if pts[i][0] > pts[j][0]:
                cnts[i] += 1
            if pts[i][1] > pts[j][1]:
                cnts[i] += 1

    lt = np.argmin(cnts)
    rb = np.argmax(cnts)

    # 计算右上角
    maxdx = -np.inf
    rt = -1
    for i in range(4):
        if i != lt and i != rb:
            dx = pts[i][0] - pts[lt][0]
            if dx > maxdx:
                maxdx = dx
                rt = i
    # 计算左下角
    for i in range(4):
        if i != lt and i != rb and i != rt:
            lb = i
    paper_p = np.array([(0,0), (PAPER_WIDTH,0), (PAPER_WIDTH, PAPER_HEIGHT), (0, PAPER_HEIGHT)]).astype(np.float32)
    pts = [pts[i].tolist() for i in range(4)]
    pts.sort()
    [lt, rt, rb, lb] = [1,3,2,0]
    sp = np.array([pts[i] for i in [lt,rt,rb,lb]]).astype(np.float32)
    #print (sp)
    M = cv2.getPerspectiveTransform(sp, paper_p)
    paper = cv2.warpPerspective(im, M, (int(PAPER_WIDTH), int(PAPER_HEIGHT)))
    result = cv2.resize(result, (1080,720))
    cv2.imshow("source", R(im))
    cv2.waitKey(0)
    cv2.imshow("lines", R(result))
    cv2.waitKey(0)
    cv2.imshow("paper", R(paper))
    cv2.waitKey(0)
    return paper
