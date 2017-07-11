from sklearn.svm import SVC
from sklearn.externals import joblib

import cv2
import numpy as np
import random
fin = open("./dataset/label.txt")

X = []
y = []
for line in fin.readlines():
    name, label = [int(d) for d in line.split(',')]
    pic = cv2.imread("./dataset/%s.png" % name)
    gray = cv2.cvtColor(pic, cv2.COLOR_BGR2GRAY)
    X.append(gray.flatten())
    y.append(label)


# create_new_data

ylen = len(y)
ylen = 0
for i in range(ylen):
    x = X[i].reshape((100, 100))
    t = y[i]
    rows, cols = x.shape[:2]
    for k in range(5):
        H = np.float32([[1, 0, random.randint(1, 20) - 10], [0, 1, random.randint(1, 20) - 10]])
        res = cv2.warpAffine(x, H, (rows, cols))
        X.append(res.flatten())
        y.append(t)
        M = cv2.getRotationMatrix2D((cols/2, rows/2), random.randint(-30, 30), 1.0)
        res = cv2.warpAffine(x, M, (rows, cols))
        X.append(res.flatten())
        y.append(t)
    print ("%d/%d" % (i, ylen))

X = np.array(X)
y = np.array(y)
print ("predict")
print (X.shape, y.shape)

clf = SVC()
clf.fit(X, y)
pred = clf.predict(X)

b = (y == pred)
print ("Accuary: %f" % (np.sum(b) * 1.0 / b.size))
print (y[~b], pred[~b])
joblib.dump(clf, "svm_model.m")
