from sklearn.svm import SVC
from sklearn.externals import joblib

clf = joblib.load("svm_model.m")
def predict(im):
    return clf.predict(im.flatten())[0]
