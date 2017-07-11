import numpy as np
from mobula import Net
from mobula.layers import FC
from mobula.layers import Sigmoid, CrossEntropy
from mobula.layers import Data

FSIZE = 20
X = np.zeros((1,1,FSIZE,FSIZE))
Y = np.zeros((1, 10))

data = Data(X, "Data", batch_size = 128, label = Y)

fc1 = FC(data, "fc1", dim_out = 25)
sig1 = Sigmoid(fc1, "sig1")
fc2 = FC(sig1, "fc2", dim_out = 10)
sig2 = Sigmoid(fc2, "sig2")
loss = CrossEntropy(sig2, "Loss", label_data = data)

net = Net()
net.setLoss(loss)

net.reshape()
fc1.W,fc1.b,fc2.W,fc2.b = np.load("./mnist_net.npy", encoding="bytes")
def predict_sub(pic):
    #pic.shape = (20, 20)
    x = pic.reshape((1,1,FSIZE,FSIZE)) / 255.0
    data.allY = x
    net.forward()
    pre = np.argmax(sig2.Y,1)
    print (sig2.Y)
    return pre[0], np.max(sig2.Y)

def predict(pz):
    #return predict_sub(pz)
    pre1, p1 = predict_sub(pz.T)
    pre2, p2 = predict_sub(pz[::-1, :].T)
    print (p1, p2, "==", pre1, pre2)
    if p1 > p2:
        return pre1
    return pre2

'''
net.reshape2()

net.lr = 0.3
for i in range(100000):
    net.forward()
    net.backward()

    if i % 100 == 0:
        print "Iter: %d, Cost: %f" % (i, loss.Y)
        old_batch_size = data.batch_size
        data.batch_size = None
        net.reshape()
        net.forward()
        pre = np.argmax(sig2.Y,1)
        pre.resize(pre.size)
        right = np.argmax(data.label, 1).reshape(pre.size)
        bs = (pre == right) 
        b = np.sum(bs)
        acc = (b * 1.0 / len(pre))
        print (pre[0:5000:50])
        print ("Accuracy: %f" % (acc))
        if b == len(pre):
            np.save("mnist_net.npy", [fc1.W, fc1.b, fc2.W, fc2.b])
            print ("Save OK")
            import sys
            sys.exit()
        data.batch_size = old_batch_size
        net.reshape()


np.save("mnist_net_over.npy", [fc1.W, fc1.b, fc2.W, fc2.b])
print ("Save Over OK")
'''
