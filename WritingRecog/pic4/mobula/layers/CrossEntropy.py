from .Layer import *

class CrossEntropy(Layer):
    def __init__(self, model, *args, **kwargs):
        Layer.__init__(self, model, *args, **kwargs)
        self.label = kwargs.get("label") # (N, C)
        self.label_data = kwargs.get("label_data")
    def __str__(self):
        return "It is a Mean Squared Error Layer"
    def reshape(self):
        self.Y = 0.
    def forward(self):
        if self.label_data:
            self.label = self.label_data.label
        self.Y = np.mean(- np.multiply(self.label, np.log(self.X)) - \
               np.multiply(1.0 - self.label, np.log(1.0 - self.X)))
    def reshape2(self):
        self.dX = np.zeros(self.X.shape) 
    def backward(self):
        self.dX = -self.label / self.X + (1.0 - self.label) / (1.0 - self.X)
