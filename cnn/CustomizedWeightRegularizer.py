class CustomizedWeightRegularizer(WeightRegularizer):
    def __init__(self, l1=0., l2=0.):
        self.l1 = K.variable(l1)
        self.l2 = K.variable(l2)
        self.uses_learning_phase = True
        self.p = None