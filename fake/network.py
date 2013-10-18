import numpy as np

class Cubic:

  def __init__(self, a=1, b=0.1, c=[1,2,3], d='name'):
    self.values = np.random.uniform(0,1,[30,30,30])

  def __str__(self):
    return "<Cubic {self.values.shape}>".format(**locals())

if __name__ == '__main__':
  cubic = Cubic()
  print( cubic )