from numpy import genfromtxt
import gzip
import pickle as cPickle

dirg = 'fold/' 

for i in range(1):

  k = i;
  name1 = dirg+'/train'+str(k)+'_x'
  name2 = dirg+'/train'+str(k)+'_y'
  train_set= genfromtxt(name1, delimiter=','),  genfromtxt(name2, delimiter=',')

  name1 = dirg+'/test'+str(k)+'_x'
  name2 = dirg+'/test'+str(k)+'_y'
  test_set = genfromtxt(name1, delimiter=','), genfromtxt(name2, delimiter=',')

  name1 = dirg+'/val'+str(k)+'_x'
  name2 = dirg+'/val'+str(k)+'_y'
  val_set = genfromtxt(name1, delimiter=','), genfromtxt(name2, delimiter=',')
 
  dataset = [train_set, test_set, test_set]
  del train_set
  del val_set
  del test_set

  name = dirg+'mpqa_spa'+str(k)+'.pkl.gz' 
  f = gzip.open(name,'wb')
  cPickle.dump(dataset, f, protocol=2)
  f.close()
