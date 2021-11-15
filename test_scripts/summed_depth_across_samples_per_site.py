import cPickle
import sys
from pprint import pprint
import numpy as np

inalign = cPickle.load(open(sys.argv[1], 'rb'))
pprint(np.sum(inalign[:, int(sys.argv[2]), 0]))
pprint(np.sum(inalign[:, int(sys.argv[2]), 1]))
pprint(np.sum(inalign[:, int(sys.argv[2]), 2]))
pprint(np.sum(inalign[:, int(sys.argv[2]), 3]))

