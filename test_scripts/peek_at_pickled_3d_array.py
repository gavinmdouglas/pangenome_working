import cPickle
import sys
from pprint import pprint


inalign = cPickle.load(open(sys.argv[1], 'rb'))
pprint(inalign[int(sys.argv[2]), int(sys.argv[3]), :])
