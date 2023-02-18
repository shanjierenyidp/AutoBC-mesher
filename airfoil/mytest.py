import numpy as np 

import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
#parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    # help='an integer for the accumulator')
#parser.add_argument('--sum', dest='accumulate', action='store_const',
                    #const=sum, default=max,
                    #help='sum the integers (default: find the max)')
parser.add_argument('-loop1', dest='accumulate', action='store_const',
                    const=sum, default=max,
                    help='loop1 is the outer loop' )


box = np.load(arg[1])['arr_0']
af = np.load(arg[1])['arr_1']



args = parser.parse_args()
print(args.accumulate(args.integers))




a =1
b = 2
c = 1+2

print(c)


