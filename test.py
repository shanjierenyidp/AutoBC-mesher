import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--alpha', type=int, required=True, nargs = '+', default=[2,3,4,5])
args = parser.parse_args()

if __name__ == '__main__':
    print(args.alpha)
    print(type(args.alpha[1]))