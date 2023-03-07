import argparse

parser = argparse.ArgumentParser(description='Process inner loops')
parser.add_argument('loops', metavar='loop', type=str, nargs='+',
                    help='a list of inner loops in the format of "points,resolution,segments"')
args = parser.parse_args()

for loop in args.loops:
    points, resolution, segments = loop.split(',')
    # process the inner loop here
    print(f"Processing inner loop with points={points}, resolution={resolution}, segments={segments}")
