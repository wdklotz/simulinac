import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('integers', metavar='N',       type=int,             nargs='+',                 help='an integer for the accumulator')
parser.add_argument('--sum',    dest='accumulate', action='store_const', const=sum,   default=max,  help='sum the integers (default: find the max)')
parser.add_argument('--p')
parser.add_argument("--hide", action="store_true",               help="hide IN/OUT scatter plots")
# parser.add_argument("--hide",             help="hide IN/OUT scatter plots")

args = parser.parse_args()
print(args.accumulate(args.integers))
print(args)