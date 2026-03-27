import argparse

def get_parser():
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter_class)
    subparsers = parser.add_subparsers(dest='command', required=True)

    p1 = subparsers.add_parser('train', help='train')
    p1.add_argument('--input', type=str, default=None, help='input data')

    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    if args.command == 'train':
        pass