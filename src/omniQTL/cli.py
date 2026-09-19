"""Command-line entry point for omniQTL."""

import argparse


def get_parser() -> argparse.ArgumentParser:
    """Build the omniQTL command-line argument parser.

    Returns:
        The configured `argparse.ArgumentParser`, with a required
        ``command`` subparser (currently only ``train``).
    """
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter_class)
    subparsers = parser.add_subparsers(dest='command', required=True)

    p1 = subparsers.add_parser('train', help='train')
    p1.add_argument('--input', type=str, default=None, help='input data')

    return parser


def main() -> None:
    """Parse command-line arguments and dispatch to the requested command."""
    parser = get_parser()
    args = parser.parse_args()
    if args.command == 'train':
        pass
