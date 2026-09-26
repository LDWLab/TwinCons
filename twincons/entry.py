#!/usr/bin/env python3
"""Entry point script for twcons command line tool."""

import sys


def main():
    """Main entry point that imports and runs TwinCons.main()."""
    from twincons.TwinCons import main as twcons_main
    sys.exit(twcons_main(sys.argv[1:]))


if __name__ == '__main__':
    main()