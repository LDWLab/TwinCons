"""Entry point for the twcons command line tool."""

import sys


def main():
    """Runs TwinCons on the command line arguments."""
    from twincons.TwinCons import main as twcons_main
    # main() returns the scores for -r, which are meant for Python callers, not an exit status.
    twcons_main(sys.argv[1:])
    return 0


if __name__ == '__main__':
    sys.exit(main())
