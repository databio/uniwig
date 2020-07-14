"""
Internal CLI for smugwig internal processes. This provides a CLI to C utilities
included with smugwig
."""

import argparse
import sys
import logmuse
import logging
import smugwigc
_LOGGER = None

def build_argparser():
    parser = argparse.ArgumentParser("Converts specific site coordinates to wiggle format.")

    parser.add_argument("--mode", choices=["exact", "smooth", "echo"],
        help = "Choose whether you want specific sites, or smoothed.")
    parser.add_argument("--step-type", default="variable", choices=["fixed", "variable"],
        help="Choose fixed or variable (default) step format.")
    parser.add_argument("--step-size", type=int, default=1,
        help="Wiggle step size. Higher numbers make smaller files with lower resolution")
    parser.add_argument("--smooth-length", type=int, default=25,
        help="Wiggle smooth size. Higher numbers smooth more")
    parser.add_argument("--chromsize", type=int,
        help="Size (length) of the chromosome being processed.")

    # CLI hooks for logmuse
    parser = logmuse.add_logging_options(parser)
    return parser

def main():
    parser = build_argparser()
    args = parser.parse_args()
    global _LOGGER
    _LOGGER = logmuse.logger_via_cli(args)
    # _LOGGER = logging.getLogger(__name__)
    _LOGGER.debug("Welcome to smugwig")

    if args.step_type == "variable":
        if args.chromsize and args.chromsize != 0:
            _LOGGER.debug("chromsize argument is unused when step-type is 'variable'")
        else:
            # set a default
            args.chromsize = 0
        if args.step_size != 1:
            _LOGGER.debug("step-size argument is unused when step-type is 'variable'")

    if args.step_type == "fixed":
        if not args.chromsize or args.chromsize == 0:
            _LOGGER.error("You must provide a chromsize if using fixed step.")
            raise Exception


    if args.mode == "echo":
        # _LOGGER.debug("Smugwig echoing...")
        smugwigc.echo()
    elif args.mode == "exact":
        # _LOGGER.debug("Smugwig building exact wiggle...")
        smugwigc.sitesToExactWig(args.step_type == "variable", args.chromsize, args.step_size)
    elif args.mode == "smooth":
        # _LOGGER.debug("Smugwig building smooth wiggle...")
        smugwigc.sitesToSmoothWig(args.step_type == "variable",
            args.chromsize, args.step_size, args.smooth_length)
    else:
        parser.print_help(sys.stderr)

if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        _LOGGER.error("Program canceled by user!")
        sys.exit(1)
