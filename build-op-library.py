import argparse
import numpy as np

parser = argparse.ArgumentParser(description="A simple CLI tool")
parser.add_argument("--radius_low", type=int, help="Smallest grain/bubble you want to include in the OP library")
parser.add_argument("--radius_high", type=int, help="Largest grain/bubble you want to include in the OP library")
parser.add_argument("--resolution", type=int, help="the step size between lower and upper radii to use when creating files.")
parser.add_argument("--ice_type", type=str, help="Type of ice to include in library: Grains, slab or both")

args = parser.parse_args()

def print_args(args):
    """
    placeholder func that will eventually create library of optical property files according to the cli args.
    currently just echoes cli args back to console in a loop
    """

    print("Creating files:")
    for i in np.arange(args.radius_low, args.radius_high+args.resolution, args.resolution):
        if args.ice_type == 'both':
            print(f"    ice_grain_{i}_um.npz")
            print(f"    ice_bubble_{i}_um.npz")
        elif args.ice_type == 'grains':
            print(f"    ice_grain_{i}_um.npz")
        elif args.ice_type == 'slab':
            print(f"    ice_bubble_{i}_um.npz")
        else:
            print("ice type not recognized, please select grains, slab or both,")
    return

print_args(args)
