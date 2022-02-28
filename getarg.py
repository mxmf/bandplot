# -*- coding: utf-8 -*-
import argparse


def get_args(argstr=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-hse",
        help="draw the band of zero weight kpoints",
        action="store_true"
    )
    parser.add_argument(
        "-yr",
        help="min and max energy value of y axis in the graph, default is [-4,6]",
        nargs=2,
        type=float,
        default=[-4, 6]
    )
    parser.add_argument(
        '-xr',
        help="min and max index of kpoints of x axis in the graph",
        nargs=2,
        type=int
    )
    parser.add_argument(
        "--fermi",
        '-f',
        help="If the parameter is an integer or a decimal, assign this parameter, if it is a file with a format similar to 'vasprun.xml', read from this file"
    )
    parser.add_argument(
        '-pband',
        help="plot projected energy band graph",
        action="store_true"
    )
    parser.add_argument(
        '-s',
        help="plot projected spin split energy band graph band for soc",
        action="store_true"
    )
    parser.add_argument(
        '-color',
        help="choose color map for proband",
        default=['hsv'],
        nargs='+'
    )
    parser.add_argument(
        "-atoms",
        help='choose which atoms you want to project',
        nargs='+',
        default=None
    )
    parser.add_argument(
        "-orbitals",
        help='choose which orbital you want to project',
        nargs='+',
        default=None
    )
    parser.add_argument(
        '--colorbar',
        '-cb',
        help="whether plot color bar",
        action="store_true",
    )
    if argstr:
        args = parser.parse_args(argstr.split())
    else:
        args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    print(args)
