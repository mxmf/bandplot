# -*- coding: utf-8 -*-
from pathlib import Path
import argparse
import matplotlib as mpl


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'rows',
        help='args.rows of subfig',
        type=int
    )
    parser.add_argument(
        'cols',
        help='number of columns of subfig',
        type=int
    )
    parser.add_argument(
        '-nk',
        help='no kpoints file',
        action='store_true',
        default=False
    )
    parser.add_argument(
        '-marker',
        help='whether marker a,b,c,d for every subgraph',
        action='store_true',
        default=False
    )
    parser.add_argument(
        '--title',
        '-t',
        help='make title for every subgraph',
        action='store_true',
        default=False
    )
    parser.add_argument(
        '-tt',
        help='Only set the title on the first line',
        action='store_true',
        default=False
    )
    parser.add_argument(
        '-xlabel',
        help='make xlabel for every subgraph',
        action='store_true',
        default=False
    )
    parser.add_argument(
        '-lb',
        help='only marker xlabel, xtickslabel, ylabel, ytickslabel on the left and bottom ',
        action='store_true',
        default=False
    )
    args = parser.parse_args()
    return args


def generate_subfig_plot(file='plot_subfig.py', args=None):
    cur_file_path = Path(__file__).resolve().parent
    with open(file, 'w') as fp:
        fp.write('# -*- coding: utf-8 -*-\n')
        fp.write('import sys\n')
        fp.write('import matplotlib\n')
        fp.write('import matplotlib.pyplot as plt\n')
        fp.write("sys.path.append('{}')\n".format(cur_file_path))
        fp.write('from plotband import PlotBand, get_args\n')
        fp.write('from readfile import ReadKpoints\n\n')

        fp.write("if __name__ == '__main__':\n")
        fp.write(
            "    matplotlib.rc_file('{}/matplotlibrc')\n\n".format(cur_file_path))

        fp.write('    fig, axs = plt.subplots(figsize=({}, {}), nrows={}, ncols={}, dpi={})\n\n'.format(
            mpl.rcParams['figure.figsize'][0]*args.cols, mpl.rcParams['figure.figsize'][1]*args.rows, args.rows, args.cols, 100/max(args.rows, args.cols)))
        fp.write("    #the h_pad and wpad are Padding around axes objects, the hspace and wspace are Space between subplot groups.\n")
        fp.write(
            "    fig.set_constrained_layout_pads(w_pad=0.04167, h_pad=0.04167, hspace=0.02, wspace=0.02)\n\n")
        if args.cols == 1 and args.rows == 1:
            fp.write('    axs=[axs]\n\n')
        else:
            fp.write('    axs=axs.flat\n\n')

        if args.nk:
            fp.write(
                "    #the numbers after {file/symbols/args} represent [row num] and [column num]. i.e. file11 represent the file at row1 and col1.\n")
        else:
            fp.write(
                "    #the numbers after {file/kpoints/args} represent [row num] and [column num]. i.e. file11 represent the file at row1 and col1.\n")

        for row in range(args.rows):
            for col in range(args.cols):
                fp.write("    file{}{} = 'vasprun.xml'\n".format(row+1, col+1))
        fp.write('\n')
        if args.nk:
            for row in range(args.rows):
                for col in range(args.cols):
                    fp.write(
                        "    symbols{}{} = ['M', 'K', '$\Gamma$', 'K\_', 'M']\n".format(row+1, col+1))
            fp.write(
                '\n    #Parameters for different subgraphs, the form is the same as the command line\n')
        else:
            for row in range(args.rows):
                for col in range(args.cols):
                    fp.write(
                        "    kpoints{}{} = ReadKpoints('KPOINTS')\n".format(row+1, col+1))
            fp.write(
                '\n    #Parameters for different subgraphs, the form is the same as the command line\n')

        for row in range(args.rows):
            for col in range(args.cols):
                plot_args = '-yr -2 2'
                fp.write("    args{}{} = '{}'\n".format(
                    row+1, col+1, plot_args))
        fp.write('\n')

        if args.marker:
            for row in range(args.rows):
                for col in range(args.cols):
                    fp.write("    axs[{}].text(-0.2, 1, '({})', transform=axs[{}].transAxes)\n".format(
                        col+args.cols*row, chr(97 + col+args.cols*row), col+args.cols*row))
            fp.write('\n')

        if args.title:
            for row in range(args.rows):
                for col in range(args.cols):
                    fp.write("    axs[{}].set_title('({})')\n".format(
                        col+args.cols*row, chr(97 + col+args.cols*row)))
            fp.write('\n')

        if args.tt:
            for col in range(args.cols):
                fp.write("    axs[{}].set_title('({})')\n".format(
                    col, chr(97 + col)))
            fp.write('\n')

        if args.xlabel:
            for row in range(args.rows):
                for col in range(args.cols):
                    fp.write("    axs[{}].set_xlabel('({})')\n".format(
                        col+args.cols*row, chr(97 + col+args.cols*row)))
            fp.write('\n')

        if args.nk:
            for row in range(args.rows):
                for col in range(args.cols):
                    fp.write("    plot{}{} = PlotBand(file=file{}{},symbols=symbols{}{}, args=get_args(args{}{}), fig=fig, ax=axs[{}])\n".format(
                        row+1, col+1, row+1, col+1, row+1, col+1, row+1, col+1, col+args.cols*row))
                    fp.write('    del plot{}{}\n'.format(row+1, col+1))
        else:
            for row in range(args.rows):
                for col in range(args.cols):
                    fp.write("    plot{}{} = PlotBand(file=file{}{}, divisions=kpoints{}{}.division,symbols=kpoints{}{}.symbols, args=get_args(args{}{}), fig=fig, ax=axs[{}])\n".format(
                        row+1, col+1, row+1, col+1, row+1, col+1, row+1, col+1, row+1, col+1, col+args.cols*row))
                    fp.write('    del plot{}{}\n'.format(row+1, col+1))
        fp.write('\n')

        if args.lb:
            for row in range(args.rows):
                for col in range(args.cols):
                    if row != args.rows-1:
                        fp.write("    axs[{}].set_xticklabels([])\n".format(
                            col+args.cols*row))
            fp.write('\n')
            for row in range(args.rows):
                for col in range(args.cols):
                    if row != args.rows-1:
                        fp.write("    axs[{}].set_xlabel(None)\n".format(
                            col+args.cols*row))
            fp.write('\n')
            for row in range(args.rows):
                for col in range(args.cols):
                    if col != 0:
                        fp.write("    axs[{}].set_yticklabels([])\n".format(
                            col+args.cols*row))
            fp.write('\n')
            for row in range(args.rows):
                for col in range(args.cols):
                    if col != 0:
                        fp.write("    axs[{}].set_ylabel(None)\n".format(
                            col+args.cols*row))
            fp.write('\n')

        fp.write("    plt.savefig('band.png')\n")
        fp.write('    plt.show()')


if __name__ == '__main__':
    path = Path('./plot_subfig.py')
    if path.exists():
        print("There is already a 'plot_subfig' file in the current directory, do you want to overwrite it?")
        a = input('yes or no?  ')
        if a == 'yes':
            args = get_args()
            generate_subfig_plot(args=args)
        else:
            pass
    else:
        args = get_args()
        generate_subfig_plot(args=args)
