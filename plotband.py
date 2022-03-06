# -*- coding: utf-8 -*-
from matplotlib import cm
import matplotlib as mpl
from matplotlib.colors import ListedColormap
from readfile import ReadVasprun
from readfile import ReadKpoints
import numpy as np
import matplotlib.pyplot as plt
from getarg import get_args
from matplotlib.collections import LineCollection
from pathlib import Path
from matplotlib.ticker import AutoMinorLocator, MultipleLocator


class PlotBand(ReadVasprun):
    def __init__(self, file='vasprun.xml', args=None, ax=None, symbols=[], divisions=None, fig=None):
        super().__init__(file)
        self.linecollection = None
        self.args = args
        self.ksymbols = symbols
        self.divisions = divisions
        if self.args.hse:
            self.use_index = np.where(self.weights == 0)[0]
        else:
            self.use_index = np.where(self.weights >= 0)[0]
        if self.args.pband:
            self.plot_proband(ax, fig)
        else:
            self.plot_band(ax, fig)

    @property
    def fermi(self):
        if self.args.fermi:
            try:
                return np.float64(self.args.fermi)
            except:
                return ReadVasprun(self.args.fermi).fermi
        else:
            fermi_energy = float(self.file.xpath(
                "/modeling/calculation/dos/i")[0].text)
            return fermi_energy

    def cal_length(self, vector1, vector2):
        result = 0
        if len(vector1) == len(vector2):
            for i in range(len(vector1)):
                result = result + (vector1[i]-vector2[i]) ** 2
            return result ** 0.5
        else:
            return 0

    def rectoreal(self, reclist, rec_cell):
        reallist = reclist.copy()
        for i in range(len(reclist)):
            reallist[i] = np.dot(reclist[i], rec_cell)
        return reallist

    def get_x(self, kpointreal):
        x = np.empty(shape=[len(kpointreal)], dtype=float)
        for i in range(len(kpointreal)):
            if i == 0:
                x[i] = 0
            else:
                x[i] = self.cal_length(kpointreal[i], kpointreal[i-1]) + x[i-1]
        return x

    def write_band_dat(self, x, y, nbands, file):
        with open(file, "w") as f:
            for num in range(nbands):
                for i, j in enumerate(x):
                    f.write("%f %f\n" % (j, y[i][num]))
                f.write("\n")

    def get_prolist(self):
        if self.args.atoms == None and self.args.orbitals == None:
            return None
        else:
            total_orbitals = ['s', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz',
                              'dx2-y2', 'fy3x2', 'fxyz', 'fyz2', 'fz3', 'fxz2', 'fzx2', 'fx3']
            prolist = [[], []]
            try:
                for proatom in self.args.atoms:
                    try:
                        prolist[0].append(int(proatom))
                    except:
                        for index, symbol in enumerate(self.symbols):
                            if symbol == proatom:
                                prolist[0].append(index)
            except TypeError:
                print(
                    'Did you forget choose atoms? I will plot all atoms for orbitals which you choose')
                prolist[0] = list(range(self.project.shape[-2]))

            try:
                for proorbital in self.args.orbitals:
                    try:
                        prolist[1].append(int(proorbital))
                    except:
                        for index, orbital in enumerate(total_orbitals):
                            if proorbital in orbital:
                                prolist[1].append(index)
            except TypeError:
                print(
                    'Did you forget choose orbitals? I will plot all orbitals for atoms which you choose')
                prolist[1] = list(range(self.project.shape[-1]))
                print(prolist)

            print('the projected atoms are {} and projected orbitals are {}'.format(
                [self.symbols[i] + ' ' + str(i) for i in prolist[0]], [total_orbitals[i] for i in prolist[1]]))
            return(prolist)

    def get_gap(self):
        withoutfermi = self.eigen[0, self.use_index, :, 0] - self.fermi
        extremum = np.array([[withoutfermi[:, i].min(), withoutfermi[:, i].max()]
                            for i in range(withoutfermi.shape[1])])
        positive_extremum = (extremum + np.absolute(extremum))/2
        negative_extremum = (extremum - np.absolute(extremum))/2
        if ((extremum[:, 0] * extremum[:, 1]) > 0).all():
            vbm_band = negative_extremum[:, 1].nonzero()[0][-1]
            cbm_band = positive_extremum[:, 0].nonzero()[0][0]
        else:
            vbm_band = negative_extremum[:, 1].nonzero()[0][-1]+1
            cbm_band = positive_extremum[:, 0].nonzero()[0][0] - 1
        vbm = withoutfermi[:, vbm_band].max()
        cbm = withoutfermi[:, cbm_band].min()
        vbm_index = np.array(np.where(withoutfermi[:, vbm_band] == vbm))
        vbm_kpoints = list(set(tuple(tuple(kpoint)
                           for kpoint in self.kpoints[self.use_index][vbm_index[0]])))
        cbm_index = np.array(np.where(withoutfermi[:, cbm_band] == cbm))
        cbm_kpoints = list(set(tuple(tuple(kpoint)
                           for kpoint in self.kpoints[self.use_index][cbm_index[0]])))

        print('========== 1th spin========')
        print('vbm locates at {} of {}th band'.format(vbm_kpoints, vbm_band + 1))
        print('cbm locates at {} of {}th band'.format(cbm_kpoints, cbm_band + 1))
        if ((extremum[:, 0] * extremum[:, 1]) > 0).all():
            print('gap is {}'.format(cbm - vbm))
        else:
            print('Metal')

        try:
            withoutfermi = self.eigen[1, self.use_index, :, 0] - self.fermi
            extremum = np.array([[withoutfermi[:, i].min(), withoutfermi[:, i].max()]
                                for i in range(withoutfermi.shape[1])])
            positive_extremum = (extremum + np.absolute(extremum))/2
            negative_extremum = (extremum - np.absolute(extremum))/2
            if ((extremum[:, 0] * extremum[:, 1]) > 0).all():
                vbm_band = negative_extremum[:, 1].nonzero()[0][-1]
                cbm_band = positive_extremum[:, 0].nonzero()[0][0]
            else:
                vbm_band = negative_extremum[:, 1].nonzero()[0][-1]+1
                cbm_band = positive_extremum[:, 0].nonzero()[0][0]-1
            vbm = withoutfermi[:, vbm_band].max()
            cbm = withoutfermi[:, cbm_band].min()
            vbm_index = np.array(np.where(withoutfermi[:, vbm_band] == vbm))
            vbm_kpoints = list(set(tuple(tuple(kpoint)
                               for kpoint in self.kpoints[self.use_index][vbm_index[0]])))
            cbm_index = np.array(np.where(withoutfermi[:, cbm_band] == cbm))
            cbm_kpoints = list(set(tuple(tuple(kpoint)
                               for kpoint in self.kpoints[self.use_index][cbm_index[0]])))

            print('========== 1th spin========')
            print('vbm locates at {} of {}th band'.format(
                vbm_kpoints, vbm_band + 1))
            print('cbm locates at {} of {}th band'.format(
                cbm_kpoints, cbm_band + 1))
            if ((extremum[:, 0] * extremum[:, 1]) > 0).all():
                print('gap is {}'.format(cbm - vbm))
            else:
                print('Metal')
        except:
            pass

    def plot_band(self, ax=None, fig=None):
        kpointlist = self.kpoints[self.use_index]
        kpointreal = self.rectoreal(kpointlist, self.rec_cell)
        xlist = self.get_x(kpointreal)
        ylist = self.eigen[:, self.use_index, :, :][:, :, :, 0] - self.fermi
        if ax is None:
            fig, ax = plt.subplots()
        ax.set_ylabel('Energy(eV)')
        ax.set_ylim(self.args.yr)
        ax.set_xlim(xlist[0], xlist[-1])
        for i in range(self.nbands):
            try:
                ax.plot(xlist, ylist[0, :, i],
                        color=self.args.color[0], zorder=1)
            except:
                ax.plot(xlist, ylist[0, :, i], color='red', zorder=1)
        try:
            for i in range(self.nbands):
                try:
                    ax.plot(xlist, ylist[1, :, i],
                            color=self.args.color[1], zorder=1)
                except:
                    ax.plot(xlist, ylist[1, :, i], color='blue', zorder=1)
        except:
            pass
        self.get_gap()
        ax.axhline(y=0, linestyle='--', color='black',
                   linewidth=mpl.rcParams['ytick.major.width'], zorder=0)

        if self.divisions is None or self.divisions == len(self.kpoints):
            xticks = [0]
            for i, x in enumerate(xlist):
                if x == xlist[i-1]:
                    xticks.append(x)
            xticks.append(xlist[-1])
        else:
            xtickindex = [0]+[(x+1)*self.divisions -
                              1 for x in range(len(xlist)//self.divisions)]
            xticks = [xlist[i] for i in xtickindex]
        for xtick in xticks:
            ax.axvline(x=xtick, linestyle='-', linewidth=mpl.rcParams['xtick.major.width'],
                       color='lightgrey', zorder=0)
        ax.set_xticks(xticks)
        if self.ksymbols == []:
            kpointlabellist = ['M', 'K', '$\Gamma$', 'K', '$\Gamma$',
                               '$\Gamma$', '$\Gamma$', '$\Gamma$', '$\Gamma$']
            print(
                'There is no label read from kpoints, use fake labels to plot fig, check it!!! ')
            ax.set_xticklabels(kpointlabellist[0:len(xticks)])
        else:
            try:
                ax.set_xticklabels(self.ksymbols)
            except ValueError:
                print('The number of kpoints label read from KPOINTS are different with the number of high symmetry kpoints, please check it!!!')
                ax.set_xticklabels(
                    (self.ksymbols+self.ksymbols)[0:len(xticks)])
        if self.args.xr:
            ax.set_xlim(xticks[self.args.xr[0]], xticks[self.args.xr[1]])

        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(bottom=False)
        return xlist, ylist

    def plot_proband(self, ax=None, fig=None):
        if self.args.atoms == ['all']:
            self.args.atoms = list(set(self.symbols))
        if self.args.orbitals == ['all']:
            self.args.orbitals = np.arange(self.project.shape[-1])
            print(self.args.atoms)

        prolist = self.get_prolist()
        kpointlist = self.kpoints[self.use_index]
        kpointreal = self.rectoreal(kpointlist, self.rec_cell)
        xlist = self.get_x(kpointreal)
        ylist = self.eigen[0, self.use_index, :, :] - self.fermi
        xlist2 = np.tile(xlist.repeat(2)[1:-1], self.nbands)
        ylist2 = ylist.repeat(2, axis=0)[1:-1, :].reshape(-1, order="F")

        if self.args.s:
            totalprolist = self.project[3,
                                        self.use_index, :, :, :][1:, :, :, :]
            proarray = (
                totalprolist[:, :, prolist[0], :][:, :, :, prolist[1]]
                .sum(axis=3)
                .sum(axis=2)
                .reshape(-1, order="F")
            )
            norm = mpl.colors.Normalize(proarray.min(), proarray.max())
        else:
            totalprolist = self.project[0,
                                        self.use_index, :, :, :][1:, :, :, :]
            proarray = (
                totalprolist[:, :, prolist[0], :][:, :, :, prolist[1]]
                .sum(axis=3)
                .sum(axis=2)
                .reshape(-1, order="F")
            )
            norm = mpl.colors.Normalize(0, proarray.max())
        if ax is None:
            fig, ax = plt.subplots()
        ax.set_ylabel("Energy(ev)")
        ax.set_yticks([-2, -1, 0, 1, 2])
        ax.set_ylim(self.args.yr)
        ax.set_xlim(xlist[0], xlist[-1])

        viridis_big = cm.get_cmap(self.args.color[0])
        newcmp = ListedColormap(viridis_big(np.linspace(0, 1, 128)))
        lc = LineCollection(
            np.array(list(zip(xlist2, ylist2))).reshape(-1, 2, 2), cmap=newcmp, norm=norm
        )

        ax.add_collection(lc)
        lc.set_array(proarray)
        if self.args.colorbar:
            if self.args.s:
                cbar = fig.colorbar(lc, ax=ax, ticks=np.linspace(
                    np.min(proarray), np.max(proarray), 5))
                cbar.ax.set_yticklabels(np.around(
                    np.linspace(np.min(proarray), np.max(proarray), 5), 2))
            else:
                cbar = fig.colorbar(lc, ax=ax, ticks=np.linspace(
                    0, np.max(proarray), 5))
                cbar.ax.set_yticklabels(np.around(
                    np.linspace(0, np.max(proarray), 5), 2))
        ax.axhline(y=0, linestyle="--", color="black",
                   linewidth=mpl.rcParams['ytick.major.width'], zorder=0)

        if self.divisions is None or self.divisions == len(self.kpoints):
            xticks = [0]
            for i, x in enumerate(xlist):
                if x == xlist[i-1]:
                    xticks.append(x)
            xticks.append(xlist[-1])
        else:
            xtickindex = [0]+[(x+1)*self.divisions -
                              1 for x in range(len(xlist)//self.divisions)]
            xticks = [xlist[i] for i in xtickindex]
        for xtick in xticks:
            ax.axvline(x=xtick, linestyle='-', linewidth=mpl.rcParams['xtick.major.width'],
                       color='lightgrey', zorder=0)
        ax.set_xticks(xticks)
        if self.ksymbols == []:
            kpointlabellist = ['M', 'K', '$\Gamma$', 'K', '$\Gamma$',
                               '$\Gamma$', '$\Gamma$', '$\Gamma$', '$\Gamma$']
            print(
                'There is no label read from kpoints, use fake labels to plot fig, check it!!! ')
            ax.set_xticklabels(kpointlabellist[0:len(xticks)])
        else:
            try:
                ax.set_xticklabels(self.ksymbols)
            except ValueError:
                print('The number of kpoints label read from KPOINTS are different with the number of high symmetry kpoints, please check it!!!')
                ax.set_xticklabels(
                    (self.ksymbols+self.ksymbols)[0:len(xticks)])
        if self.args.xr:
            ax.set_xlim(xticks[self.args.xr[0]], xticks[self.args.xr[1]])
        ax.tick_params(bottom=False)
        self.linecollection = lc
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(bottom=False)
        return lc


if __name__ == '__main__':
    cur_file_path = Path(__file__).resolve().parent
    try:
        mpl.rc_file('{}/matplotlibrc'.format(cur_file_path))
    except:
        pass
    fig, ax = plt.subplots()
    args = get_args()
    kpoints = ReadKpoints('KPOINTS')
    plot = PlotBand(file='vasprun.xml', divisions=kpoints.division,
                    symbols=kpoints.symbols, args=args, fig=fig, ax=ax)
    plt.savefig('band.png')
    plt.show()
