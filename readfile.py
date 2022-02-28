# -*- coding: utf-8 -*-
from timeit import timeit
import numpy as np
from lxml.etree import parse
import re
import time


class ReadVasprun(object):
    def __init__(self, file="vasprun.xml"):
        self.file = parse(file)

    @property
    def symbols(self):
        symbolstr = self.file.xpath(
            "/modeling/atominfo/array[@name='atoms']/set")[0].xpath('string(.)')
        symbols = symbolstr.split()[0:-1:2]
        return symbols

    @property
    def nedos(self):
        nedos = int(self.file.xpath(
            "/modeling/parameters/separator[@name='dos']/i[@name='NEDOS']")[0].text)
        return nedos

    @property
    def ionnum(self):
        ionnum = int(self.file.xpath(
            "/modeling/atominfo/atoms")[0].text)
        return ionnum

    @property
    def positions(self):
        positionstr = self.file.xpath(
            "/modeling/structure[@name='finalpos']/varray[@name='positions']")[0].xpath('string(.)')
        positionarray = np.array(
            positionstr.split(), dtype=float).reshape(-1, 3)
        return positionarray

    @property
    def real_cell(self):
        basis = self.file.xpath(
            "/modeling/structure[@name='finalpos']/crystal/varray[@name='basis']")[0].xpath('string(.)')
        cellarray = np.array(basis.split(), dtype=float).reshape(3, 3)
        return cellarray

    @property
    def rec_cell(self):
        basis = self.file.xpath(
            "/modeling/structure[@name='finalpos']/crystal/varray[@name='rec_basis']")[0].xpath('string(.)')
        cellarray = np.array(basis.split(), dtype=float).reshape(3, 3)
        return cellarray

    @property
    def fermi(self):
        fermi_energy = float(self.file.xpath(
            "/modeling/calculation/dos/i")[0].text)
        return fermi_energy

    @property
    def kpoints(self):
        kpoints = self.file.xpath(
            "/modeling/kpoints/varray[@name='kpointlist']")[0].xpath('string(.)')
        kpointsarray = np.array(kpoints.split(), dtype=float).reshape(-1, 3)
        return kpointsarray

    @property
    def weights(self):
        weights = self.file.xpath(
            "/modeling/kpoints/varray[@name='weights']")[0].xpath('string(.)')
        weightsarray = np.array(weights.split(), dtype=float)
        return weightsarray

    @property
    def nbands(self):
        nbands = int(self.file.xpath(
            "/modeling/parameters/separator[@name='electronic']/i[@name='NBANDS']")[0].text)
        return nbands

    @property
    def eigen(self):
        """a four-dimensional array about eigenvalues

        Returns
        -------
        <class 'numpy.ndarray'>
            The first dimension represents the different spins

            The second demension represents the different kpoints

            The third demension represents the different bands

            The fourth demension has two columns and they represent eigenenergy and occupied number respectively.
        """
        eigen = self.file.xpath(
            "/modeling/calculation/eigenvalues/array/set")[0].xpath('string(.)')
        eigenarray = np.array(eigen.split(
        ), dtype=float).reshape(-1, len(self.kpoints), self.nbands, 2)
        return eigenarray

    @property
    def total_dos(self):
        """a three-dimensional array about density of states

        Returns
        -------
        <class 'numpy.ndarray'>
            The first dimension represents the different spins

            The second demension represents the different points of energy

            The third demension has three columns and they represent energy, total and integrated respectively.

        """
        dos = self.file.xpath(
            "/modeling/calculation/dos/total/array/set")[0].xpath('string(.)')
        dosarray = np.array(
            dos.split(), dtype=float).reshape(-1, self.nedos, 3)
        return dosarray

    @property
    def partial_dos(self):
        """a four-dimensional array about partial dos

        Returns
        -------
        <class 'numpy.ndarray'>
            The first dimension represents the different ions

            The second demension represents the different spins

            The third demension represents the different points of energy, the totol is nedos

            The fourth demension represents the energy and different orbitals, the index 0 is energy and orbitals start index from 1
        """
        dos = self.file.xpath(
            "/modeling/calculation/dos/partial/array/set")[0].xpath('string(.)')
        orbitals = self.file.xpath(
            "/modeling/calculation/dos/partial/array/field")
        dosarray = np.array(dos.split(), dtype=float).reshape(
            self.ionnum, -1, self.nedos, len(orbitals))
        return dosarray

    @property
    def project(self):
        """a five-dimensional array about partial dos

        Returns
        -------
        <class 'numpy.ndarray'>
            The first dimension represents the different spins

            The second demension represents the different kpoints

            The third demension represents the different bands

            The fourth demension represents the different ions

            The fifth demension represents  different orbitals, the orbitals start index from 0
        """
        projectstring = self.file.xpath(
            "/modeling/calculation/projected/array/set")[0].xpath('string(.)')
        orbitals = self.file.xpath(
            "/modeling/calculation/projected/array/field")
        projectarray = np.array(projectstring.split(
        ), dtype=float).reshape(-1, len(self.kpoints), self.nbands, self.ionnum, len(orbitals))
        return(projectarray)


class ReadKpoints(object):
    def __init__(self, file="KPOINTS"):
        self.file = file
        self.lines = self.readlines()

    def readlines(self):
        try:
            with open(self.file) as fk:
                lines = fk.readlines()
        except:
            print(
                'open the Kponts file you specified failed , try to  read the "KPOINTS" in current directory')
            try:
                with open('KPOINTS') as fk:
                    lines = fk.readlines()
            except:
                print('failed to read the "KPOINTS" in current directory')
                lines = None
        return lines

    @property
    def division(self):
        try:
            return int(self.lines[1])
        except:
            return None

    @property
    def symbols(self):
        symbollist = []
        inputlist = []
        try:
            for line in self.lines:
                if r'!' in line:
                    symbol = re.findall(r' \! (.*)', line)[0]
                    if inputlist == []:
                        inputlist.append(symbol)
                    else:
                        if symbol != inputlist[-1]:
                            inputlist.append(symbol)
                    if '$' in symbol:
                        symbol = r'\rm {}'.format(symbol)
                    else:
                        symbol = r'$\rm {}$'.format(symbol)
                    if symbollist == []:
                        symbollist.append(symbol)
                    else:
                        if symbol != symbollist[-1]:
                            symbollist.append(symbol)
            print(
                'The high symmetry point used in the calculation is: {}'.format(inputlist))
        except:
            print('The high symmetry point used in the calculation is: []')
        return symbollist


if __name__ == "__main__":
    start = time.perf_counter()
    data = ReadVasprun('mnbt.xml')
    a = data.partial_dos
    print(a)
    end = time.perf_counter()
    print(end-start)
    # print(a)
    # kpoints = ReadKpoints()
    # symbols = kpoints.symbols
    # print(kpoints.division)
