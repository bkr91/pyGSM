# standard library imports
# local application imports
import sys
from utilities import manage_xyz
import os
from os import path
import subprocess
import re

# third party
import numpy as np

sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

try:
    from .base_lot import Lot
except:
    from base_lot import Lot


class DFTB(Lot):

    def __init__(self, options):
        super(DFTB, self).__init__(options)
        self.iii = 0
        os.system('rm -f dftb_jobs.txt')
        print(" making folder scratch/{}".format(self.node_id))
        os.system('mkdir -p scratch/{}'.format(self.node_id))
        os.system('cp {} scratch/{}'.format(self.lot_inp_file, self.node_id))
        if options['DFTB_lattice_file']:
            self.periodic = True
            self.DFTB_lattice_file = options['DFTB_lattice_file']
        else:
            self.periodic = False

    def run(self, geom, multiplicity, state):
        # print("In run():")
        # print("multiplicity=",multiplicity)
        # print("state=", state)
        # print("")
        owd = os.getcwd()
        manage_xyz.write_xyz(
            'scratch/{}/tmp.xyz'.format(self.node_id), geom, scale=1.0)
        os.system('xyz2gen {} scratch/{}/tmp.xyz'.format('-l ' +
                  self.DFTB_lattice_file if self.periodic else '', self.node_id))
        #os.system('./xyz2gen scratch/{}/tmp.xyz'.format(self.node_id))
        os.chdir('scratch/{}'.format(self.node_id))
        os.system('pwd')
        cmd = "dftb+"
        proc = subprocess.Popen(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                )
        stdout, stderr = proc.communicate()
        # print(stdout)
        # print(stderr)
        with open('dftb.out', 'a') as out:
            print("Stdout:", file=out)
            print(stdout.decode('utf-8'), file=out)
        with open('dftb.err', 'a') as out:
            print("Stderr:", file=out)
            print(stderr.decode('utf-8'), file=out)

        ofilepath = "detailed.out"
        with open(ofilepath, 'r') as ofile:
            olines = ofile.readlines()

        #self.Energies = {}
        temp = 0
        tmpgrada = []
        tmpgrad = []
        pattern = re.compile(r"Total energy:\s+[-+]?[0-9]*\.?[0-9]+ H")
        for line in olines:
            # print(line)
            for match in re.finditer(pattern, line):
                # print("match")
                tmpline = line.split()
                # self.E.append((1, 0, float(tmpline[2])))
                self._Energies[(multiplicity, state)] = self.Energy(
                    float(tmpline[2]), 'Hartree')
                self.Energies[(multiplicity, state)
                              ] = self._Energies[(multiplicity, state)]
            if line.strip() == "Total Forces":
                temp += 1
                # print(line)
            elif temp > 0:
                # print(line)
                tmpline = line.split()
                tmpgrad.append([float(i) for i in tmpline])
                temp += 1
            if temp > len(self.atoms):
                break
        # print("printing tmpgrad...")
        # print(tmpgrad)
        # print("")
        # print("printing modded tmpgrad...")
        # print([row[1:] for row in tmpgrad])
        # print("")
        # print("printing modded tmpgrad as array...")
        # arr = np.array([row[1:] for row in tmpgrad])
        # print(arr)
        # print("")
        # print(arr[0][0])
        # print("")

        tmpgrada.append(tmpgrad)

        self.grada = []
        for i, state in enumerate(self.states):
            if state[0] == 1:
                self.grada.append((1, state[1], tmpgrada[i]))
            if state[0] == 3:
                self.grada.append((3, state[1], tmpgrada[i]))
        self.hasRanForCurrentCoords = True
        # print(self.grada)

        # self.Gradients={}
        # for multiplicity,state,g in self.grada:
        #     self._Gradients[(multiplicity,state)]=self.Gradient(np.array(g),"Hartree/Bohr")
        # print(self.Gradients)
        # print("multiplicity=",multiplicity)
        # print("state=", state)
        # print("")
        # print("self._Energies")
        # print(self._Energies)

        # IMPORTANT: multiply forces from dftb by (-1) to get gradient as expected by program
        self._Gradients[(1, 0)] = self.Gradient(
            (-1)*np.array([row[1:] for row in tmpgrad]), "Hartree/Bohr")
        self.Gradients[(1, 0)] = self._Gradients[(1, 0)]

        self.iii += 1
        if self.node_id == 0 and self.iii == 1:
            print(self.Gradients)
            sys.exit()

        os.chdir(owd)
        # print(self.Energies)
        # print(self._Gradients)
        # print("DEBUG: Exiting...")
        # sys.exit(1)
        return


if __name__ == '__main__':
    filepath = "../../data/ethylene.xyz"
    dftb = DFTB.from_options(
        states=[(1, 0)], fnm=filepath, lot_inp_file='../../data/dftb_in.hsd')
    geom = manage_xyz.read_xyz(filepath)
    xyz = manage_xyz.xyz_to_np(geom)
    print(dftb.get_energy(xyz, 1, 0))
    print(dftb.get_gradient(xyz, 1, 0))

    xyz = xyz + np.random.rand(xyz.shape[0], xyz.shape[1])*0.1
    print(dftb.get_energy(xyz, 1, 0))
    print(dftb.get_gradient(xyz, 1, 0))
