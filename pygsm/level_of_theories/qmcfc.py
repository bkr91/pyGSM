# standard library imports
# local application imports
import sys
from utilities import manage_xyz
import os
from os import path
import subprocess
import re
import shutil

# third party
import numpy as np

sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))


try:
    from .base_lot import Lot
except:
    from base_lot import Lot


# for this class to work, the path to the qmcfc executable `self.qmcfc_exec` needs to be set appropriately.
# self.lot_inp_file corresponds to the normal qmcfc run-xx.in / md-xx.in file
# additionally a correct rst file for the calculated system (i.e. correct moldescriptor data) needs to be provided, however the geometry/forces/velocities can be arbitrary.


class QMCFC(Lot):

    def __init__(self, options):
        super(QMCFC, self).__init__(options)
        self.qmcfc_exec = "/home/bkr/qmcfc_dissociative_dev_20201111/qmcfc"
        print(f" making folder scratch/{self.node_id}")
        os.system(f'mkdir -p scratch/{self.node_id}')
        os.system(f'cp {self.lot_inp_file} scratch/{self.node_id}')
        os.system(
            f'cp {" ".join(self.important_files())} scratch/{self.node_id}')
        self.write_traj_list(f"scratch/{self.node_id}")
        self.info_file, self.energy_file, self.force_file = [self.get_filename(
            name) for name in ("info_file", "energy_file", "force_file")]
        if self.info_file == "" or self.energy_file == "" or self.force_file == "":
            raise Exception(
                f"The QMCFC input file {self.lot_inp_file} needs to specify info_file, energy_file and force_file")

    def important_files(self):
        fixed_names = ["intra_nonbonded.dat", "guff.dat",
                       "moldescriptor.dat", "dftb_energy.template", "dftb_force.template", "dftb_in.template"]

        return [self.get_filename(name) for name in ("start_file", "parameter_file", "topology_file")] + fixed_names + [self.lot_inp_file]

    def run(self, geom, multiplicity, state):
        owd = os.getcwd()
        os.chdir(f'scratch/{self.node_id}')
        manage_xyz.write_xyz(
            'geom.xyz', geom, scale=1.0)
        os.system('pwd')

        proc = subprocess.Popen([self.qmcfc_exec, self.lot_inp_file],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                )
        stdout, stderr = proc.communicate()

        with open('result.out', 'a') as out:
            print("Stdout:", file=out)
            print(stdout.decode('utf-8'), file=out)
        with open('result.err', 'a') as out:
            print("Stderr:", file=out)
            print(stderr.decode('utf-8'), file=out)

        tmpgrada = []
        tmpgrad = []

        print(os.listdir())
        qm_en_index, mm_en_index = self.get_energy_file_indices()
        with open(self.energy_file, "r") as f:
            s = f.readline().split()
            qm_energy = float(s[qm_en_index])
            mm_energy = float(s[mm_en_index])

        # convert energy from qmcfc (kcal/mol) to ORCA (Eh)
        total_energy = (qm_energy + mm_energy) * 0.0015936010974213599

        self._Energies[(multiplicity, state)] = self.Energy(
            total_energy, 'Hartree')
        self.Energies[(multiplicity, state)
                      ] = self._Energies[(multiplicity, state)]

        with open(self.force_file, "r") as f:
            # skip first 2 lines: n-atoms line + comment line
            for line in f.readlines()[2:]:
                if line.lower().split()[0].strip() == "x":  # x-particle
                    continue
                tmpgrad.append([float(i) for i in line.split()[1:]])

        tmpgrada.append(tmpgrad)

        self.grada = []
        for i, state in enumerate(self.states):
            if state[0] == 1:
                self.grada.append((1, state[1], tmpgrada[i]))
            if state[0] == 3:
                self.grada.append((3, state[1], tmpgrada[i]))
        self.hasRanForCurrentCoords = True

        # convert forces from qmcfc( kcal/mol/A) to ORCA(Eh/a0), also forces => gradient
        # *(-8.4329744e-4)
        self._Gradients[(1, 0)] = self.Gradient(
            (-8.4329744e-4)*np.array([row for row in tmpgrad]), "Hartree/Bohr")
        self.Gradients[(1, 0)] = self._Gradients[(1, 0)]

        for ofile in ["result_qmcfc.xyz", "result_qmcfc.out", "result_qmcfc.info", "result_qmcfc.forces", "result_qmcfc.en", "result_qmcfc.chrg", "result_qmcfc.vel"]:
            if path.exists(ofile):
                # name, ext = path.splitext(ofile)
                # shutil.move(ofile, f"{name}_{self.iii}{ext}")
                os.remove(ofile)

        os.chdir(owd)
        return

    def write_traj_list(self, path="."):
        with open(f"{path}/traj_list.dat", "w") as f:
            print("geom.xyz", file=f)

    def get_filename(self, query):
        with open(self.lot_inp_file, "r") as f:
            for line in f.readlines():
                if line.strip().startswith(query):
                    return line.split("=")[1].strip(" \t;\n")
        return ""

    def get_energy_file_indices(self):  # 0-based
        qm_en_index, mm_en_index = -1, -1
        with open(self.info_file, "r") as f:
            for i, line in enumerate(f.readlines()[3:-1]):
                if line.startswith("--") or line.strip() == '':
                    continue
                parts = line.strip("\n     |").split()
                if parts[0].strip() == "E(QM)":
                    qm_en_index = 2*i
                elif parts[0].strip() == "E(MM)":
                    mm_en_index = 2*i
                if len(parts) > 2 and parts[2].strip() == "E(QM)":
                    qm_en_index = 2*i+1
                elif len(parts) > 2 and parts[2].strip() == "E(MM)":
                    mm_en_index = 2*i+1

        if qm_en_index == -1 or mm_en_index == -1:
            raise Exception("Could not parse info file!")
        return qm_en_index, mm_en_index


if __name__ == '__main__':
    filepath = "../../data/ethylene.xyz"
    qmcfc = QMCFC.from_options(
        states=[(1, 0)], fnm=filepath, lot_inp_file='../../data/dftb_in.hsd')
    geom = manage_xyz.read_xyz(filepath)
    xyz = manage_xyz.xyz_to_np(geom)
    print(qmcfc.get_energy(xyz, 1, 0))
    print(qmcfc.get_gradient(xyz, 1, 0))

    xyz = xyz + np.random.rand(xyz.shape[0], xyz.shape[1])*0.1
    print(qmcfc.get_energy(xyz, 1, 0))
    print(qmcfc.get_gradient(xyz, 1, 0))
