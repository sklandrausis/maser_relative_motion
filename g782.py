import sys
import os
import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator
import numpy as np


def main(ispec_files_dir, input_files_dir):
    minorLocatorvel = MultipleLocator(1)
    ispec_files = os.listdir(ispec_files_dir)
    input_files = os.listdir(input_files_dir)

    file_pairs = []
    for i in range(0, len(ispec_files)):
        for j in range(0, len(input_files)):
            name_ispec_file = ispec_files[i].split("_")[0].lower()
            name_input_file = input_files[j].split("_")[0].lower()
            if name_ispec_file == name_input_file:
                file_pairs.append((ispec_files[i], input_files[j]))
    print(file_pairs)

    fig, ax = plt.subplots(nrows=1, ncols=5)
    for index in range(0, len(file_pairs)):
        title = file_pairs[index][0].split("_")[0].lower()
        ispec_file = ispec_files_dir + "/" + file_pairs[index][0]
        input_file = input_files_dir + "/" + file_pairs[index][1]

        ch, v10, i1, i2, x1, x2, x3, y1, y2, y3 = np.loadtxt(input_file, unpack=True)
        v1 = v10 / 1000.
        dv = (v1.max() - v1.min())
        vm = v1.min()
        vx = v1.max()

        nu, v, s = np.loadtxt(ispec_file, unpack=True)
        v = v / 1000.0
        sm = s.max()

        for i in range(len(v) - 1):
            if v[i] < vm or v[i] > vx:
                c = (0, 0, 0)
            #     print v[i]
            else:
                c = cm.jet((v[i] - vm) / dv, 1)

            ax[index].plot((v[i], v[i + 1]), (s[i], s[i + 1]), c=c, lw=2)
            ax[index].set_xlim(-12, -2)
            ax[index].xaxis.set_minor_locator(minorLocatorvel)
            ax[index].set_title(title, size=12)
            ax[index].set_xlabel('$V_{\\rm LSR}$ (km s$^{-1}$)', fontsize=12)
            ax[index].set_ylabel('Flux density (Jy)', fontsize=12)

    plt.subplots_adjust(top=0.95)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot ispec')
    parser.add_argument('ispec_dir', type=str, help='ispec files directory.')
    parser.add_argument('input_file_dir', type=str, help='Input files directory')
    args = parser.parse_args()
    main(args.ispec_dir, args.input_file_dir)
    sys.exit(0)
