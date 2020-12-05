import sys
import os
import argparse
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator
import numpy as np


def main(ispec_files_dir, input_files_dir):
    rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=12)
    minorLocatorx = MultipleLocator(20)
    minorLocatory = MultipleLocator(20)
    minorLocatorvel = MultipleLocator(1)

    ispec_files = os.listdir(ispec_files_dir)
    input_files = os.listdir(input_files_dir)

    file_pairs = []
    for i in range(0, len(ispec_files)):
        for j in range(0, len(input_files)):
            name_ispec_file = ispec_files[i].split("_")[0].lower()
            name_input_file = input_files[j].split(".")[0].lower()
            if name_ispec_file == name_input_file:
                file_pairs.append((ispec_files[i], input_files[j]))
    print(file_pairs)

    fig, ax = plt.subplots(nrows=2, ncols=5)

    for index in range(0, len(file_pairs)):
        title = file_pairs[index][0].split("_")[0].lower()
        ispec_file = ispec_files_dir + "/" + file_pairs[index][0]
        input_file = input_files_dir + "/" + file_pairs[index][1]

        ch, v10, i1, i2, ra, dec = np.loadtxt(input_file, unpack=True)
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
            else:
                c = cm.jet((v[i] - vm) / dv, 1)

            ax[0][index].plot((v[i], v[i + 1]), (s[i], s[i + 1]), c=c, lw=2)
            ax[0][index].set_xlim(-12, -2)
            ax[0][index].xaxis.set_minor_locator(minorLocatorvel)
            ax[0][index].set_title(title, size=12)
            ax[0][index].set_xlabel('$V_{\\rm LSR}$ (km s$^{-1}$)', fontsize=12)
            ax[0][index].set_ylabel('Flux density (Jy)', fontsize=12)

        rel = []
        for i in range(len(ra)):
            el = Circle((ra[i], dec[i]), radius=10 * np.sqrt(i1[i]), angle=0, lw=2)
            ax[1][index].add_artist(el)
            c = cm.jet((v1[i] - vm) / dv, 1)
            el.set_facecolor(c)
            rel.append([ra[i], dec[i], v1[i]])
        ax[1][index].add_artist(
            Circle((300, -200), radius=10, angle=0, edgecolor='black', facecolor='white', alpha=0.9))
        ax[1][index].annotate('1 Jy beam$^{-1}$', [290, -200], fontsize=12)
        ax[1][index].set_aspect("equal")
        ax[1][index].set_xlim(ra.min() - 100, ra.max() + 100)
        ax[1][index].set_ylim(dec.min() - 100, dec.max() + 100)
        ax[1][index].set_xlabel('$\\Delta$ RA (mas)', fontsize=12)
        ax[1][index].set_ylabel('$\\Delta$ Dec (mas)', fontsize=12)
        ax[1][index].xaxis.set_minor_locator(minorLocatorx)
        ax[1][index].yaxis.set_minor_locator(minorLocatory)

    plt.subplots_adjust(top=0.95)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot ispec')
    parser.add_argument('ispec_dir', type=str, help='ispec files directory.')
    parser.add_argument('input_file_dir', type=str, help='Input files directory')
    args = parser.parse_args()
    main(args.ispec_dir, args.input_file_dir)
    sys.exit(0)
