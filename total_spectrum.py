import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator
import numpy as np

from parsers.configparser_ import ConfigParser


def get_configs(section, key):
    """

    :param section: configuration file secti
    :param key: configuration file sections
    :return: configuration file section key
    """
    config_file_path = "config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def main():
    rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=12)
    minor_locator_x = MultipleLocator(20)
    minor_locator_y = MultipleLocator(20)
    minor_locator_vel = MultipleLocator(1)

    ispec_files_dir = get_configs("paths", "ispec")
    data_files_dir = get_configs("paths", "dataFiles")

    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    ispec_files = []
    data_files = []
    file_pairs = []
    for file in file_order:
        ispec_files.append(file.split(".")[0].upper() + "_FRING.ISPEC")
        data_files.append(file)
        file_pairs.append((file.split(".")[0].upper() + "_FRING.ISPEC", file))

    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}

    fig, ax = plt.subplots(nrows=2, ncols=len(data_files), figsize=(16, 16), gridspec_kw={'height_ratios': [2, 2]})

    max_ra = []
    min_ra = []
    min_dec = []
    max_dec = []
    v_s = []
    vms = []
    vxs = []
    dvs = []
    ss = []
    ras = []
    decs = []
    v1s = []
    i1s = []
    for index in range(0, len(file_pairs)):
        ispec_file = ispec_files_dir + "/" + file_pairs[index][0]
        data_file = data_files_dir + "/" + file_pairs[index][1]

        ch, v10, i1, i2, ra, dec = np.loadtxt(data_file, unpack=True)
        v1 = v10 / 1000.
        dv = (v1.max() - v1.min())
        vm = v1.min()
        vx = v1.max()

        nu, v, s = np.loadtxt(ispec_file, unpack=True)
        v = v / 1000.0
        v_s.append(v)
        vms.append(vm)
        vxs.append(vx)
        dvs.append(dv)
        ss.append(s)
        ras.append(ra)
        decs.append(dec)
        v1s.append(v1)
        i1s.append(i1)
        max_ra.append(np.max(ra))
        min_ra.append(np.min(ra))
        min_dec.append(np.min(dec))
        max_dec.append(np.max(dec))

    for index in range(0, len(file_pairs)):
        v = v_s[index]
        vm = vms[index]
        vx = vxs[index]
        dv = dvs[index]
        s = ss[index]
        ra = ras[index]
        dec = decs[index]
        v1 = v1s[index]
        i1 = i1s[index]
        title = file_pairs[index][0].split("_")[0].upper() + "-" + dates[file_pairs[index][1].split(".")[0]]
        ax[0][0].set_ylabel('Flux density (Jy)', fontsize=12)
        for i in range(len(v) - 1):
            if v[i] < vm or v[i] > vx:
                c = (0, 0, 0)
            else:
                c = cm.jet((v[i] - vm) / dv, 1)

            ax[0][index].plot((v[i], v[i + 1]), (s[i], s[i + 1]), c=c, lw=2)
            ax[0][index].set_xlim(-12, -2)
            ax[0][index].xaxis.set_minor_locator(minor_locator_vel)
            ax[0][index].set_title(title, size=12)
            ax[0][index].set_xlabel('$V_{\\rm LSR}$ (km s$^{-1}$)', fontsize=12)

        rel = []
        ax[1][0].set_ylabel('$\\Delta$ Dec (mas)', fontsize=12)
        for i in range(len(ra)):
            el = Circle((ra[i], dec[i]), radius=10 * np.sqrt(i1[i]), angle=0, lw=2)
            ax[1][index].add_artist(el)
            c = cm.jet((v1[i] - vm) / dv, 1)
            el.set_facecolor(c)
            rel.append([ra[i], dec[i], v1[i]])
        ax[1][index].add_artist(
            Circle((285, -200), radius=10, angle=0, edgecolor='black', facecolor='white', alpha=0.9))
        ax[1][index].annotate('1 Jy beam$^{-1}$', [275, -200], fontsize=12)
        ax[1][index].set_aspect("equal", adjustable='box')
        coord_range = max(max(max_ra) - min(min_ra), max(max_dec) - min(min_dec))
        ax[1][index].set_xlim(np.mean((max(max_ra), min(min_ra))) - (coord_range / 2) - 10*np.sqrt(max(i1)),
                              np.mean((max(max_ra), min(min_ra))) + (coord_range / 2) + 10*np.sqrt(max(i1)))
        ax[1][index].set_ylim(np.mean((max(max_dec), min(min_dec))) - (coord_range / 2) - 10*np.sqrt(max(i1)),
                              np.mean((max(max_dec), min(min_dec))) + (coord_range / 2) + 10*np.sqrt(max(i1)))
        ax[1][index].set_xlabel('$\\Delta$ RA (mas)', fontsize=12)
        ax[1][index].xaxis.set_minor_locator(minor_locator_x)
        ax[1][index].yaxis.set_minor_locator(minor_locator_y)
        ax[1][index].invert_xaxis()

    plt.tight_layout()
    plt.subplots_adjust(top=0.97, bottom=0, wspace=0.18, hspace=0, left=0.05, right=0.99)
    plt.show()


if __name__ == "__main__":
    main()
    sys.exit(0)
