import sys
import argparse

from matplotlib import rc, cm
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
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


def main(group_numbers):
    rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=12)
    minorLocatorx = MultipleLocator(20)
    minorLocatory = MultipleLocator(20)
    minorLocatorvel = MultipleLocator(1)

    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    input_files = []

    for file in file_order:
        input_files.append(file)

    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}

    fig, ax = plt.subplots(nrows=2, ncols=len(input_files), figsize=(16, 16))

    data_dict = dict()
    for index in range(0, len(input_files)):
        input_file = "groups/" + "/" + input_files[index].split(".")[0] + ".groups"
        velocity = np.empty(0)
        intensity = np.empty(0)
        ra = np.empty(0)
        dec = np.empty(0)
        group_tmp, channel_tmp, velocity_tmp, intensity_tmp, integral_intensity_tmp, ra_tmp, dec_tmp = np.loadtxt(
            input_file, unpack=True)
        for j in group_numbers:
            if j not in data_dict.keys():
                coord_ranges = []
                velocitys = []
                vms = []
                vxs = []
                dvs = []
                intensitys = []
                ras = []
                decs = []
                avgs_ra = []
                avgs_dec = []
                max_ra = []
                min_ra = []
                min_dec = []
                max_dec = []
                data_dict[j] = [coord_ranges, velocitys, vms, vxs, dvs, intensitys, ras, decs,
                                avgs_ra, avgs_dec, max_ra, min_ra, min_dec, max_dec]
            for i in range(0, len(channel_tmp)):
                if group_tmp[i] == int(j):
                    velocity = np.append(velocity, velocity_tmp[i])
                    intensity = np.append(intensity, intensity_tmp[i])
                    ra = np.append(ra, ra_tmp[i])
                    dec = np.append(dec, dec_tmp[i])

            dv = (max(velocity) - min(velocity))
            vm = min(velocity)
            vx = max(velocity)

            coord_range = max(abs(abs(max(ra)) - abs(min(ra))), abs(abs(max(dec)) - abs(min(dec))))
            data_dict[j][0].append(coord_range)
            data_dict[j][1].append(velocity)
            data_dict[j][2].append(vm)
            data_dict[j][3].append(vx)
            data_dict[j][4].append(dv)
            data_dict[j][5].append(intensity)
            data_dict[j][6].append(ra)
            data_dict[j][7].append(dec)

            data_dict[j][8].append(np.mean(ra))
            data_dict[j][9].append(np.mean(dec))
            data_dict[j][10].append(np.max(ra))
            data_dict[j][11].append(np.min(ra))
            data_dict[j][12].append(np.min(dec))
            data_dict[j][13].append(np.max(dec))

    symbols = ["o", "*", "v", "^", "<", ">", "1", "2", "3", "4"]
    for index in range(0, len(input_files)):
        for j in group_numbers:
            symbol = symbols[group_numbers.index(j)]
            velocity = data_dict[j][1][index]
            vm = data_dict[j][2][index]
            vx = data_dict[j][3][index]
            dv = data_dict[j][4][index]
            intensity = data_dict[j][5][index]
            dec = data_dict[j][6][index]
            ra = data_dict[j][7][index]
            title = input_files[index].split(".")[0].upper() + "-" + dates[input_files[index].split(".")[0]]
            ax[0][0].set_ylabel('Flux density (Jy)', fontsize=12)

            for i in range(len(velocity) - 1):
                if velocity[i] < vm or velocity[i] > vx:
                    c = (0, 0, 0)
                else:
                    c = cm.jet((velocity[i] - vm) / dv, 1)

                ax[0][index].scatter((velocity[i], velocity[i + 1]), (intensity[i], intensity[i + 1]), color=c, lw=2,
                                     marker=symbol)
                ax[0][index].set_xlim(min(velocity) - 0.5, max(velocity) + 0.5)
                ax[0][index].xaxis.set_minor_locator(minorLocatorvel)
                ax[0][index].set_title(title, size=12)
                ax[0][index].set_xlabel('$V_{\\rm LSR}$ (km s$^{-1}$)', fontsize=12)

                rel = []
                ax[1][0].set_ylabel('$\\Delta$ Dec (mas)', fontsize=12)
                #ax[1][index].scatter(ra, dec, marker=symbol, s=0.1 * np.sqrt(intensity[i]))

                for k in range(0, len(ra)):
                    el = Circle((ra[k], dec[k]),radius=0.1 * np.sqrt(intensity[i]), angle=0, lw=2, hatch=symbol)

                    ax[1][index].add_artist(el)
                    c = cm.jet((velocity[k] - vm) / dv, 1)
                    el.set_facecolor(c)
                    rel.append([ra[k], dec[k], velocity[k]])

                coord_range = max(max(data_dict[j][10]) - min(data_dict[j][11]), max(data_dict[j][12]) - min(data_dict[j][13]))
                #coord_range = max(max(max_ra) - min(min_ra), max(max_dec) - min(min_dec))
                ax[1][index].set_aspect("equal", adjustable='box')
                #ax[1][index].set_xlim(np.mean((max(data_dict[j][10]), min(data_dict[j][11]))) - (coord_range/2) - 0.5, np.mean((max(data_dict[j][10]), min(data_dict[j][11]))) + (coord_range/2) + 0.5)
                #ax[1][index].set_ylim(np.mean((max(data_dict[j][12]), min(data_dict[j][13]))) - (coord_range/2) - 0.5, np.mean((max(data_dict[j][12]), min(data_dict[j][13]))) + (coord_range/2) + 0.5)
                ax[1][index].set_xlim(min(ra) - 0.5, max(ra) + 0.5)
                ax[1][index].set_ylim(min(dec) - 0.5, max(dec) + 0.5)
                ax[1][index].set_xlabel('$\\Delta$ RA (mas)', fontsize=12)
                ax[1][index].xaxis.set_minor_locator(minorLocatorx)
                ax[1][index].yaxis.set_minor_locator(minorLocatory)
                ax[1][index].invert_xaxis()

    plt.tight_layout()
    plt.subplots_adjust(top=0.97, bottom=0, wspace=0.18, hspace=0, left=0.05, right=0.99)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot group')
    parser.add_argument('group_numbers', type=int, help='group numbers',  nargs='+',)
    args = parser.parse_args()
    main(args.group_numbers)
    sys.exit(0)
