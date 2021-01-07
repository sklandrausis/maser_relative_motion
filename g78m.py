import sys
import argparse

from matplotlib import rc, cm
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


def check_if_group_is_in_file(file, group):
    input_file = "groups/" + "/" + file
    group_nr = np.loadtxt(input_file, unpack=True, usecols=0)

    if group not in group_nr:
        return False
    else:
        return True


def main(group_numbers):
    rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=12)
    minorLocatorx = MultipleLocator(20)
    minorLocatory = MultipleLocator(20)
    minorLocatorvel = MultipleLocator(1)

    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    input_files = []

    for file in file_order:
        input_files.append(file)

    v_max = []
    v_min = []
    for index in range(0, len(input_files)):
        input_file = "groups/" + "/" + input_files[index].split(".")[0] + ".groups"
        group_tmp, channel_tmp, velocity_tmp, intensity_tmp, integral_intensity_tmp, ra_tmp, dec_tmp = np.loadtxt(
            input_file, unpack=True)

        v_max.append(max(velocity_tmp))
        v_min.append(min(velocity_tmp))

    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}

    bad_epoch_dict = {}
    for g in group_numbers:
        for file in input_files:
            if not check_if_group_is_in_file(file.split(".")[0] + ".groups", g):
                if file not in bad_epoch_dict.keys():
                    bad_epoch_dict[file] = 1
                else:
                    bad_epoch_dict[file] = +1

    for file in bad_epoch_dict.keys():
        if bad_epoch_dict[file] == len(group_numbers):
            input_files.remove(file)
            del dates[file.split(".")[0]]

    fig, ax = plt.subplots(nrows=2, ncols=len(input_files), figsize=(16, 16))

    data_dict = dict()
    for index in range(0, len(input_files)):
        input_file = "groups/" + "/" + input_files[index].split(".")[0] + ".groups"

        group_tmp, channel_tmp, velocity_tmp, intensity_tmp, integral_intensity_tmp, ra_tmp, dec_tmp = np.loadtxt(
            input_file, unpack=True)
        for j in group_numbers:
            velocity = np.empty(0)
            intensity = np.empty(0)
            ra = np.empty(0)
            dec = np.empty(0)
            if j not in data_dict.keys():
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
                data_dict[j] = [velocitys, vms, vxs, dvs, intensitys, ras, decs,
                                avgs_ra, avgs_dec, max_ra, min_ra, min_dec, max_dec]

            for i in range(0, len(channel_tmp)):
                if group_tmp[i] == int(j):
                    velocity = np.append(velocity, velocity_tmp[i])
                    intensity = np.append(intensity, intensity_tmp[i])
                    ra = np.append(ra, ra_tmp[i])
                    dec = np.append(dec, dec_tmp[i])

            if len(velocity) == 0:
                vm = 0
                vx = 0
            else:
                vm = min(velocity)
                vx = max(velocity)

            data_dict[j][0].append(velocity)
            data_dict[j][1].append(vm)
            data_dict[j][2].append(vx)
            data_dict[j][3].append(intensity)
            data_dict[j][4].append(ra)
            data_dict[j][5].append(dec)
            if len(ra) == 0:
                data_dict[j][6].append(0)
                data_dict[j][7].append(0)
                data_dict[j][8].append(0)
                data_dict[j][9].append(0)
                data_dict[j][10].append(0)
                data_dict[j][11].append(0)
            else:
                data_dict[j][6].append(np.mean(ra))
                data_dict[j][7].append(np.mean(dec))
                data_dict[j][8].append(abs(np.max(ra)))
                data_dict[j][9].append(abs(np.min(ra)))
                data_dict[j][10].append(abs(np.min(dec)))
                data_dict[j][11].append(abs(np.max(dec)))

    symbols = ["o", "*", "v", "^", "<", ">", "1", "2", "3", "4"]

    vm_min = []
    vx_max = []
    for j in group_numbers:
        vm_min.append(min(data_dict[j][1]))
        vx_max.append(max(data_dict[j][2]))

    for index in range(0, len(input_files)):
        vms = []
        vxs = []
        coord_ranges = []
        ra_maxs = []
        ra_mins = []
        dec_maxs = []
        dec_mins = []

        for j in group_numbers:
            symbol = symbols[group_numbers.index(j)]
            velocity = data_dict[j][0][index]
            vm = data_dict[j][1][index]
            vms.append(vm)
            vx = data_dict[j][2][index]
            vxs.append(vx)
            intensity = data_dict[j][3][index]
            ra = data_dict[j][4][index]
            dec = data_dict[j][5][index]
            coord_range = max(max(data_dict[j][8]) - min(data_dict[j][9]),
                              max(data_dict[j][11]) - min(data_dict[j][10]))
            coord_ranges.append(coord_range)
            if len(ra) == 0:
                ra_maxs.append(0)
                ra_mins.append(0)
                dec_maxs.append(0)
                dec_mins.append(0)
            else:
                ra_maxs.append(max(ra))
                ra_mins.append(min(ra))
                dec_maxs.append(max(dec))
                dec_mins.append(min(dec))
            title = input_files[index].split(".")[0].upper() + "-" + dates[input_files[index].split(".")[0]]
            ax[0][0].set_ylabel('Flux density (Jy)', fontsize=12)

            for i in range(len(velocity) - 1):
                if velocity[i] < min(v_min) or velocity[i] > max(v_max):
                    c = (0, 0, 0)
                else:
                    c = cm.jet((velocity[i] - min(v_min)) / (max(v_max) - min(v_min)), 1)

                ax[0][index].scatter((velocity[i], velocity[i + 1]), (intensity[i], intensity[i + 1]), color=c, lw=2,
                                     marker=symbol)
                ax[0][index].xaxis.set_minor_locator(minorLocatorvel)
                ax[0][index].set_title(title, size=12)
                ax[0][index].set_xlabel('$V_{\\rm LSR}$ (km s$^{-1}$)', fontsize=12)

                ax[1][0].set_ylabel('$\\Delta$ Dec (mas)', fontsize=12)
                ax[1][index].scatter(ra[i], dec[i], s=max(coord_ranges) * np.sqrt(intensity[i]), color=c, lw=2, marker=symbol)
                ax[1][index].set_aspect("equal", adjustable='box')
                ax[1][index].set_xlabel('$\\Delta$ RA (mas)', fontsize=12)
                ax[1][index].xaxis.set_minor_locator(minorLocatorx)
                ax[1][index].yaxis.set_minor_locator(minorLocatory)
                ax[1][index].invert_xaxis()

        ax[0][index].set_xlim(min(vms) - 0.5, max(vxs) + 0.5)
        coord_range_max = max(coord_ranges) + 125
        ax[1][index].set_xlim(np.mean((max(ra_maxs), min(ra_mins))) - (coord_range_max * 2) - 0.5, np.mean((max(ra_maxs), min(ra_mins))) + (coord_range_max * 2) + 0.5)
        ax[1][index].set_ylim(np.mean((max(dec_maxs), min(dec_mins))) - (coord_range_max * 2) - 0.5, np.mean((max(dec_maxs), min(dec_mins))) + (coord_range_max * 2) + 0.5)

    plt.tight_layout()
    plt.subplots_adjust(top=0.97, bottom=0, wspace=0.18, hspace=0, left=0.05, right=0.99)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot group')
    parser.add_argument('group_numbers', type=int, help='group numbers',  nargs='+',)
    args = parser.parse_args()
    main(args.group_numbers)
    sys.exit(0)