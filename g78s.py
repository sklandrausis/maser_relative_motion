import sys
import argparse

from matplotlib import cm
from matplotlib import rcParams
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from parsers.configparser_ import ConfigParser


def gauss(x, *p):
    a, b, c = p
    return a*np.exp(-(x-b)**2*np.log(2)/(c**2))


def gauss2(x, *p):
    a1, b1, c1, a2, b2, c2 = p
    return a1*np.exp(-(x-b1)**2*np.log(2)/c1**2) + \
           a2*np.exp(-(x-b2)**2*np.log(2)/c2**2)


def get_configs(section, key):
    """
    :param section: configuration file secti
    :param key: configuration file sections
    :return: configuration file section key
    """
    config_file_path = "config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def get_configs_items():
    """
    :return: None
    """
    config_file_path = "config/plot.cfg"
    config = ConfigParser(config_file_path)
    return config.get_items("main")


def check_if_group_is_in_file(file, group):
    input_file = "groups/" + "/" + file
    group_nr = np.loadtxt(input_file, unpack=True, usecols=0)

    if group not in group_nr:
        return False
    else:
        return True


def main(group_number):
    configuration_items = get_configs_items()
    for key, value in configuration_items.items():
        rcParams[key] = value

    minorLocatorx = MultipleLocator(20)
    minorLocatory = MultipleLocator(20)
    minorLocatorvel = MultipleLocator(1)

    gauss2_list = get_configs("parameters", "gauss").split(";")
    gauss2_dict = dict()

    for epoch in gauss2_list:
        gauss2_dict[epoch.split(":")[0]] = epoch.split(":")[1].split(",")

    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    input_files = []

    for file in file_order:
        input_files.append(file)

    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}

    v_max = []
    v_min = []
    for index in range(0, len(input_files)):
        input_file = "groups/" + "/" + input_files[index].split(".")[0] + ".groups"
        group_tmp, channel_tmp, velocity_tmp, intensity_tmp, integral_intensity_tmp, ra_tmp, dec_tmp = np.loadtxt(
            input_file, unpack=True)

        v_max.append(max(velocity_tmp))
        v_min.append(min(velocity_tmp))

    bad_files = []
    for file in input_files:
        if not check_if_group_is_in_file(file.split(".")[0] + ".groups", group_number):
            bad_files.append(file)
            del dates[file.split(".")[0]]
    input_files = [file for file in input_files if file not in bad_files]

    fig, ax = plt.subplots(nrows=2, ncols=len(input_files), figsize=(16, 16), dpi=90)

    velocitys = []
    intensitys = []
    ras = []
    decs = []
    max_ra = []
    min_ra = []
    min_dec = []
    max_dec = []
    intensitys_max = []
    intensitys_min = []
    for index in range(0, len(input_files)):
        input_file = "groups/" + "/" + input_files[index].split(".")[0] + ".groups"
        velocity = np.empty(0)
        intensity = np.empty(0)
        ra = np.empty(0)
        dec = np.empty(0)
        group_tmp, channel_tmp, velocity_tmp, intensity_tmp, integral_intensity_tmp, ra_tmp, dec_tmp = np.loadtxt(
            input_file, unpack=True)
        for i in range(0, len(channel_tmp)):
            if group_tmp[i] == int(group_number):
                velocity = np.append(velocity, velocity_tmp[i])
                intensity = np.append(intensity, intensity_tmp[i])
                ra = np.append(ra, ra_tmp[i])
                dec = np.append(dec, dec_tmp[i])

        if len(intensity) == 0:
            intensity = [0]

        if len(velocity) == 0:
            velocity = [0]

        velocitys.append(velocity)
        intensitys.append(intensity)
        ras.append(ra)
        decs.append(dec)

        if len(ra) != 0:
            intensitys_max.append(max(intensity))
            intensitys_min.append(min(intensity))
            max_ra.append(np.max(ra))
            min_ra.append(np.min(ra))
            min_dec.append(np.min(dec))
            max_dec.append(np.max(dec))

    coord_range = max(max(max_ra) - min(min_ra), max(max_dec) - min(min_dec))
    for index in range(0, len(input_files)):
        velocity = velocitys[index]
        intensity = intensitys[index]
        dec = decs[index]
        ra = ras[index]
        title = input_files[index].split(".")[0].upper() + "-" + dates[input_files[index].split(".")[0]]
        if len(velocity) >= 3:
            p1 = [max(intensity), min(velocity) + 0.5 * (max(velocity) - min(velocity)), 0.2]
            p2 = [max(intensity), min(velocity) + 0.5 * (max(velocity) - min(velocity)), 0.3,
                  max(intensity) / 4, min(velocity) + 0.5 * (max(velocity) - min(velocity)), 0.1]
            q = np.linspace(min(velocity), max(velocity), 1000)

            gauss2_groups_for_epoch = gauss2_dict[input_files[index].split(".")[0].upper()]
            if str(group_number) in gauss2_groups_for_epoch:
                try:
                    coeff, var_matrix = curve_fit(gauss2, velocity, intensity, p0=p2, maxfev=100000)
                    hist_fit = gauss2(q, *coeff)
                    ax[0][index].plot(q, hist_fit, 'k')
                except:
                    pass
            else:
                try:
                    coeff, var_matrix = curve_fit(gauss, velocity, intensity, p0=p1, maxfev=100000)
                    hist_fit = gauss(q, *coeff)
                    ax[0][index].plot(q, hist_fit, 'k')
                except:
                    pass

        rel = []
        for i in range(len(velocity) - 1):
            if velocity[i] < min(v_min) or velocity[i] > max(v_max):
                c = (0, 0, 0)
            else:
                c = cm.jet((velocity[i] - min(v_min)) / (max(v_max) - min(v_min)), 1)

            ax[0][index].scatter((velocity[i], velocity[i + 1]), (intensity[i], intensity[i + 1]), color=c, lw=2)
            ax[0][index].set_xlim(min(velocity) - 0.5, max(velocity) + 0.5)
            ax[0][index].xaxis.set_minor_locator(minorLocatorvel)
            ax[0][index].set_title(title)
            ax[0][index].set_xlabel('$V_{\\rm LSR}$ (km s$^{-1}$)')

            el = Circle((ra[i], dec[i]), radius=0.1 * np.sqrt(intensity[i]), angle=0, lw=2)
            ax[1][index].add_artist(el)
            c = cm.jet((velocity[i] - min(v_min)) / (max(v_max) - min(v_min)), 1)
            el.set_facecolor(c)
            rel.append([ra[i], dec[i], velocity[i]])

        ax[0][index].set_ylim((min(intensitys_min)) - 0.1, (max(intensitys_max) + 0.1))
        ax[1][index].set_aspect("equal", adjustable='box')
        ax[1][index].set_xlim(np.mean((max(max_ra), min(min_ra))) - (coord_range/2) - 0.5, (np.mean((max(max_ra), min(min_ra))) - (coord_range/2) - 0.5) + 12)
        ax[1][index].set_ylim(np.mean((max(max_dec), min(min_dec))) - (coord_range/2) - 0.5, np.mean((max(max_dec), min(min_dec))) + (coord_range/2) + 0.5 + 12)
        ax[1][index].set_xlabel('$\\Delta$ RA (mas)')
        ax[1][index].xaxis.set_minor_locator(minorLocatorx)
        ax[1][index].yaxis.set_minor_locator(minorLocatory)
        ax[1][index].invert_xaxis()

    ax[0][0].set_ylabel('Flux density (Jy)')
    ax[1][0].set_ylabel('$\\Delta$ Dec (mas)')
    plt.tight_layout()
    plt.subplots_adjust(top=0.97, bottom=0, wspace=0.18, hspace=0, left=0.05, right=0.99)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot group')
    parser.add_argument('group_number', type=int, help='group number')
    args = parser.parse_args()
    main(args.group_number)
    sys.exit(0)
