import sys
import argparse

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm, rcParams
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit

from parsers.configparser_ import ConfigParser


def check_if_group_is_in_file(file, group):
    input_file = "groups/" + "/" + file
    group_nr = np.loadtxt(input_file, unpack=True, usecols=0)

    if group not in group_nr:
        return False
    else:
        return True


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


def gauss(x, *p):
    a, b, c = p
    return a * np.exp(-(x - b) ** 2 * np.log(2) / (c ** 2))


def gauss2(x, *p):
    a1, b1, c1, a2, b2, c2 = p
    return a1 * np.exp(-(x - b1) ** 2 * np.log(2) / c1 ** 2) + a2 * np.exp(-(x - b2) ** 2 * np.log(2) / c2 ** 2)


def firs_exceeds(array, value):
    index = -1
    for i in range(0, len(array)):
        if abs(array[i]) > value:
            index = i
            break
    return index


def main(group_number):
    matplotlib.use('TkAgg')
    configuration_items = get_configs_items()
    for key, value in configuration_items.items():
        rcParams[key] = value

    minor_locatorx = MultipleLocator(20)
    minor_locatory = MultipleLocator(20)
    minor_locator_level = MultipleLocator(1)

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

    bad_files = []
    for file in input_files:
        if not check_if_group_is_in_file(file.split(".")[0] + ".groups", group_number):
            bad_files.append(file)
            del dates[file.split(".")[0]]
    input_files = [file for file in input_files if file not in bad_files]

    data = dict()
    max_intensity = []
    for index in range(0, len(input_files)):
        epoch = input_files[index].split(".")[0]
        data[epoch] = dict()
        input_file = "groups/" + "/" + input_files[index].split(".")[0] + ".groups"
        intensity = np.empty(0)
        channels = np.empty(0)
        ra = np.empty(0)
        dec = np.empty(0)
        velocity = np.empty(0)
        group_tmp, channel_tmp, velocity_tmp, intensity_tmp, ra_tmp, dec_tmp = \
            np.loadtxt(input_file, unpack=True, usecols=(0, 1, 2, 3, 5, 6))
        for i in range(0, len(channel_tmp)):
            if group_tmp[i] == int(group_number):
                intensity = np.append(intensity, intensity_tmp[i])
                channels = np.append(channels, channel_tmp[i])
                ra = np.append(ra, ra_tmp[i])
                dec = np.append(dec, dec_tmp[i])
                velocity = np.append(velocity, velocity_tmp[i])

        max_intensity.append(max(intensity))
        data[epoch]["index_for_max_intensity"] = np.where(intensity == max(intensity))[0][0]
        data[epoch]["intensity"] = intensity
        data[epoch]["channels"] = channels
        data[epoch]["velocity"] = velocity
        data[epoch]["ra"] = ra
        data[epoch]["dec"] = dec

    max_max_intensity = max(max_intensity)
    epoch_index_for_max_intensity = max_intensity.index(max_max_intensity)
    epochs = list(data.keys())
    epoch_with_max_intensity = epochs[epoch_index_for_max_intensity]
    print("epoch with max intensity", epoch_with_max_intensity)
    reference_index = data[epoch_with_max_intensity]["index_for_max_intensity"]
    references_ra = data[epoch_with_max_intensity]["ra"][reference_index]
    references_dec = data[epoch_with_max_intensity]["dec"][reference_index]
    references_velocity = data[epoch_with_max_intensity]["velocity"][reference_index]
    print("references ra", references_ra, "references dec", references_dec, "references velocity", references_velocity)

    velocity_max = []
    velocity_min = []
    intensity_max = []
    intensity_min = []
    ra_max = []
    ra_min = []
    dec_max = []
    dec_min = []
    for epoch in data:

        if epoch != epoch_with_max_intensity:
            closet_velocity_index_to_references_velocity = \
                (np.abs(data[epoch]["velocity"] - references_velocity)).argmin()
            data[epoch]["ra"] -= data[epoch]["ra"][closet_velocity_index_to_references_velocity]
            data[epoch]["dec"] -= data[epoch]["dec"][closet_velocity_index_to_references_velocity]

        else:
            data[epoch]["ra"] -= references_ra
            data[epoch]["dec"] -= references_dec

        velocity_max.append(max(data[epoch]["velocity"]))
        velocity_min.append(min(data[epoch]["velocity"]))
        intensity_max.append(max(data[epoch]["intensity"]))
        intensity_min.append(min(data[epoch]["intensity"]))
        ra_max.append(max(data[epoch]["ra"]))
        ra_min.append(min(data[epoch]["ra"]))
        dec_max.append(max(data[epoch]["dec"]))
        dec_min.append(min(data[epoch]["dec"]))

    fig, ax = plt.subplots(nrows=2, ncols=len(input_files), figsize=(16, 16), dpi=90)
    coord_range = max(max(ra_max) - min(ra_min), max(dec_max) - min(dec_min))
    output = []
    for epoch in data:
        print("epoch", epoch)
        velocity = data[epoch]["velocity"]
        intensity = data[epoch]["intensity"]
        index = list(data.keys()).index(epoch)
        ra = data[epoch]["ra"]
        dec = data[epoch]["dec"]
        title = epoch.upper() + "-" + dates[epoch]

        if len(velocity) >= 3:
            firs_exceeds_tmp = firs_exceeds(np.diff(velocity), 0.5)
            split_index = firs_exceeds_tmp + 1
            if firs_exceeds_tmp != -1:
                a = intensity[0:split_index]
                b = intensity[split_index:len(intensity)]
                c = velocity[0:split_index]
                d = velocity[split_index:len(velocity)]

                velocity_tmp = [c, d]
                intensity_tmp = [a, b]

            else:
                velocity_tmp = [velocity]
                intensity_tmp = [intensity]

            size = []
            max_intensity_index = np.array(intensity).argmax()
            for j in range(0, len(velocity)):
                for k in range(j + 1, len(velocity)):
                    dist = np.sqrt((ra[j] - ra[k]) ** 2 + (dec[j] - dec[k]) ** 2)
                    size.append(dist)

            for gauss_nr in range(0, len(velocity_tmp)):
                if len(velocity_tmp[gauss_nr]) >= 3:
                    '''
                    p1 = [max(intensity_tmp[gauss_nr]), min(velocity_tmp[gauss_nr]) +
                          0.5 * (max(velocity_tmp[gauss_nr]) - min(velocity_tmp[gauss_nr])), 0.2]
                    p2 = [max(intensity_tmp[gauss_nr]), min(velocity_tmp[gauss_nr]) +
                          0.5 * (max(velocity_tmp[gauss_nr]) - min(velocity_tmp[gauss_nr])), 0.3,
                          max(intensity_tmp[gauss_nr]) / 4, min(velocity_tmp[gauss_nr]) +
                          0.5 * (max(velocity_tmp[gauss_nr]) - min(velocity_tmp[gauss_nr])), 0.1] 
                    '''

                    amplitude = max(intensity_tmp[gauss_nr])
                    centre_of_peak_index = list(intensity_tmp[gauss_nr]).index(amplitude)
                    centre_of_peak = velocity_tmp[gauss_nr][centre_of_peak_index]
                    second_largest_amplitude_index = (-intensity_tmp[gauss_nr]).argsort()[1]
                    second_largest_amplitude = intensity_tmp[gauss_nr][second_largest_amplitude_index]
                    second_largest_centre_of_peak = velocity_tmp[gauss_nr][second_largest_amplitude_index]
                    standard_deviation = np.std(intensity_tmp[gauss_nr])
                    p1 = [amplitude, centre_of_peak, standard_deviation]
                    p2 = [amplitude, centre_of_peak, standard_deviation,
                          second_largest_amplitude, second_largest_centre_of_peak, standard_deviation]
                    q = np.linspace(min(velocity_tmp[gauss_nr]), max(velocity_tmp[gauss_nr]), 10000)

                    perrs = []
                    coeffs = []
                    try:
                        coeff, var_matrix = curve_fit(gauss2, velocity_tmp[gauss_nr], intensity_tmp[gauss_nr],
                                                      p0=p2, maxfev=100000, method="lm")
                        perr = np.sqrt(np.diag(var_matrix))
                        perr = perr[~np.isnan(perr)]
                        perrs.append(sum(perr) / len(perr))
                        coeffs.append(coeff)
                    except:
                        pass

                    try:
                        coeff, var_matrix = curve_fit(gauss, velocity_tmp[gauss_nr], intensity_tmp[gauss_nr],
                                                      p0=p1, maxfev=100000, method="lm")
                        perr = np.sqrt(np.diag(var_matrix))
                        perr = perr[~np.isnan(perr)]
                        perrs.append(sum(perr) / len(perr))
                        coeffs.append(coeff)
                    except:
                        pass

                    if len(perrs) > 0:
                        coeff_index = perrs.index(min(perrs))
                        coeff = coeffs[coeff_index]

                        if len(coeff) == 6:
                            hist_fit = gauss2(q, *coeff)
                            ax[0][index].plot(q, hist_fit, 'k')

                            print("{\\it %d} & %.3f & %.3f & %.1f & %.2f & %.2f & %.3f & %.3f & %.2f & %.2f & %.3f & "
                                  "%.1f(%.1f) & %.3f( ""%.3f)\\\\" %
                                  (group_number, ra[max_intensity_index], dec[max_intensity_index],
                                   velocity[max_intensity_index],
                                   coeff[1], coeff[2] * 2, intensity[max_intensity_index], coeff[0],
                                   coeff[4], coeff[5] * 2, coeff[3],
                                   max(size), max(size) * 1.64, (velocity[0] - velocity[len(velocity) - 1]) /
                                   max(size), (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64)))

                        elif len(coeff) == 3:
                            hist_fit = gauss(q, *coeff)
                            ax[0][index].plot(q, hist_fit, 'k')

                            print("{\\it %d} & %.3f & %.3f & %.1f & %.2f & %.2f & %.3f & %.3f & %.1f(%.1f) & %.3f("
                                  "%.3f)\\\\" %
                                  (group_number, ra[max_intensity_index], dec[max_intensity_index],
                                   velocity[max_intensity_index],
                                   coeff[1], coeff[2] * 2, intensity[max_intensity_index], coeff[0],
                                   max(size), max(size) * 1.64, (velocity[0] - velocity[len(velocity) - 1]) /
                                   max(size), (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64)))
                else:
                    if len(size) > 0:
                        print("{\\it %d} & %.3f & %.3f & %.1f & %s & %s & %.3f & %s & %.1f(%.1f) & %.3f(%.3f)\\\\" %
                              (group_number, ra[max_intensity_index], dec[max_intensity_index],
                               velocity[max_intensity_index], "-",
                               "-", intensity[max_intensity_index], "-", max(size), max(size) * 1.64,
                               (velocity[0] - velocity[len(velocity) - 1]) / max(size),
                               (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64)))

                    else:
                        print("{\\it %d} & %.3f & %.3f & %.1f & %s & %s & %.3f & %s & %s & %s\\\\" %
                              (group_number, ra[max_intensity_index], dec[max_intensity_index],
                               velocity[max_intensity_index], "-", "-", intensity[max_intensity_index], "-", "-", "-"))

        for o in range(0, len(velocity)):
            output.append([epoch, velocity[o], intensity[o], ra[o], dec[o]])

        np.savetxt(str(group_number) + ".txt", output, delimiter=" ", fmt='%s')

        for i in range(len(velocity) - 1):
            if velocity[i] < min(velocity_min) or velocity[i] > max(velocity_max):
                c = (0, 0, 0)
            else:
                c = cm.turbo((velocity[i] - min(velocity_min)) / (max(velocity_max) - min(velocity_min)), 1)

            ax[0][index].scatter((velocity[i], velocity[i + 1]), (intensity[i], intensity[i + 1]), color=c, lw=2)

            el = Circle((ra[i], dec[i]), radius=0.05 * np.log(intensity[i] * 1000), angle=0, lw=2)
            el.set_facecolor(c)
            ax[1][index].add_artist(el)

        ax[0][index].set_xlim(min(velocity_min) - 0.5, max(velocity_max) + 0.5)
        ax[0][index].set_ylim((min(intensity_min)) - 0.7, (max(intensity_max) + 0.7))
        ax[0][index].xaxis.set_minor_locator(minor_locator_level)
        ax[0][index].set_title(title)
        ax[0][index].set_xlabel('$V_{\\rm LSR}$ (km s$^{-1}$)')
        ax[1][index].set_aspect("equal", adjustable='box')
        ax[1][index].set_xlim(np.mean((max(ra_max), min(ra_min))) - (coord_range / 2) - 0.5,
                              np.mean((max(ra_max), min(ra_min))) + (coord_range / 2) + 0.5)
        ax[1][index].set_ylim(np.mean((max(dec_max), min(dec_min))) - (coord_range / 2) - 0.5,
                              np.mean((max(dec_max), min(dec_min))) + (coord_range / 2) + 0.5)
        ax[1][index].set_xlabel('$\\Delta$ RA (mas)')
        ax[1][index].xaxis.set_minor_locator(minor_locatorx)
        ax[1][index].yaxis.set_minor_locator(minor_locatory)
        ax[1][index].invert_xaxis()
        ax[1][index].set_yscale('linear')
        ax[1][index].set_yscale('linear')

    ax[0][0].set_ylabel('Flux density (Jy)')
    ax[1][0].set_ylabel('$\\Delta$ Dec (mas)')
    plt.tight_layout()
    plt.subplots_adjust(top=0.947, bottom=0.085, left=0.044, right=0.987, hspace=0.229, wspace=0.182)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot group')
    parser.add_argument('group_number', type=int, help='group number')
    args = parser.parse_args()
    main(args.group_number)
    sys.exit(0)
