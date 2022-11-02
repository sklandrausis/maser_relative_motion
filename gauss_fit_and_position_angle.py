import sys
import argparse

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm, rcParams
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit
from astropy import units as u
from astropy.coordinates import SkyCoord

from parsers.configparser_ import ConfigParser


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def check_if_group_is_in_file(file, group):
    input_file = file
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


def main(group_number, ddddd):
    cloudlet_dir = get_configs("paths", "cloudlet")
    groups_file_path = get_configs("paths", "groups")
    matplotlib.use('TkAgg')
    configuration_items = get_configs_items()
    for key, value in configuration_items.items():
        rcParams[key] = value

    minor_locator_x = MultipleLocator(20)
    minor_locator_y = MultipleLocator(20)
    minor_locator_level = MultipleLocator(1)

    file_order = [file.strip() for file in get_configs("parameters", "fileOrder").split(",")]
    input_files = []

    for file in file_order:
        input_files.append(file)

    dates = {file.split("-")[0].strip(): file.split("-")[1].strip() for file in
             get_configs("parameters", "dates").split(",")}

    bad_files = []
    for file in input_files:
        if not check_if_group_is_in_file(groups_file_path + file.split(".")[0] + ".groups", group_number):
            bad_files.append(file)
            del dates[file.split(".")[0]]
    input_files = [file for file in input_files if file not in bad_files]

    data = dict()
    max_intensity = []
    for index in range(0, len(input_files)):
        epoch = input_files[index].split(".")[0]
        data[epoch] = dict()
        input_file = groups_file_path + input_files[index].split(".")[0] + ".groups"
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
    output2 = []
    for epoch in data:
        print("epoch", epoch)
        velocity = data[epoch]["velocity"]
        intensity = data[epoch]["intensity"]
        index = list(data.keys()).index(epoch)
        ra = data[epoch]["ra"]
        dec = data[epoch]["dec"]
        title = dates[epoch]

        max_separation = {"r": 0, "d": -1, "separation": 0}
        sky_coords = [SkyCoord(ra[coord], dec[coord], unit=u.arcsec) for coord in range(0, len(ra))]
        for r in range(0, len(ra)):
            for d in range(0, len(dec)):
                if r != d:
                    separation = sky_coords[r].separation(sky_coords[d])
                    if separation > max_separation["separation"]:
                        max_separation["r"] = r
                        max_separation["d"] = d
                        max_separation["separation"] = separation

        m, b = np.polyfit([ra[max_separation["r"]], ra[max_separation["d"]]],
                          [dec[max_separation["r"]], dec[max_separation["d"]]], 1)
        if ddddd:
            ax[1][index].plot([ra[max_separation["r"]], ra[max_separation["d"]]],
                              [m * ra[max_separation["r"]] + b, m * ra[max_separation["d"]] + b], "k--")

        position_angle = 90 + np.degrees(np.arctan(m))
        print("position angle is ", position_angle)

        if len(velocity) >= 3:
            firs_exceeds_tmp = firs_exceeds(np.diff(velocity), 0.9)
            split_index = firs_exceeds_tmp + 1
            if firs_exceeds_tmp != -1:
                a = intensity[0:split_index]
                b = intensity[split_index:len(intensity)]
                c = velocity[0:split_index]
                d = velocity[split_index:len(velocity)]
                e = ra[0:split_index]
                f = ra[split_index:len(velocity)]
                g = dec[0:split_index]
                h = dec[split_index:len(velocity)]

                velocity_tmp = [c, d]
                intensity_tmp = [a, b]
                ra_tmp = [e, f]
                dec_tmp = [g, h]

                print(split_index, len(intensity))

            else:
                velocity_tmp = [velocity]
                intensity_tmp = [intensity]
                ra_tmp = [ra]
                dec_tmp = [dec]

            print("number of gauss", len(velocity_tmp))
            for gauss_nr in range(0, len(velocity_tmp)):
                size = []
                max_intensity_index = np.array(intensity_tmp[gauss_nr]).argmax()
                for j in range(0, len(velocity_tmp[gauss_nr])):
                    for k in range(j + 1, len(velocity_tmp[gauss_nr])):
                        dist = np.sqrt((ra[j] - ra[k]) ** 2 + (dec[j] - dec[k]) ** 2)
                        size.append(dist)

                if len(velocity_tmp[gauss_nr]) >= 3:

                    amplitude = max(intensity_tmp[gauss_nr])
                    centre_of_peak_index = list(intensity_tmp[gauss_nr]).index(amplitude)
                    centre_of_peak = velocity_tmp[gauss_nr][centre_of_peak_index]
                    second_largest_amplitude_index = (-intensity_tmp[gauss_nr]).argsort()[1]
                    second_largest_amplitude = intensity_tmp[gauss_nr][second_largest_amplitude_index]
                    second_largest_centre_of_peak = velocity_tmp[gauss_nr][second_largest_amplitude_index]
                    standard_deviation = np.std(intensity_tmp[gauss_nr])
                    ps = [[amplitude, centre_of_peak, standard_deviation],
                          [amplitude, centre_of_peak, standard_deviation, second_largest_amplitude,
                           second_largest_centre_of_peak, standard_deviation], [0.9, -6.45, 0.2],
                          [0.9, -6.45, 0.2, 0.32, -5.43, 0.1], [0.361, -6.98, 0.2, 0.149, -6.489, 0.2],
                          [2.2, -6.9, 0.2, 23.6, -6.22, 0.2], [1.99, -6.977, 0.05, 0.6, -7.3, 0.05],
                          [0.035, -7.75, 0.001]]

                    q = np.linspace(min(velocity_tmp[gauss_nr]), max(velocity_tmp[gauss_nr]), 10000)
                    perrs = []
                    coeffs = []
                    for p in ps:
                        print(p)
                        # if epoch == "ea063":
                        #    p = [0.035, -7.75, 0.001]
                        try:
                            if len(p) == 3:
                                coeff, var_matrix = curve_fit(gauss, velocity_tmp[gauss_nr], intensity_tmp[gauss_nr],
                                                              p0=p, method="lm")
                            else:
                                coeff, var_matrix = curve_fit(gauss2, velocity_tmp[gauss_nr], intensity_tmp[gauss_nr],
                                                              p0=p, method="lm")

                            perr = np.sqrt(np.diag(var_matrix))
                            perr = perr[~np.isnan(perr)]
                            perrs.append(np.mean(perr) / len(perr))
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
                                  (gauss_nr, ra_tmp[gauss_nr][max_intensity_index],
                                   dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index],
                                   coeff[1], coeff[2] * 2, intensity[max_intensity_index], coeff[0],
                                   coeff[4], coeff[5] * 2, coeff[3],
                                   max(size), max(size) * 1.64, (velocity[0] - velocity[len(velocity) - 1]) /
                                   max(size), (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64)))

                            output2.append([epoch, gauss_nr, ra_tmp[gauss_nr][max_intensity_index],
                                            dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index],
                                            coeff[1], coeff[2] * 2, intensity[max_intensity_index], coeff[0],
                                            coeff[4], coeff[5] * 2, coeff[3],
                                            max(size), max(size) * 1.64, (velocity[0] - velocity[len(velocity) - 1]) /
                                            max(size), (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64),
                                            position_angle])

                        elif len(coeff) == 3:
                            hist_fit = gauss(q, *coeff)
                            ax[0][index].plot(q, hist_fit, 'k')
                            print("{\\it %d} & %.3f & %.3f & %.1f & %.2f & %.2f & %.3f & %.3f & %.1f(%.1f) & %.3f("
                                  "%.3f)\\\\" %
                                  (gauss_nr, ra_tmp[gauss_nr][max_intensity_index],
                                   dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index],
                                   coeff[1], coeff[2] * 2, intensity[max_intensity_index], coeff[0],
                                   max(size), max(size) * 1.64, (velocity[0] - velocity[len(velocity) - 1]) /
                                   max(size), (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64)))

                            output2.append([epoch, gauss_nr, ra_tmp[gauss_nr][max_intensity_index],
                                            dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index],
                                            coeff[1], coeff[2] * 2, intensity[max_intensity_index], coeff[0],
                                            "-", "-", "-", max(size), max(size) * 1.64,
                                            (velocity[0] - velocity[len(velocity) - 1]) / max(size),
                                            (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64),
                                            position_angle])
                else:
                    if len(size) > 0:
                        print("{\\it %d} & %.3f & %.3f & %.1f & %s & %s & %.3f & %s & %.1f(%.1f) & %.3f(%.3f)\\\\" %
                              (gauss_nr, ra_tmp[gauss_nr][max_intensity_index],
                               dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index], "-",
                               "-", intensity[max_intensity_index], "-", max(size), max(size) * 1.64,
                               (velocity[0] - velocity[len(velocity) - 1]) / max(size),
                               (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64)))

                        output2.append([epoch, gauss_nr, ra_tmp[gauss_nr][max_intensity_index],
                                        dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index],
                                        "-", "-", intensity[max_intensity_index], "-", "-", "-", "-", max(size),
                                        max(size) * 1.64, (velocity[0] - velocity[len(velocity) - 1]) / max(size),
                                        (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64),
                                        position_angle])

                    else:
                        print("{\\it %d} & %.3f & %.3f & %.1f & %s & %s & %.3f & %s & %s & %s\\\\" %
                              (gauss_nr, ra_tmp[gauss_nr][max_intensity_index],
                               dec_tmp[gauss_nr][max_intensity_index],
                               velocity[max_intensity_index], "-", "-", intensity[max_intensity_index], "-", "-", "-"))

                        output2.append([epoch, gauss_nr, ra_tmp[gauss_nr][max_intensity_index],
                                        dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index], "-",
                                        "-", intensity[max_intensity_index], "-", "-", "-", "-", "-", "-", "-",
                                        position_angle])

        print("\n")

        for o in range(0, len(velocity)):
            output.append([epoch, velocity[o], intensity[o], ra[o], dec[o], position_angle])

        for i in range(len(velocity)):
            if velocity[i] < min(velocity_min) or velocity[i] > max(velocity_max):
                c = (0, 0, 0)
            else:
                c = cm.turbo((velocity[i] - min(velocity_min)) / (max(velocity_max) - min(velocity_min)), 1)

            ax[0][index].scatter((velocity[i]), (intensity[i],), color=c, lw=2)

            el = Circle((ra[i], dec[i]), radius=0.05 * np.log(intensity[i] * 1000), angle=0, lw=2)
            el.set_facecolor(c)
            ax[1][index].add_artist(el)

        ax[0][index].set_xlim(min(velocity_min) - 0.2, max(velocity_max) + 0.5)
        ax[0][index].set_ylim((min(intensity_min)) - 0.5, (max(intensity_max) + 0.5))
        ax[0][index].xaxis.set_minor_locator(minor_locator_level)
        ax[0][index].set_title(title)
        ax[0][index].set_xlabel('$V_{\\rm LSR}$ (km s$^{-1}$)')
        ax[1][index].set_aspect("equal", adjustable='box')
        ax[1][index].set_xlim(np.mean((max(ra_max), min(ra_min))) - (coord_range / 2) - 0.5,
                              np.mean((max(ra_max), min(ra_min))) + (coord_range / 2) + 0.5)
        ax[1][index].set_ylim(np.mean((max(dec_max), min(dec_min))) - (coord_range / 2) - 0.5,
                              np.mean((max(dec_max), min(dec_min))) + (coord_range / 2) + 0.5)
        ax[1][index].set_xlabel('$\\Delta$ RA (mas)')
        ax[1][index].xaxis.set_minor_locator(minor_locator_x)
        ax[1][index].yaxis.set_minor_locator(minor_locator_y)
        ax[1][index].invert_xaxis()
        ax[1][index].set_yscale('linear')
        ax[1][index].set_yscale('linear')

    header1 = ["epoch", "velocity", "intensity", "ra", "dec", "position_angle"]
    header2 = ["epoch", "gauss_nr", "ra", "dec", "velocity", "coeff1", "coeff2_*_2", "max_intensity", "coeff0",
               "coeff4", "coeff5_*_2", "coeff3", "max_distance", "max_distance_au", "gradient", "gradient_au",
               "position_angle"]

    np.savetxt(cloudlet_dir + "cloudlet_" + str(group_number) + "._coords.csv", np.array(output, dtype=object),
               delimiter=", ", fmt='%s', header=",".join(header1))
    np.savetxt(cloudlet_dir + "cloudlet_" + str(group_number) + "._sats.csv", np.array(output2, dtype=object),
               delimiter=", ", fmt='%s', header=",".join(header2))
    ax[0][0].set_ylabel('Flux density (Jy)')
    ax[1][0].set_ylabel('$\\Delta$ Dec (mas)')
    plt.tight_layout()
    plt.subplots_adjust(top=0.947, bottom=0.085, left=0.044, right=0.987, hspace=0.229, wspace=0.182)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot group')
    parser.add_argument('group_number', type=int, help='group number')
    parser.add_argument('--d', type=str2bool, help='plot line', default=True)
    args = parser.parse_args()
    main(args.group_number, args.d)
    sys.exit(0)
