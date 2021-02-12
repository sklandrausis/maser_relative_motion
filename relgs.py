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
from scipy.stats import stats, pearsonr

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


def main(group_number, epoch, ddddd):
    splits_index = []
    groups = []

    matplotlib.use('TkAgg')
    configuration_items = get_configs_items()
    for key, value in configuration_items.items():
        rcParams[key] = value

    minor_locatorx = MultipleLocator(20)
    minor_locatory = MultipleLocator(20)
    minor_locator_level = MultipleLocator(1)

    input_file = "groups/" + epoch + ".groups"
    date = {date.split("-")[0].strip():
            date.split("-")[1].strip() for date in get_configs("parameters", "dates").split(",")}[epoch]

    if check_if_group_is_in_file(input_file, group_number):
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

        max_max_intensity = max(intensity)
        reference_index = np.where(intensity == max_max_intensity)[0][0]
        references_ra = ra[reference_index]
        references_dec = dec[reference_index]
        references_velocity = velocity[reference_index]
        print("references ra", references_ra, "references dec", references_dec,
              "references velocity", references_velocity)
        ra -= references_ra
        dec -= references_dec

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(16, 16), dpi=90)
        coord_range = max(max(ra) - min(ra), max(dec) - min(dec))

        color = []
        for v in range(len(velocity)):
            if velocity[v] < min(velocity) or velocity[v] > max(velocity):
                c = (0, 0, 0)
            else:
                c = cm.turbo((velocity[v] - min(velocity)) / (max(velocity) - min(velocity)), 1)

            color.append(c)

            el = Circle((ra[v], dec[v]), radius=0.05 * np.log(intensity[v] * 1000), angle=0, lw=2)
            el.set_facecolor(c)
            ax[1].add_artist(el)

        ax[0].scatter(velocity, intensity, color=color, lw=2, picker=True)

        slope, intercept, r_value, p_value, std_err = stats.linregress(ra, dec)
        line = slope * ra + intercept
        ax[1].plot(ra, line, 'r', label='y={:.2f}x+{:.2f} p value {:.2f} std err {:.2f}'.
                   format(slope, intercept, p_value, std_err))

        position_angle2 = 90 + np.degrees(np.arctan(slope))
        print("position angle from linear fit is ", position_angle2)

        print("Distance between fit and points", line-dec)

        print("Pearsonr correlation", pearsonr(ra, line))

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
            ax[1].plot([ra[max_separation["r"]], ra[max_separation["d"]]],
                       [m * ra[max_separation["r"]] + b, m * ra[max_separation["d"]] + b], "k--")

        position_angle = 90 + np.degrees(np.arctan(m))
        print("position angle is ", position_angle)

        if len(velocity) >= 3:
            firs_exceeds_tmp = firs_exceeds(np.diff(velocity), 0.5)
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

            else:
                velocity_tmp = [velocity]
                intensity_tmp = [intensity]
                ra_tmp = [ra]
                dec_tmp = [dec]

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
                        if epoch == "ea063":
                            p = [0.035, -7.75, 0.001]
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
                            ax[0].plot(q, hist_fit, 'k')

                            print("{\\it %d} & %.3f & %.3f & %.1f & %.2f & %.2f & %.3f & %.3f & %.2f & %.2f & %.3f & "
                                  "%.1f(%.1f) & %.3f( ""%.3f)\\\\" %
                                  (gauss_nr, ra_tmp[gauss_nr][max_intensity_index],
                                   dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index],
                                   coeff[1], coeff[2] * 2, intensity[max_intensity_index], coeff[0],
                                   coeff[4], coeff[5] * 2, coeff[3],
                                   max(size), max(size) * 1.64, (velocity[0] - velocity[len(velocity) - 1]) /
                                   max(size), (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64)))

                        elif len(coeff) == 3:
                            hist_fit = gauss(q, *coeff)
                            ax[0].plot(q, hist_fit, 'k')

                            print("{\\it %d} & %.3f & %.3f & %.1f & %.2f & %.2f & %.3f & %.3f & %.1f(%.1f) & %.3f("
                                  "%.3f)\\\\" %
                                  (gauss_nr, ra_tmp[gauss_nr][max_intensity_index],
                                   dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index],
                                   coeff[1], coeff[2] * 2, intensity[max_intensity_index], coeff[0],
                                   max(size), max(size) * 1.64, (velocity[0] - velocity[len(velocity) - 1]) /
                                   max(size), (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64)))

                else:
                    if len(size) > 0:
                        print("{\\it %d} & %.3f & %.3f & %.1f & %s & %s & %.3f & %s & %.1f(%.1f) & %.3f(%.3f)\\\\" %
                              (gauss_nr, ra_tmp[gauss_nr][max_intensity_index],
                               dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index], "-",
                               "-", intensity[max_intensity_index], "-", max(size), max(size) * 1.64,
                               (velocity[0] - velocity[len(velocity) - 1]) / max(size),
                               (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64)))

                    else:
                        print("{\\it %d} & %.3f & %.3f & %.1f & %s & %s & %.3f & %s & %s & %s\\\\" %
                              (gauss_nr, ra_tmp[gauss_nr][max_intensity_index],
                               dec_tmp[gauss_nr][max_intensity_index],
                               velocity[max_intensity_index], "-", "-", intensity[max_intensity_index], "-", "-", "-"))

        ax[1].invert_xaxis()
        ax[1].legend(loc='upper left', bbox_to_anchor=(1.05, 1))
        ax[0].set_xlim(min(velocity) - 0.2, max(velocity) + 0.5)
        ax[0].set_ylim((min(intensity)) - 0.5, (max(intensity) + 0.5))
        ax[0].xaxis.set_minor_locator(minor_locator_level)
        ax[0].set_title(date)
        ax[1].set_aspect("equal", adjustable='box')
        ax[1].set_xlim(np.mean((max(ra), min(ra))) - (coord_range / 2) - 0.5,
                       np.mean((max(ra), min(ra))) + (coord_range / 2) + 0.5)
        ax[1].set_ylim(np.mean((max(dec), min(dec))) - (coord_range / 2) - 0.5,
                       np.mean((max(dec), min(dec))) + (coord_range / 2) + 0.5)
        ax[0].set_ylabel('Flux density (Jy)')
        ax[1].set_ylabel('$\\Delta$ Dec (mas)')
        ax[0].set_xlabel('$V_{\\rm LSR}$ (km s$^{-1}$)')
        ax[1].set_xlabel('$\\Delta$ RA (mas)')

        def onpick1(event):
            ind = event.ind[0]
            if ind not in splits_index:
                splits_index.append(ind)
            max_index = len(velocity) - 1
            if len(groups) == 0:
                groups.append([0, ind + 1])
                groups.append([ind + 1, max_index + 1])

            else:
                for g in groups:
                    if ind in g:
                        new_group_max = max(g)
                        new_group_min = min(g)
                        groups.append([new_group_min, ind + 1])
                        groups.append([ind + 1, new_group_max])
                        groups.remove(g)
                        break

            for g in groups:
                if g[0] < g[1]:
                    index1 = g[0]
                    index2 = g[1]
                elif g[1] < g[0]:
                    index1 = g[1]
                    index2 = g[0]
                else:
                    index1 = g[0]
                    index2 = g[1]

                x = velocity[index1:index2]
                y = intensity[index1:index2]
                if len(x) >= 3:
                    amplitude = max(y)
                    centre_of_peak_index = list(y).index(amplitude)
                    centre_of_peak = x[centre_of_peak_index]
                    standard_deviation = np.std(y)
                    p = [amplitude, centre_of_peak, standard_deviation],
                    coeff, var_matrix = curve_fit(gauss, x, y, p0=p, method="lm", maxfev=100000)
                    q = np.linspace(min(x), max(x), 10000)
                    hist_fit = gauss(q, *coeff)
                    ax[0].plot(q, hist_fit, 'r')
                    event.canvas.draw()

        fig.canvas.mpl_connect('pick_event', onpick1)
        plt.tight_layout()
        plt.subplots_adjust(top=0.947, bottom=0.085, left=0.044, right=0.987, hspace=0.229, wspace=0.182)
        plt.show()
    else:
        print("group is not in epoch")
        sys.exit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot group')
    parser.add_argument('group_number', type=int, help='group number')
    parser.add_argument('epoch', type=str, help='epoch name', choices=["el032", "em064c", "em064d", "es066e", "ea063"])
    parser.add_argument('--d', type=str2bool, help='plot line', default=True)
    args = parser.parse_args()
    main(args.group_number, args.epoch, args.d)
    sys.exit(0)
