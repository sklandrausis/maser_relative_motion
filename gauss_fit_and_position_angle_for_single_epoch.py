import sys
import argparse

import numpy as np
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


def get_configs_items(config_file_path, section):
    """
    :return: None
    """
    config = ConfigParser(config_file_path)
    return config.get_items(section)

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
    cloudlet_dir = get_configs("paths", "cloudlet")
    groups_file_path = get_configs("paths", "groups")
    groups = [[int(g.split(",")[0]), int(g.split(",")[1])] for g in get_configs("groups", epoch).split(";")]
    output = []

    configuration_items = get_configs_items("config/plot.cfg", "main")
    for key, value in configuration_items.items():
        rcParams[key] = value

    minor_locatorx = MultipleLocator(20)
    minor_locatory = MultipleLocator(20)
    minor_locator_level = MultipleLocator(1)

    input_file = groups_file_path + epoch + ".groups"
    date = {date.split("-")[0].strip():
            date.split("-")[1].strip() for date in get_configs("parameters", "dates").split(",")}[epoch]

    if check_if_group_is_in_file(input_file, group_number):
        group_tmp, velocity_tmp, intensity_tmp, ra_tmp, dec_tmp = \
            np.loadtxt(input_file, unpack=True, usecols=(0, 2, 3, 5, 6))

        dtype = [('group_nr', int), ('velocity', float), ('intensity', float), ("ra", float), ("dec", float)]
        values = [(group_tmp[ch], velocity_tmp[ch], intensity_tmp[ch], ra_tmp[ch], dec_tmp[ch])
                  for ch in range(0, len(group_tmp))]
        data = np.array(values, dtype=dtype)
        data = np.sort(data, order=['group_nr', 'velocity'])
        data = data[data["group_nr"] == group_number]

        max_intensity = max(data["intensity"])
        reference_index = np.where(data["intensity"] == max_intensity)[0][0]
        references_ra = data["ra"][reference_index]
        references_dec = data["dec"][reference_index]

        velocity = data["velocity"]
        vel_max = max(velocity)
        vel_min = min(velocity)
        intensity = data["intensity"]
        ra = data["ra"]
        dec = data["dec"]
        ra -= references_ra
        dec -= references_dec

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(16, 16), dpi=100, subplot_kw=dict(box_aspect=1))
        fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=90)
        coord_range = max(max(ra) - min(ra), max(dec) - min(dec))

        slope, intercept, r_value, p_value, std_err = stats.linregress(ra, dec)
        if epoch in get_configs_items("config/config.cfg", "ra_dec_fit_indexes"):
            ra_dec_fit_indexes = [int(index) for index in get_configs("ra_dec_fit_indexes", epoch).split(",")]
            ra_dec_fit_index_a = ra_dec_fit_indexes[0]
            ra_dec_fit_index_b = ra_dec_fit_indexes[1]
            m, b = np.polyfit(ra[ra_dec_fit_index_a:ra_dec_fit_index_b], dec[ra_dec_fit_index_a:ra_dec_fit_index_b], 1)

        else:
            m, b = np.polyfit(ra, dec, 1)

        position_angle = 90 + np.degrees(np.arctan(m))
        position_angle2 = 90 + np.degrees(np.arctan(slope))
        ra_prime = ra * np.cos(np.radians(np.degrees(np.arctan(m)))) + dec * np.sin(
            np.radians(np.degrees(np.arctan(m))))
        dec_prime = dec * np.cos(np.radians(np.degrees(np.arctan(m)))) - ra * np.sin(
            np.radians(np.degrees(np.arctan(m))))

        max_intensity_prime = np.max(intensity)
        reference_index_prime = np.where(intensity == max_intensity_prime)[0][0]
        references_ra_prime = ra_prime[reference_index_prime]
        references_dec_prime = dec_prime[reference_index_prime]
        ra_prime -= references_ra_prime
        dec_prime -= references_dec_prime

        dist_from_prime = np.sqrt((ra - ra_prime) ** 2 + (dec - dec_prime) ** 2)
        references_spot_index = np.where(dist_from_prime == np.min(dist_from_prime))
        distance_between_fit_and_points = np.sqrt((ra_prime[references_spot_index] - ra_prime) ** 2 +
                                                  (dec_prime[references_spot_index] - dec_prime) ** 2)

        for i in range(0, len(distance_between_fit_and_points)):
            if ra_prime[i] < 0:
                distance_between_fit_and_points[i] = (-1) * distance_between_fit_and_points[i]

        m2, b2 = np.polyfit(distance_between_fit_and_points, velocity, 1)

        color = []
        for v in range(0, len(velocity)):
            if velocity[v] < min(velocity) or velocity[v] > max(velocity):
                c = (0, 0, 0)
            else:
                c = cm.turbo((velocity[v] - min(velocity)) / (max(velocity) - min(velocity)), 1)

            color.append(c)

            ra_dec_point = Circle((ra[v], dec[v]), radius=0.05 * np.log(intensity[v] * 1000), angle=0, lw=2,
                                  facecolor=c)
            ax[1].add_artist(ra_dec_point)

        ax[0].scatter(velocity, intensity, color=color, lw=2)
        ax[2].scatter(distance_between_fit_and_points, velocity, color=color, lw=2,
                             s=np.pi * (0.5 * np.log(intensity * 1000)) ** 2)

        if len(groups) == 1:
            ax[2].plot(distance_between_fit_and_points, distance_between_fit_and_points * m2 + b2, "k-")

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
        p_out = []
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
                        #if epoch == "ea063":
                            #p = [0.035, -7.75, 0.001]
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
                        p_out.append(coeff)

                        if len(coeff) == 6:
                            hist_fit = gauss2(q, *coeff)
                            ax[0].plot(q, hist_fit, 'k')

                            print("{\\it %d} & %.3f & %.3f & %.1f & %.2f & %.2f & %.3f & %.3f & %.2f & %.2f & %.3f & "
                                  "%.1f(%.1f) & %.3f( ""%.3f)\\\\" %
                                  (gauss_nr, ra_tmp[gauss_nr][max_intensity_index] + references_ra,
                                   dec_tmp[gauss_nr][max_intensity_index] + references_dec,
                                   velocity[max_intensity_index], coeff[1], coeff[2] * 2,
                                   intensity[max_intensity_index], coeff[0], coeff[4], coeff[5] * 2, coeff[3],
                                   max(size), max(size) * 1.64, (velocity[0] - velocity[len(velocity) - 1]) /
                                   max(size), (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64)))

                            output.append([-1, ra_tmp[gauss_nr][max_intensity_index],
                                           dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index],
                                           coeff[1], coeff[2] * 2, intensity[max_intensity_index], coeff[0],
                                           coeff[4], coeff[5] * 2, coeff[3],
                                           max(size), max(size) * 1.64, (velocity[0] - velocity[len(velocity) - 1]) /
                                           max(size), (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64),
                                           position_angle, position_angle2])

                        elif len(coeff) == 3:
                            hist_fit = gauss(q, *coeff)
                            ax[0].plot(q, hist_fit, 'k',)

                            print("{\\it %d} & %.3f & %.3f & %.1f & %.2f & %.2f & %.3f & %.3f & %.1f(%.1f) & %.3f("
                                  "%.3f)\\\\" %
                                  (gauss_nr, ra_tmp[gauss_nr][max_intensity_index] + references_ra,
                                   dec_tmp[gauss_nr][max_intensity_index] + references_dec, velocity[max_intensity_index],
                                   coeff[1], coeff[2] * 2, intensity[max_intensity_index], coeff[0],
                                   max(size), max(size) * 1.64, (velocity[0] - velocity[len(velocity) - 1]) /
                                   max(size), (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64)))

                            output.append([-1, ra_tmp[gauss_nr][max_intensity_index] + references_ra,
                                           dec_tmp[gauss_nr][max_intensity_index] + references_dec, velocity[max_intensity_index],
                                           coeff[1], coeff[2] * 2, intensity[max_intensity_index], coeff[0], "-",
                                           "-", "-", max(size), max(size) * 1.64,
                                           (velocity[0] - velocity[len(velocity) - 1]) / max(size),
                                           (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64),
                                           position_angle, position_angle2])
                else:
                    if len(size) > 0:
                        print("{\\it %d} & %.3f & %.3f & %.1f & %s & %s & %.3f & %s & %.1f(%.1f) & %.3f(%.3f)\\\\" %
                              (gauss_nr, ra_tmp[gauss_nr][max_intensity_index] + references_ra,
                               dec_tmp[gauss_nr][max_intensity_index] + references_dec, velocity[max_intensity_index],
                               "-", "-", intensity[max_intensity_index], "-", max(size), max(size) * 1.64,
                               (velocity[0] - velocity[len(velocity) - 1]) / max(size),
                               (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64)))

                        output.append([-1, ra_tmp[gauss_nr][max_intensity_index],
                                       dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index],
                                       "-", "-", intensity[max_intensity_index], "-", "-", "-", "-", max(size),
                                       max(size) * 1.64, (velocity[0] - velocity[len(velocity) - 1]) / max(size),
                                       (velocity[0] - velocity[len(velocity) - 1]) / (max(size) * 1.64), position_angle,
                                       position_angle2])

                    else:
                        print("{\\it %d} & %.3f & %.3f & %.1f & %s & %s & %.3f & %s & %s & %s\\\\" %
                              (gauss_nr, ra_tmp[gauss_nr][max_intensity_index] + references_ra,
                               dec_tmp[gauss_nr][max_intensity_index] + references_dec,
                               velocity[max_intensity_index], "-", "-", intensity[max_intensity_index], "-", "-", "-"))

                        output.append([-1, ra_tmp[gauss_nr][max_intensity_index],
                                       dec_tmp[gauss_nr][max_intensity_index], velocity[max_intensity_index],
                                       "-", "-", intensity[max_intensity_index], "-", "-", "-", "-", "-", "-", "-",
                                       position_angle, position_angle2])

        hist_fits = list()
        hist_fits2 = list()
        hist_fits3 = list()
        q2 = np.linspace(min(velocity), max(velocity), 10000)
        r_y = -6.7
        r_y2 = -6.7
        sub_group_nr = 0
        for g in groups:
            index1 = g[0]
            index2 = g[1]

            x = velocity[index1:index2]
            y = intensity[index1:index2]
            q = np.linspace(min(x), max(x), 10000)
            q2 = np.linspace(min(velocity), max(velocity), 10000)
            if len(p_out) == 1:
                p = p_out[0]
                if len(p) == 3:
                    pass
                    # hist_fit = gauss(q2, *p)
                elif len(p) == 6:
                    p_tmp = [p[0:3], p[3:6]]
                    p_tmp_tmp = p_tmp[groups.index(g)]
                    hist_fit = gauss(q2, *p_tmp_tmp)
                    ax[0].plot(q2, hist_fit, '--', c="gray")
            elif len(p_out) == 2:
                p = p_out[groups.index(g)]
                hist_fit = gauss(q2, *p)
                ax[0].plot(q2, hist_fit, '--', c="gray")

            ra_tmp = ra[index1:index2]
            dec_tmp = dec[index1:index2]
            slope, intercept, r_value, p_value, std_err = stats.linregress(ra_tmp, dec_tmp)
            if epoch in get_configs_items("config/config.cfg", "ra_dec_fit_indexes"):
                ra_dec_fit_indexes = [int(index) for index in get_configs("ra_dec_fit_indexes", epoch).split(",")]
                ra_dec_fit_index_a = ra_dec_fit_indexes[0]
                ra_dec_fit_index_b = ra_dec_fit_indexes[1]
                m, b = np.polyfit(ra_tmp[ra_dec_fit_index_a:ra_dec_fit_index_b],
                                  dec_tmp[ra_dec_fit_index_a:ra_dec_fit_index_b], 1)
                ax[1].plot(ra_tmp[ra_dec_fit_index_a:ra_dec_fit_index_b],
                                  ra_tmp[ra_dec_fit_index_a:ra_dec_fit_index_b] * m + b, "k-",
                                  label=str(groups.index(g)))
            else:
                m, b = np.polyfit(ra_tmp, dec_tmp, 1)
                ax[1].plot(ra_tmp, ra_tmp * m + b, "k-", label=str(groups.index(g)))
            ax[1].text(-8, r_y, "r = " + "+ %.3f" % pearsonr(ra_tmp, dec_tmp)[0])

            position_angle = 90 + np.degrees(np.arctan(m))
            position_angle2 = 90 + np.degrees(np.arctan(slope))

            ra_prime_ = ra_tmp * np.cos(np.radians(np.degrees(np.arctan(m)))) + dec_tmp * np.sin(
                np.radians(np.degrees(np.arctan(m))))
            dec_prime_ = dec_tmp * np.cos(np.radians(np.degrees(np.arctan(m)))) - ra_tmp * np.sin(
                np.radians(np.degrees(np.arctan(m))))

            max_intensity_prime_ = np.max(y)
            reference_index_prime_ = np.where(y == max_intensity_prime_)[0][0]
            references_ra_prime_ = ra_prime_[reference_index_prime_]
            references_dec_prime_ = dec_prime_[reference_index_prime_]
            ra_prime_ -= references_ra_prime_
            dec_prime_ -= references_dec_prime_

            dist_from_prime_ = np.sqrt((ra_tmp - ra_prime_) ** 2 + (dec_tmp - dec_prime_) ** 2)
            references_spot_index_ = np.where(dist_from_prime_ == np.min(dist_from_prime_))
            distance_between_fit_and_points_ = np.sqrt((ra_prime_[references_spot_index_] - ra_prime_) ** 2 +
                                                       (dec_prime[references_spot_index_] - dec_prime_) ** 2)

            for i in range(0, len(distance_between_fit_and_points_)):
                if ra_prime_[i] < 0:
                    distance_between_fit_and_points_[i] = (-1) * distance_between_fit_and_points_[i]

            m2, b2 = np.polyfit(distance_between_fit_and_points[index1:index2], x, 1)
            if len(groups) > 1:
                ax[2].plot(distance_between_fit_and_points[index1:index2],
                                  distance_between_fit_and_points[index1:index2] * m2 + b2, "k-")

            ax[2].text(0, r_y2, "r = " + "+ %.3f" % pearsonr(distance_between_fit_and_points[index1:index2], x)[0])

            r_y += 1.5
            r_y2 += 0.5

            max_separation = {"r": 0, "d": -1, "separation": 0}
            sky_coords = [SkyCoord(ra_tmp[coord], dec_tmp[coord], unit=u.arcsec)
                          for coord in range(0, len(ra_tmp))]
            size = []
            max_intensity_index = np.array(y).argmax()

            for j in range(0, len(x)):
                for k in range(j + 1, len(x)):
                    dist = np.sqrt((ra_tmp[j] - ra_tmp[k]) ** 2 + (dec_tmp[j] - dec_tmp[k]) ** 2)
                    size.append(dist)

            for r in range(0, len(ra_tmp)):
                for d in range(0, len(dec_tmp)):
                    if r != d:
                        separation = sky_coords[r].separation(sky_coords[d])
                        if separation > max_separation["separation"]:
                            max_separation["r"] = r
                            max_separation["d"] = d
                            max_separation["separation"] = separation

            output.append([sub_group_nr, ra_tmp[max_intensity_index], dec_tmp[max_intensity_index],
                           x[max_intensity_index], coeff[1], coeff[2] * 2, y[max_intensity_index], coeff[0],
                           "-", "-", "-", max(size), max(size) * 1.64, (x[0] - x[len(x) - 1]) / max(size),
                           (x[0] - x[len(x) - 1]) / (max(size) * 1.64), position_angle, position_angle2])

        sub_group_nr += 1
        q2 = np.linspace(min(velocity), max(velocity), 10000)
        ax2.plot(velocity, intensity - sum(hist_fits2), "k-")
        ax2.plot(velocity, intensity - sum(hist_fits2), "k.", markersize=20)
        ax[0].set_xlim(vel_min - 0.1, vel_max + 0.1)
        ax[0].set_ylim((min(intensity)) - 0.5, (max(intensity) + 0.5))
        ax[0].xaxis.set_minor_locator(minor_locator_level)
        ax[0].set_title(date)
        ax[1].set_aspect("equal", adjustable='box')
        ax[1].set_xlim(np.mean((max(ra), min(ra))) - (coord_range / 2) - 0.5,
                       np.mean((max(ra), min(ra))) + (coord_range / 2) + 0.5)
        ax[1].set_ylim(np.mean((max(dec), min(dec))) - (coord_range / 2) - 0.5,
                       np.mean((max(dec), min(dec))) + (coord_range / 2) + 0.5)
        ax[1].invert_xaxis()
        ax[2].invert_xaxis()
        ax[0].set_ylabel('Flux density (Jy)')
        ax[1].set_ylabel('$\\Delta$ Dec (mas)')
        ax[2].set_ylabel('$V_{\\rm LSR}$ (km s$^{-1}$)')
        ax[0].set_xlabel('$V_{\\rm LSR}$ (km s$^{-1}$)')
        ax[1].set_xlabel('$\\Delta$ RA (mas)')
        ax[2].set_xlabel("offset major axis (mas)")
        ax[1].xaxis.set_minor_locator(minor_locatorx)
        ax[1].yaxis.set_minor_locator(minor_locatory)
        ax2.set_title("Residuals for spectre")
        top = 1.0
        bottom = 0.0
        left = 0.035
        right = 0.97
        hspace = 0.0
        wspace = 0.11
        fig.subplots_adjust(top=top, bottom=bottom, left=left, right=right, hspace=hspace, wspace=wspace)
        plt.show()

        header = ["sub_group_nr", "ra", "dec", "velocity", "vel_fit", "sigma", "max_intensity", "fit_amp", "vel_fit2",
                   "sigma2", "fit_amp2", "max_distance", "max_distance_au", "gradient", "gradient_au",
                   "position_angle", "position_angle2"]
        np.savetxt(cloudlet_dir + "_" + epoch + "_" + str(group_number) + "._sats.csv",
                   np.array(output, dtype=object), delimiter=", ", fmt='%s', header=",".join(header))
    else:
        print("group is not in epoch")
        sys.exit()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot group')
    parser.add_argument('group_number', type=int, help='group number')
    choices = [file.strip().split(".")[0] for file in get_configs("parameters", "fileOrder").split(",")]
    parser.add_argument('epoch', type=str, help='epoch name', choices=choices)
    parser.add_argument('--d', type=str2bool, help='plot line', default=True)
    args = parser.parse_args()
    main(args.group_number, args.epoch, args.d)
    sys.exit(0)
