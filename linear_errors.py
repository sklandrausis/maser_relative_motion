import sys
from datetime import datetime
from astropy.io import ascii
from astropy.time import Time
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from parsers.configparser_ import ConfigParser


def get_configs(section, key):
    """

    :param section: configuration file section
    :param key: configuration file sections
    :return: configuration file section key
    """
    config_file_path = "config/config.cfg"
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def convert_datetime_object_to_mjd(time):
    """

    :param time: datetime object
    :return: MJD
    """
    time = time.isoformat()
    tt_ = Time(time, format='isot')
    return tt_.mjd


def main():
    output_file = "output/positionanglemotion_linearity.dat"
    output_data = ascii.read(output_file)
    output_data_headers = output_data.keys()
    velocity = output_data["vel"]
    mjd = np.array([convert_datetime_object_to_mjd(datetime.strptime(date.split("-")[1].strip(), '%d.%m.%Y')) for date
                    in get_configs("parameters", "dates").split(",")])
    mjd = mjd - mjd[0]
    x = np.array([output_data[header] for header in output_data_headers if "x" in header]).T
    y = np.array([output_data[header] for header in output_data_headers if "y" in header]).T
    i = np.array([output_data[header] for header in output_data_headers if "i" in header]).T

    print('PM relative to centres:')
    for epoch in range(0, len(x)):
        x_mean = x[epoch].mean()
        y_mean = y[epoch].mean()
        print("ra dec epoch " + str(epoch + 1) + " mean", x_mean, y_mean)
        x[epoch] = x[epoch] - x_mean
        y[epoch] = y[epoch] - y_mean

    def fun(x, a, b):
        return a * x + b

    f1, a1 = plt.subplots(len(x), 2, sharex="all", squeeze=False)
    f1.subplots_adjust(hspace=0.0, top=0.95, bottom=0.05, left=0.05, right=0.95)
    f1.set_figheight(25, forward=True)
    f1.set_figwidth(30, forward=True)
    f1.set_dpi(80)

    lsvel = []
    ls = []
    lstex = []

    for r in range(0, len(x)):
        a1[r][0].plot(mjd, x[r], ls="", marker="o")
        a1[r][1].plot(mjd, y[r], ls="", marker="o")
        c, m = curve_fit(fun, mjd, x[r])
        a1[r][0].plot(mjd, fun(mjd, c[0], c[1]), lw=1, c="g")
        cdec, mdec = curve_fit(fun, mjd, y[r])
        a1[r][0].plot(mjd[-1], mjd[-1] * c[0] + c[1], ls="", marker="x")
        a1[r][1].plot(mjd[-1], mjd[-1] * cdec[0] + cdec[1], ls="", marker="x")
        a1[r][0].plot(mjd[-1], mjd[-1] * (c[0] + np.sqrt(np.diag(m)[0])) + c[1], ls="", marker="x", color="grey")
        a1[r][1].plot(mjd[-1], mjd[-1] * (cdec[0] + np.sqrt(np.diag(mdec)[0])) + cdec[1], ls="", marker="x", color="grey")
        a1[r][0].plot(mjd[-1], mjd[-1] * (c[0] - np.sqrt(np.diag(m)[0])) + c[1], ls="", marker="x", color="grey")
        a1[r][1].plot(mjd[-1], mjd[-1] * (cdec[0] - np.sqrt(np.diag(mdec)[0])) + cdec[1], ls="", marker="x", color="grey")
        a1[r][1].plot(mjd, fun(mjd, cdec[0], cdec[1]), lw=1, c="g")
        a1[r][0].text(100, np.max(x[r]), "Vlsr %.3f   a_RA %.6f   err_a_RA: %.6f: " % (velocity[r], c[0], np.sqrt(np.diag(m)[0])))
        a1[r][1].text(100, np.min(y[r]), "Vlsr %.3f   a_Dec %.6f  err_a_Dec: %.6f: " % (velocity[r], cdec[0], np.sqrt(np.diag(mdec)[0])))
        a1[r][0].text(3000, np.max(x[r]), "Feature %i" % (r + 1))

        lsvel.append([np.sqrt(((mjd[-1] * c[0] + c[1] - x[r]) / (3886.0 / 365.0)) ** 2 +
                              ((mjd[-1] * cdec[0] + cdec[1] - y[r]) / (3886.0 / 365.0)) ** 2)])

        ls.append([velocity[r], i[r][0], x[r][0], y[r][0], mjd[-1] * c[0] + c[1], mjd[-1] * cdec[0] + cdec[1], mjd[-1] *
                   (c[0] + np.sqrt(np.diag(m)[0])) + c[1], mjd[-1] * (cdec[0] + np.sqrt(np.diag(mdec)[0])) + cdec[1],
                   mjd[-1] * (c[0] - np.sqrt(np.diag(m)[0])) + c[1], mjd[-1] * (cdec[0] - np.sqrt(np.diag(mdec)[0])) +
                   cdec[1], 20 * mjd[-1] * c[0] + c[1], 20 * mjd[-1] * cdec[0] + cdec[1], 20 * mjd[-1] *
                   (c[0] + np.sqrt(np.diag(m)[0])) + c[1], 20 * mjd[-1] * (cdec[0] + np.sqrt(np.diag(mdec)[0])) +
                   cdec[1], 20 * mjd[-1] * (c[0] - np.sqrt(np.diag(m)[0])) + c[1], 20 * mjd[-1] *
                   (cdec[0] - np.sqrt(np.diag(mdec)[0])) + cdec[1]])

        lstex.append([velocity[r], x[r][0], y[r][0], (mjd[-1] * c[0] + c[1] - x[r][0]) / (3886.0 / 365.0),
                      (mjd[-1] * (c[0] + np.sqrt(np.diag(m)[0])) + c[1] - x[r][0] - (mjd[-1] * c[0] + c[1] - x[r][0])) /
                      (3886.0 / 365.0), (mjd[-1] * cdec[0] + cdec[1] - y[r][0]) / (3886.0 / 365.0),
                      (mjd[-1] * (cdec[0] + np.sqrt(np.diag(mdec)[0])) + cdec[1] - y[r][0] -
                       (mjd[-1] * cdec[0] + cdec[1] - y[r][0])) /
                      (3886.0 / 365.0), i[r][0], i[r][1], i[r][2], i[r][3], i[r][4]])

    plt.xlabel("Days")
    a1[3][0].set_ylabel("Shifts in RA [mas]")
    a1[3][1].set_ylabel("Shifts in Dec [mas]")
    a1[0][0].set_title("G78 fit y=a*x+b to RA and Dec shifts")
    a1[0][1].set_title("shifts relative to brightest feature")

    als = np.array(lsvel)
    print("max vel", als.max(), "mas/yr", (als.max() * 1.64 * 150e6) / (365 * 24 * 3600), "km/s")
    print("max vel", als.min(), "mas/yr", (als.min() * 1.64 * 150e6) / (365 * 24 * 3600), "km/s")

    lstexsort = sorted(lstex, key=lambda lstex: lstex[0])
    np.savetxt("output/linearity_errors_fitted_cm.dat", np.array(ls))
    np.savetxt("output/linearity_errors_fitted_tex_cm.dat", np.array(lstex))
    np.savetxt("output/linearity_errors_fitted_tex_sort.dat", np.array(lstexsort))

    plt.show()
    sys.exit(0)


if __name__ == "__main__":
    main()
