import sys
from datetime import datetime
from astropy.io import ascii
from astropy.time import Time
from scipy.stats import linregress
import matplotlib.pyplot as plt
import numpy as np

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
    mjd = np.array(
        [convert_datetime_object_to_mjd(datetime.strptime(date.split("-")[1].strip(), '%d.%m.%Y')) for date
         in get_configs("parameters", "dates").split(",")] )
    mjd = mjd - mjd[0]
    ras = []
    decs = []
    fluxes = []
    output_file = "output/output.dat"
    output_data = ascii.read(output_file)
    output_data_headers = output_data.keys()
    velocity = output_data["vel"]
    number_of_points = len(get_configs("parameters", "fileOrder").split(","))
    tmp = 1
    for index in range(1, number_of_points * 3 + 1):
        tmp_data = output_data[output_data_headers[index]]

        if tmp == 1:
            ras.append(tmp_data)
            tmp = 2
        elif tmp == 2:
            decs.append(tmp_data)
            tmp = 3
        else:
            fluxes.append(tmp_data)
            tmp = 1

    reference_ras = ras[0]
    reference_decs = decs[0]
    groups_indexies = []
    group_tmp = []

    for group_index in range(0, len(reference_ras) - 1):
        ra_diff = reference_ras[group_index] - reference_ras[group_index + 1]
        dec_diff = reference_decs[group_index] - reference_decs[group_index + 1]

        if abs(ra_diff) <= 1.27 and abs(dec_diff) <= 1.27:
            group_tmp.append(group_index)
            group_tmp.append(group_index + 1)

        else:
            if sorted(group_tmp) not in groups_indexies:
                groups_indexies.append(sorted(group_tmp))
            group_tmp = []

    for group_index in range(len(reference_ras)-1, 0, -1):
        ra_diff = reference_ras[group_index] - reference_ras[group_index - 1]
        dec_diff = reference_decs[group_index] - reference_decs[group_index - 1]

        if abs(ra_diff) <= 1.27 and abs(dec_diff) <= 1.27:
            group_tmp.append(group_index)
            group_tmp.append(group_index - 1)

        else:
            if sorted(group_tmp) not in groups_indexies:
                groups_indexies.append(sorted(group_tmp))
            group_tmp = []

    ra_avg = {}
    dec_avg = {}
    for group in groups_indexies:
        key = str(group)
        for index in range(0, len(ras)):
            ra_tmp = sum([ras[index][gi] for gi in group])/len(group)
            dec_tmp = sum([decs[index][gi] for gi in group]) / len(group)
            if key not in ra_avg.keys():
                ra_avg[key] = []
                ra_avg[key].append(ra_tmp)
                dec_avg[key] = []
                dec_avg[key].append(dec_tmp)
            else:
                ra_avg[key].append(ra_tmp)
                dec_avg[key].append(dec_tmp)

    plt.figure(figsize=(14, 10), dpi=100)
    sub_plots = [plt.subplot(i, aspect='auto') for i in [121, 122]]
    sub_plots[0].set_title("ra fit")
    sub_plots[1].set_title("dec fit")
    sub_plots[0].set_xlabel("mjd")
    sub_plots[1].set_xlabel("mjd")
    sub_plots[0].set_ylabel("ra")
    sub_plots[1].set_ylabel("dec")
    keys = list(ra_avg.keys())
    alfas_ra = []
    betas_ra = []
    alfas_dec = []
    betas_dec = []
    for epoch in range(0, number_of_points):
        coefficients_ra = linregress(mjd, ra_avg[keys[epoch]])
        coefficients_dec = linregress(mjd, dec_avg[keys[epoch]])
        sub_plots[0].plot(mjd, mjd * coefficients_ra.slope + coefficients_ra.intercept, '--k')
        sub_plots[1].plot(mjd, mjd * coefficients_dec.slope + coefficients_dec.intercept, '--k')
        sub_plots[0].scatter(mjd, ra_avg[keys[epoch]])
        sub_plots[1].scatter(mjd, dec_avg[keys[epoch]])
        alfas_ra.append(coefficients_ra.slope)
        betas_ra.append(coefficients_ra.intercept)
        alfas_dec.append(coefficients_dec.slope)
        betas_dec.append(coefficients_dec.intercept)
        print("coefficients ra for epoch " + str(epoch + 1), coefficients_ra)
        print("coefficients dec for epoch " + str(epoch + 1), coefficients_dec)

    print("\n\n average alfa ra, average beta ra, average alfa dec, average beta dec", np.average(alfas_ra), np.average(betas_ra), np.average(alfas_dec), np.average(betas_dec))
    plt.show()
    sys.exit(0)


if __name__ == "__main__":
    main()
