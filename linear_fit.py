import sys
from astropy.io import ascii
from scipy.stats import linregress
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


def main():
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

    plt.figure(figsize=(14, 10), dpi=100)
    sub_plots = [plt.subplot(i, aspect='equal') for i in [121, 122]]
    sub_plots[0].set_title("ra fit")
    sub_plots[1].set_title("dec fit")
    for epoch in range(1, number_of_points):
        coefficients_ra = linregress(reference_ras, ras[epoch])
        coefficients_dec = linregress(reference_decs, decs[epoch])
        print("coefficients ra for epoch " + str(epoch),  coefficients_ra)
        print("coefficients dec for epoch " + str(epoch), coefficients_dec)
        sub_plots[0].plot(reference_ras, reference_ras * coefficients_ra.slope + coefficients_ra.intercept, '--k')
        sub_plots[1].plot(reference_ras, reference_ras * coefficients_dec.slope + coefficients_dec.intercept, '--k')

    plt.show()
    sys.exit(0)


if __name__ == "__main__":
    main()
