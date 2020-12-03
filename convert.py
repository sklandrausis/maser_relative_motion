import sys
import argparse
import numpy as np


def main(input_file, output_file, ref_vel):
    ref_vel = float(ref_vel)
    input_data = np.loadtxt(input_file, unpack=True)
    data_length = input_data[0].shape[0]
    chan = input_data[0]
    vlsr = input_data[1]
    intensity = input_data[2]
    integral = input_data[3]

    reference_index = (np.abs(vlsr - ref_vel)).argmin()
    ra = (input_data[6] - input_data[6][reference_index]) * np.cos(np.radians(input_data[7][reference_index])) * 15000
    dec = (input_data[9] - input_data[9][reference_index]) * 1000

    output_data = []
    for i in range(0, data_length):
        output_data.append([chan[i], vlsr[i], intensity[i], integral[i], ra[i], dec[i]])

    np.savetxt(output_file, output_data)

    sys.exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Visualize pulsar data')
    parser.add_argument('input_file_name', type=str, help='Input file name')
    parser.add_argument('output_file_name', type=str, help='Output file name')
    parser.add_argument('--rev_vel', type=float, help='Reference maser velocity', default=-4944.5)
    args = parser.parse_args()
    main(args.input_file_name, args.output_file_name, args.rev_vel)
