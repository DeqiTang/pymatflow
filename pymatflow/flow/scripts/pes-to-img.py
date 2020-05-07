import copy
import argparse
import numpy as np
import matplotlib.pyplot as plt


def scale(data):
    """
    :param data: an array
    :return out: an array, value transformed to [0, 255]
    """
    out = copy.deepcopy(data)
    div = np.max(data) - np.min(data)
    minimum = np.min(data)
    for i in range(len(out)):
        out[i] = (out[i] - minimum) / div * 255
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str, default=None, required=True,
        help="the pes data file")

    parser.add_argument("-o", "--output", type=str, default="pes.png",
        help="the output image file name")

    parser.add_argument("--shape", type=int, nargs=2, required=True,
        help="size of the image")

    args = parser.parse_args()

    data = np.loadtxt(args.input)
    img = scale(data[:, 2]).reshape((args.shape[0], args.shape[1]))

    plt.imshow(img, cmap="gray")
    plt.autoscale()
    plt.tight_layout()
    #plt.show()
    plt.savefig(args.output)
    