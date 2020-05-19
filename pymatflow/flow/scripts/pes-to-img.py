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

    parser.add_argument("--shape", type=int, nargs=2, #required=True,
        help="size of the image")

    # ============================================================================
    args = parser.parse_args()

    data = np.loadtxt(args.input)

    # if args.shape is not specified, we guess it from the data
    if args.shape != None:
        shape = args.shape
    else:
        shape = [None, None]
        shape[0] = len(set(data[:, 0]))
        shape[1] = len(set(data[:, 1]))
        print("=============================================\n")
        print("you are not specifying the data shape\n")
        print('we will try to guess it from the data it self\n')


    # ----------------
    # gray scale image
    # ----------------
    img = scale(data[:, 2]).reshape((shape[0], shape[1]))
    cset = plt.imshow(img, cmap="gray")
    plt.colorbar(cset)
    plt.autoscale()
    plt.tight_layout()
    #plt.show()
    plt.savefig(args.output+".gray.png")
    plt.close()

    # ----------------
    # 2D contour plot
    #-----------------
    X = data[:, 0].reshape((shape[0], shape[1]))
    Y = data[:, 1].reshape((shape[0], shape[1]))
    Z = data[:, 2].reshape((shape[0], shape[1]))
    # fill color, three color are divided into three layer(6)
    # cmap = plt.cm.hot means using thermostat plot(graduated red yellow)
    cset = plt.contourf(X, Y, Z, 6, cmap=plt.cm.hot)
    contour = plt.contour(X, Y, Z, [20,40], colors='k')
    plt.colorbar(cset)
    #plt.autoscale()
    #plt.tight_layout()
    #plt.show()
    plt.savefig(args.output+".2d-contour.png")
    plt.close()
    
    # -----------------
    # 3D surface
    # -----------------
    X = data[:, 0].reshape((shape[0], shape[1]))
    Y = data[:, 1].reshape((shape[0], shape[1]))
    Z = data[:, 2].reshape((shape[0], shape[1]))

    ax = plt.axes(projection='3d')
    cset = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='rainbow')
    plt.colorbar(cset)
    #plt.autoscale()
    #plt.tight_layout()
    #plt.show()
    plt.savefig(args.output+".3d-surface.png")
    plt.close()    

    # -----------------
    # 3D Contour
    # -----------------
    X = data[:, 0].reshape((shape[0], shape[1]))
    Y = data[:, 1].reshape((shape[0], shape[1]))
    Z = data[:, 2].reshape((shape[0], shape[1]))

    ax = plt.axes(projection='3d')
    cset = ax.contour3D(X, Y, Z, 50, cmap='rainbow')
    plt.colorbar(cset)
    #plt.autoscale()
    #plt.tight_layout()
    #plt.show()
    plt.savefig(args.output+".3d-contour.png")
    plt.close()    