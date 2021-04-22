import numpy as np
import matplotlib.pyplot as plt

def save_data(filename, arr):
    """
    Saves a 2D floating point array in a simple format.

    * The first 4 bytes are an unsigned integer with the number of rows.
    * The second 4 bytes are an unsigned integer with the number of columns.
    * The remaining data are rows*cols 32-bit floats.

    Parameters
    ----------
    filename: string or path
       The filename to write to
    arr: np.ndarray
       Array of float32s to write out.

    """

    # Open a file handle for binary writing
    with open(filename, "wb") as outfile:

        # Get the arrays shape as unsigned integers and write to file
        sh = np.asarray(arr.shape).astype(np.uint32)
        sh.tofile(outfile)

        # Convert the data to float32 and write it to file
        arr = arr.astype(np.float32)
        arr.tofile(outfile)


def load_data(filename):
    """
    Reads a 2D floating point array in a simple format.

    * The first 4 bytes are an unsigned integer with the number of rows.
    * The second 4 bytes are an unsigned integer with the number of columns.
    * The remaining data are rows*cols 32-bit floats.

    Parameters
    ----------
    filename: string or path
       The filename to write to

    Returns
    -------
    arr: np.ndarray
       Array of float32s read in.

    """

    # Open a file handle for binary reading
    with open(filename, "rb") as infile:

        # Read 2 unsigned integers to get the array's shape
        sh = np.fromfile(infile, dtype=np.uint32, count=2)

        # Read the remaining data
        arr = np.fromfile(infile, dtype=np.float32, count=np.prod(sh))

        # Reshape the array to the expected shape
        arr.shape = sh

    return arr

c = []
c.append(load_data("computed_0.000000.arr"))
c.append(load_data("computed_0.499217.arr"))
c.append(load_data("computed_0.999849.arr"))
c.append(load_data("computed_1.499066.arr"))
c.append(load_data("computed_1.999698.arr"))
c.append(load_data("computed_2.498916.arr"))


t = []
t.append(load_data("true_0.000000.arr"))
t.append(load_data("true_0.499217.arr"))
t.append(load_data("true_0.999849.arr"))
t.append(load_data("true_1.499066.arr"))
t.append(load_data("true_1.999698.arr"))
t.append(load_data("true_2.498916.arr"))

d = []
for cc, tt, in zip(c, t):
    d.append(cc - tt)

for i in range(len(t)):

    plt.figure()
    plt.subplot(131)
    plt.imshow(c[i], clim=(-1,1))
    plt.colorbar()
    plt.subplot(132)
    plt.imshow(t[i], clim=(-1,1))
    plt.colorbar()
    plt.subplot(133)
    plt.imshow(d[i], clim=(-1e-5, 1e-5))
    plt.colorbar()

# computeds = []
# trues = []
# diffs = []

# for k in range(1, 84):

#     c = load_data(f"computed_{k:03}.arr")
#     t = load_data(f"true_{k:03}.arr")

#     computeds.append(c)
#     trues.append(t)
#     diffs.append(c - t)

plt.show()


