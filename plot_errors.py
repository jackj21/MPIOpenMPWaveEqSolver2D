import numpy as np
import matplotlib.pyplot as plt

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
        sh = np.fromfile(infile, dtype=np.uint32, count=1)

        # Read the remaining data
        arr = np.fromfile(infile, dtype=np.float32, count=np.prod(sh))

        # Reshape the array to the expected shape
        arr.shape = sh

    return arr

files = ["error_00000251.dat", 
         "error_00000501.dat",
         "error_00001001.dat",
         "error_00002501.dat" ]

for file in files:

    error = load_data(file)

    times = np.linspace(0., 2.5, error.shape[0])

    plt.semilogy(times, error, label=file)

plt.legend()

plt.show()