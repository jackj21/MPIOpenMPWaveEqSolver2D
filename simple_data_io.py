import numpy as np

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


# Setup a quick test to see that we can save and load
ny = 5
nx = 10

# Create some floating point data
x = np.arange(ny*nx).astype(np.float32).reshape((ny, nx))

# Save it to a file
save_data("test.data", x)

# Load it from a file
y = load_data("test.data")

print(x)
print(y)

# Check that they match
assert np.all(x == y)
