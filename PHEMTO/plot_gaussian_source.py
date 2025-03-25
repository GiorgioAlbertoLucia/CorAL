import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D

def read_and_create_histogram(filename):
    with open(filename, 'r') as file:
        data = file.readlines()
    
    # Initialize parameters
    dx = dy = dz = xoffset = yoffset = zoffset = None
    ndata = nx = ny = nz = None
    ydy_datablock = []
    
    # Parse file
    reading_data_block = False
    for line in data:
        tokens = line.strip().split()
        if not tokens:
            continue
        if tokens[0] == 'double' and tokens[1] == 'dx':
            dx = float(tokens[2])
        elif tokens[0] == 'double' and tokens[1] == 'dy':
            dy = float(tokens[2])
        elif tokens[0] == 'double' and tokens[1] == 'dz':
            dz = float(tokens[2])
        elif tokens[0] == 'int' and tokens[1] == 'ndata':
            ndata = int(tokens[2])
        elif tokens[0] == 'int' and tokens[1] == 'nx':
            nx = int(tokens[2])
        elif tokens[0] == 'int' and tokens[1] == 'ny':
            ny = int(tokens[2])
        elif tokens[0] == 'int' and tokens[1] == 'nz':
            nz = int(tokens[2])
        elif tokens[0] == 'double' and tokens[1] == 'xoffset':
            xoffset = float(tokens[2])
        elif tokens[0] == 'double' and tokens[1] == 'yoffset':
            yoffset = float(tokens[2])
        elif tokens[0] == 'double' and tokens[1] == 'zoffset':
            zoffset = float(tokens[2])
        elif tokens[0] == 'matrix_double' and tokens[1] == 'ydy_datablock':
            reading_data_block = True
            continue
        elif reading_data_block:
            if tokens[0] == '}':
                reading_data_block = False
                continue
            ydy_datablock.append([float(tokens[0]), float(tokens[1])])

    # Ensure that all required values have been read
    if None in [dx, dy, dz, ndata, nx, ny, nz, xoffset, yoffset, zoffset]:
        raise ValueError("File is missing some parameters.")
    
    # Create histogram array
    histogram = np.zeros((nz, ny, nx))  # Shape (nz, ny, nx)

    # Fill the histogram using the data block
    for i, (value, error) in enumerate(ydy_datablock):
        # Calculate the 3D index from the flat index i
        iz = i // (ny * nx)
        iy = (i % (ny * nx)) // nx
        ix = i % nx
        histogram[iz, iy, ix] = value  # Store only the value; errors can be stored separately if needed

    return histogram

def plot_histogram_to_pdf(histogram, output_pdf='histogram_3d_with_projections.pdf'):
    nz, ny, nx = histogram.shape

    # Calculate 1D projections along each axis
    x_projection = np.sum(histogram, axis=(0, 1))  # Sum over y and z to get x-axis projection
    y_projection = np.sum(histogram, axis=(0, 2))  # Sum over x and z to get y-axis projection
    z_projection = np.sum(histogram, axis=(1, 2))  # Sum over x and y to get z-axis projection

    # Set up PDF document for saving the 3D scatter plot and projections
    with PdfPages(output_pdf) as pdf:
        # Plot the 3D scatter plot
        fig = plt.figure(figsize=(14, 10))
        ax3d = fig.add_subplot(221, projection='3d')

        # Get non-zero bins for 3D scatter plot
        z_indices, y_indices, x_indices = np.nonzero(histogram)
        values = histogram[z_indices, y_indices, x_indices]

        # 3D scatter plot with color indicating particle count
        scatter = ax3d.scatter(x_indices, y_indices, z_indices, c=values, cmap='viridis', s=5)
        colorbar = fig.colorbar(scatter, ax=ax3d, label="Particle Count")
        
        # Label the 3D plot axes
        ax3d.set_xlabel('x')
        ax3d.set_ylabel('y')
        ax3d.set_zlabel('z')
        ax3d.set_title("3D Histogram Scatter Plot of Particle Counts")

        # 1D projection plot for x-axis
        ax_x = fig.add_subplot(222)
        ax_x.plot(np.arange(nx), x_projection, color='b', marker='o', linestyle='-')
        ax_x.set_xlabel('x')
        ax_x.set_ylabel('Particle Count')
        ax_x.set_title('Projection onto x-axis')

        # 1D projection plot for y-axis
        ax_y = fig.add_subplot(223)
        ax_y.plot(np.arange(ny), y_projection, color='g', marker='o', linestyle='-')
        ax_y.set_xlabel('y')
        ax_y.set_ylabel('Particle Count')
        ax_y.set_title('Projection onto y-axis')

        # 1D projection plot for z-axis
        ax_z = fig.add_subplot(224)
        ax_z.plot(np.arange(nz), z_projection, color='r', marker='o', linestyle='-')
        ax_z.set_xlabel('z')
        ax_z.set_ylabel('Particle Count')
        ax_z.set_title('Projection onto z-axis')

        # Save all plots to the PDF
        pdf.savefig(fig)
        plt.close(fig)  # Close the figure to free memory

if __name__ == '__main__':
    filename = 'output/phemto_output_source.dat'
    histogram = read_and_create_histogram(filename)
    plot_histogram_to_pdf(histogram, 'output/histogram.pdf')
