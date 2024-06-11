import matplotlib.pyplot as plt
import numpy as np
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D

class VisualizeTrajectory:
    def __init__(self, data, planet, data_list, labels):
        """
        Initialize the VisualizeTrajectory object with the trajectory data.
        
        Parameters:
            data (pandas.DataFrame): Dataframe containing the trajectory data with columns 'x', 'y', 'z', and 'time'.
        """
        self.data = data
        self.planet = planet
        self.data_list = data_list
        self.labels = labels

    def plot_trajectory(self, scale):
        """
        Plot the x, y, and z coordinates of the spacecraft over time on separate subplots within one figure.
        """
        # Create a figure with three subplots
        fig, axs = plt.subplots(4, 1, figsize=(10, 15))

        # Plot X coordinate over time
        axs[0].plot(self.data['time'], self.data['x'] / scale, label='X, км', color='blue')
        # axs[0].set_title('X')
        axs[0].set_xlabel('Время, с')
        axs[0].grid(True)
        axs[0].set_ylabel('X, км')
        # axs[0].legend()

        # Plot Y coordinate over time
        axs[1].plot(self.data['time'], self.data['y'] / scale, label='Y, км', color='red')
        # axs[1].set_title('Y Coordinate Over Time')
        axs[1].set_xlabel('Время, с')
        axs[1].set_ylabel('Y, км')
        # axs[1].legend()
        axs[1].grid(True)

        # Plot Z coordinate over time
        axs[2].plot(self.data['time'], self.data['z']  / scale, label='Z, км', color='green')
        # axs[2].set_title('Z Coordinate Over Time')
        axs[2].set_xlabel('Время, с')
        axs[2].set_ylabel('Z, км')
        # axs[2].legend()
        axs[2].grid(True)
        
        # Plot Z coordinate over time
        axs[3].plot(self.data['time'], self.data['alt'] / scale, label='Высота, км', color='black')
        # axs[3].set_title('Z Coordinate Over Time')
        axs[3].set_xlabel('Время, с')
        axs[3].set_ylabel('Высота, км')
        # axs[3].legend()
        axs[3].grid(True)

        # Improve layout to prevent labels from overlapping
        plt.tight_layout()

        # Display the plot
        plt.show()
        
    
    
    def plot_trajectory_over_planet(self):
        image = mpimg.imread('map_Earth.jpg')
        
        plt.figure(figsize=(10, 6))
        for data, label, color in zip(self.data_list, self.labels, ['red', 'black']):
            plt.scatter(data['longitude'], data['latitude'], label=f'{label}', s=1, color=color)
        
        plt.xlabel('Долгота (град)')
        plt.ylabel('Широта (град)')
        plt.title('Траектория КА вокруг Земли')
        plt.grid(True)
        plt.legend()
        plt.imshow(image, extent=[-180, 180, -90, 90])
        plt.show()
    
    def plot_3d_trajectory(self):
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Earth parameters
        earth_radius = 6371  # in km, approximate radius of Earth

        # Create a sphere (Earth)
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = earth_radius * np.outer(np.cos(u), np.sin(v))
        y = earth_radius * np.outer(np.sin(u), np.sin(v))
        z = earth_radius * np.outer(np.ones(np.size(u)), np.cos(v))
        
        # Plot Earth
        ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b', alpha=0.5)  # semi-transparent

        # Normalize time for color mapping
        normalized_time = (self.data['time'] - self.data['time'].min()) / (self.data['time'].max() - self.data['time'].min())
        
        # Create a scatter plot of the trajectory
        scatter = ax.scatter(self.data['x'], self.data['y'], self.data['z'], c=normalized_time, cmap='viridis', s=10)

        ax.set_xlabel('X км')
        ax.set_ylabel('Y км')
        ax.set_zlabel('Z км')
        ax.set_title('3D Траектория')
        ax.set_aspect('equal')
        # Add a colorbar
        cbar = fig.colorbar(scatter, ax=ax, pad=0.1)
        cbar.set_label('Нормализованное время')

        plt.legend()
        plt.show()

    

