import netCDF4 as nc
import matplotlib.pyplot as pl

# Function to view the contents of a netCDF file
def view_netcdf(file_path):
    # Open the netCDF file
    dataset = nc.Dataset(file_path, 'r')
    
    # Print the dimensions
    print("Dimensions:")
    for dim in dataset.dimensions.values():
        print(dim)
    
    # Print the variables
    print("\nVariables:")
    for var in dataset.variables.keys():
        print(var)
    
    temp = dataset.variables['TEMP'][:]
    pressure = dataset.variables['PRES'][:]

    # Close the dataset
    dataset.close()
    return pressure,temp

# Example usage
file_path = 'R5905714_238.nc'  # Replace with your netCDF file path
pressure,temp = view_netcdf(file_path)