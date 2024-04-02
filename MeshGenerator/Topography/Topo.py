# %%
import rasterio as rio
import os
import requests
import math
import numpy as np
import pyvista as pv
from scipy.interpolate import griddata

# =============================================================================
# Choosing location and basic information
# =============================================================================
#  More Hall 47.652900105910135, -122.30527631518189
#  Gas works park 47.64547951377896, -122.33632371514638
center_lat       =  47.64547951377896; # this is the latitude of the desired location
center_long      =  -122.33632371514638 # this is the longitude of the desired location

East_dist       = 50; # meters; # This is the istnce we want to extract the topography to the east
West_dist       = 30; # meters; # This is the istnce we want to extract the topography to the west
North_dist      = 50; # meters; # This is the istnce we want to extract the topography to the north
South_dist      = 70; # meters; # This is the istnce we want to extract the topography to the south

base = 0.0 ; # meters; # This is the elevetion of the base of the mesh we use to extrude the topography

PMLThickness      = 5.0                               ; # thickness of the each PML layer
numPMLLayers      = 2                                 ; # number of PML layers
PMLTotalThickness = PMLThickness * numPMLLayers       ; # total thickness of the PML layers
DrmThickness      = 5                                ; # thickness of the DRM layer
numDrmLayers      = 1                                 ; # number of DRM layers
DRMTotalThickness = DrmThickness * numDrmLayers       ; # total thickness of the DRM layers
padLayers         = numPMLLayers + numDrmLayers       ; # number of layers to pad the meshgrid
padThickness      = PMLTotalThickness + DrmThickness  ; # thickness of the padding layers


eleNumZ = 50
eleNumx = 50
eleNumy = 50

# =============================================================================
# helpher functions
# =============================================================================
usgs_url   = "https://portal.opentopography.org/API/usgsdem?datasetName={}&south={}&north={}&west={}&east={}&outputFormat=GTiff&API_Key=0d1cc4f67fd73bb02ab825d895fd27a1"
global_url = "https://portal.opentopography.org/API/globaldem?demtype={}&south={}&north={}&west={}&east={}&outputFormat=GTiff&API_Key=0d1cc4f67fd73bb02ab825d895fd27a1"

def get_OT_GlobalDEM(demtype, bounds, out_fn=None):
    if out_fn is None:
        out_fn = '{}.tif'.format(demtype)
    
    if not os.path.exists(out_fn):
        #Prepare API request url
        #Bounds should be [minlon, minlat, maxlon, maxlat]
        # base_url=global_url
        base_url=usgs_url
        west, south, east, north = bounds
        bnds = (south, north, west, east)
        url = base_url.format(demtype, *bnds)
        print(url)
        #Get
        response = requests.get(url)
        print(response.status_code)
        if response.status_code == 200:
            
            #Write to disk
            open(out_fn, 'wb').write(response.content)



def haversine_distance(lat1, lon1, lat2, lon2):
    R = 6371  # Radius of the Earth in kilometers

    # Convert latitude and longitude from degrees to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    # Calculate differences
    delta_lat = lat2_rad - lat1_rad
    delta_lon = lon2_rad - lon1_rad

    # Haversine formula
    a = math.sin(delta_lat / 2)**2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(delta_lon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    distance = R * c

    return distance






def destination_point(lat, lon, distance, bearing):
    R = 6371  # Radius of the Earth in kilometers

    # Convert latitude and longitude from degrees to radians
    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)

    # Convert bearing from degrees to radians
    bearing_rad = math.radians(bearing)

    # Calculate the new latitude
    new_lat = math.asin(math.sin(lat_rad) * math.cos(distance / R) +
                        math.cos(lat_rad) * math.sin(distance / R) * math.cos(bearing_rad))


    # Calculate the new longitude
    new_lon = lon_rad + math.atan2(math.sin(bearing_rad) * math.sin(distance / R) * math.cos(lat_rad),
                                   math.cos(distance / R) - math.sin(lat_rad) * math.sin(new_lat))

    # Convert new latitude and longitude from radians to degrees
    new_lat_deg = math.degrees(new_lat)
    new_lon_deg = math.degrees(new_lon)

    return new_lat_deg, new_lon_deg





# =============================================================================
# calculate the lat long bounds and get the raster DEM
# =============================================================================
# calculate the bounds
_,maxlong = destination_point(center_lat, center_long, East_dist/1000, 90);   # calculate the maximum longitude
_,minlong = destination_point(center_lat, center_long, West_dist/1000, 270);  # calculate the minimum longitude
maxlat,_  = destination_point(center_lat, center_long, North_dist/1000, 0);   # calculate the maximum latitude
minlat,_  = destination_point(center_lat, center_long, South_dist/1000, 180); # calculate the minimum latitude

bounds  = (minlong, minlat, maxlong, maxlat) ; # creating a bouds to download the topography using API
demtype = "USGS1m"                           ; # the type of the topography we want to download
out_fn  = 'data/{}.tif'.format(demtype)      ; # the name of the file we want to save the topography to     
get_OT_GlobalDEM(demtype, bounds, out_fn)    ; # download the raster file


# =============================================================================
# read the DEM
# =============================================================================
src = rio.open(out_fn)

# create x y z meshgrid for the raster
x    = np.arange(src.width)
y    = np.arange(src.height)
Z    = src.read(1)
X,Y  = np.meshgrid(x,y)         ; # create the meshgrid
Z    = Z[::-1]                  ; # reverse the rows order


# =============================================================================
# Smooth the Z values
# =============================================================================
points = np.vstack([X.ravel(), Y.ravel()]).T ; # create the points for the interpolation

# crete a new grid for the interpolation
x     = np.linspace(x.min(), x.max(), eleNumx)
y     = np.linspace(y.min(), y.max(), eleNumy)
z     = np.linspace(0, 1.0, eleNumZ)
X, Y  = np.meshgrid(x, y)

ZZ = griddata(points, Z.ravel(), (X, Y), method='cubic') ; # interpolate the Z values


# =============================================================================
# padding the meshgrid to incoporate the PML and DRM layers
# =============================================================================
# padding x and y for PML and DRM layers
x  = np.pad(x, (numDrmLayers,numDrmLayers), "linear_ramp", end_values=(x.min()-DRMTotalThickness, x.max()+DRMTotalThickness))
y  = np.pad(y, (numDrmLayers,numDrmLayers), "linear_ramp", end_values=(y.min()-DRMTotalThickness, y.max()+DRMTotalThickness))

# padding the x and y for PML and PML layers
x  = np.pad(x, (numPMLLayers,numPMLLayers), "linear_ramp", end_values=(x.min()-PMLTotalThickness, x.max()+PMLTotalThickness))
y  = np.pad(y, (numPMLLayers,numPMLLayers), "linear_ramp", end_values=(y.min()-PMLTotalThickness, y.max()+PMLTotalThickness))

# padding the Z for PML and DRM layers
ZZ = np.pad(ZZ, (padLayers,padLayers), "edge")

# =============================================================================
# Extrude the topography and create structured grid
# =============================================================================
X, Y, Z = np.meshgrid(x, y, z)

# modify the Z 
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        Z[i,j,:] = np.linspace(base, ZZ[i,j], eleNumZ).ravel()


grid = pv.StructuredGrid(X, Y, Z)
grid["Elevation"] = Z.ravel(order="F")
grid.plot(show_edges=True)









