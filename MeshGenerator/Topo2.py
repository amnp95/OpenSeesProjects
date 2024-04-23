# %%
import rasterio as rio
import os
import requests
import math
import numpy as np
import pyvista as pv
from scipy.interpolate import griddata
# from Partition import partition
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import time
os.chdir(os.path.dirname(os.path.abspath(__file__)))
# chekc if the data folder exists
if not os.path.exists("data"):
    os.makedirs("data")
time1 = time.time()
# =============================================================================
# Choosing location and basic information
# =============================================================================
#  More Hall 47.652900105910135, -122.30527631518189
#  Gas works park 47.64547951377896, -122.33632371514638
center_lat       =   47.64547951377896; # this is the latitude of the desired location
center_long      =  -122.33632371514638 # this is the longitude of the desired location

East_dist       = 50; # meters; # This is the istnce we want to extract the topography to the east
West_dist       = 30; # meters; # This is the istnce we want to extract the topography to the west
North_dist      = 50; # meters; # This is the istnce we want to extract the topography to the north
South_dist      = 70; # meters; # This is the istnce we want to extract the topography to the south

base = 0.0 ; # meters; # This is the elevetion of the base of the mesh we use to extrude the topography

PMLThickness      = 2.0                               ; # thickness of the each PML layer
numPMLLayers      = 2                                 ; # number of PML layers
PMLTotalThickness = PMLThickness * numPMLLayers       ; # total thickness of the PML layers
DrmThickness      = 2.0                               ; # thickness of the DRM layer
numDrmLayers      = 1                                 ; # number of DRM layers
DRMTotalThickness = DrmThickness * numDrmLayers       ; # total thickness of the DRM layers
padLayers         = numPMLLayers + numDrmLayers       ; # number of layers to pad the meshgrid
padThickness      = PMLTotalThickness + DrmThickness  ; # thickness of the padding layers


eleNumZ = 200
eleNumx = 200
eleNumy = 200

reg_num_cores = 7
DRM_num_cores = 3
PML_num_cores = 3

info = {
    "RegularDomain": 1,
    "DRMDomain":     2,
    "PMLDomain":     3,
}


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
# create interploation function instead
f = RegularGridInterpolator((x,y),Z.T,method="linear") ; # create the interpolation function
# %%
# =============================================================================
# Smooth the Z values
# =============================================================================
points = np.vstack([X.ravel(), Y.ravel()]).T ; # create the points for the interpolation
ZeleThick = (Z.max() - Z.min()) / eleNumZ
# crete a new grid for the interpolation
x     = np.linspace(x.min(), x.max(), eleNumx)
y     = np.linspace(y.min(), y.max(), eleNumy)
z     = np.linspace(Z.min()-ZeleThick, Z.max(), eleNumZ)
X, Y  = np.meshgrid(x, y)

ZZ = griddata(points, Z.ravel(), (X, Y), method='linear') ; # interpolate the Z values
# %%

# =============================================================================

X,Y,Z = np.meshgrid(x,y,z) ; # create the meshgrid

mesh = pv.StructuredGrid(X,Y,Z) ;
mesh.plot(show_edges=True) # create the mesh


# %%
# =============================================================================

cellcenters = mesh.cell_centers(); # get the cell centers
xy = cellcenters.points[:,[0,1]] ; # get the xy coordinates
z = cellcenters.points[:,2] ; # get the z coordinates
# %%
ztopo = f(xy) ; # get the elevation at the cell centers
# loop through cells and get the elevation
indexes = []
# %%
for i in range(mesh.n_cells):
    if z[i] > ztopo[i]:
        indexes.append(i)
# %%
# remove the cells that are below the topography
 # remove the cells that are below the topography
indcies = np.ones(mesh.n_cells, dtype=bool)
indcies[indexes] = False
mesh = mesh.extract_cells(indcies) ; # remove the cells that are below the topography        
mesh.plot(show_edges=True) # create the mesh

# %%

mesh.plot(show_edges=False) # create the mesh
