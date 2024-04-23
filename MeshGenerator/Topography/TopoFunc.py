import os
import requests
import math

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
