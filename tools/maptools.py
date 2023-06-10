# numpy for array manipulations
import numpy as np
# the workhorse, provides routines to download tiles and do some
# coordinate transforms
import cartopy.io.img_tiles as cartotiles
# image manipulating and saving
from PIL import Image
# polygons to define regions
from shapely.geometry.polygon import Polygon
# more coordinate support
import sys
sys.path.insert(1, "../tools")
import globalmaptiles as gmaptiles
# for working with metadata in csv
import pandas as pd
# for parallel downloads
import concurrent.futures



def tile_coords(lon, lat, zoom):
    """ 
    Converts a lat, lon coordinate to the (x,y,z) coordinates
    of the tile containing lat, lon at the given zoom level
    for calculations, see: https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames
    
    Parameters
    ----------
    lon (float):
        The longitude of the point in degrees (decimal)
    lat (float):
        The latitude of the point in degrees (decimal), between -85.0511 and +85.0511
    zoom (int):
        The zoom level, should be between 0 and 22 for mapbox, 0 and 18 for Openstreetmap
        
    Returns
    -------
    slipmapcoordinates (Tuple):
        Tuple (x,y,z) of ints corresponding to the tile coordinates on a slipmap at zoom level z
    """
    # conversion factors
    DEG_TO_RAD = np.pi/180
    n = 2**zoom
    
    xtile = n * ((lon+180)/360)
    ytile = n * (1 - (np.log(np.tan(lat*DEG_TO_RAD) + 1/np.cos(lat*DEG_TO_RAD))/np.pi))/2
    
    return (int(np.floor(xtile)), int(np.floor(ytile)), zoom)



def area_tiles(lon1, lat1, lon2, lat2, zoom):
    """ 
    Calculates the (x,y,z) coordinates of all tiles overlapping with the rectangle
    definide by (lon1,lat1), (lon2,lat2) at zoom level zoom
    
    Parameters
    lon1, lon2 (float):
        longitude of coordinates defining rectangle in degrees
    lat1, lat2 (float):
        latitiude of coordinates defining rectangle in degrees
    zoom (int):
        zoomlevel
    """
    lon_min = min(lon1,lon2)
    lon_max = max(lon1,lon2)
    lat_min = min(lat1,lat2)
    lat_max = max(lat1,lat2)
    
    # get tile coordinates for (min,min) and (max,max)
    # coordinates start from north pole!
    x_min, y_min, z = tile_coords(lon_min, lat_max, zoom)
    x_max, y_max, z = tile_coords(lon_max, lat_min, zoom)
    
    # construct all tile coords and collect them in matrix
    x_coords = np.arange(x_min, x_max+1, 1)
    y_coords = np.arange(y_min, y_max+1, 1)
    
    return (x_coords, y_coords)



def download_tiles(x_coords, y_coords, zoom):
    """ 
    Downloads tiles in bbox determined by area_tiles
    """
    # fetch all tiles
    tiles = []
    tile_imgs = []
    print("rows: %d" %len(y_coords))
    print("tiles per row: %d" %len(x_coords))
    for y_coord in y_coords:
        row = []
        row_imgs = []
        print("row: %d" %(y_coord-y_coords[0]+1))
        for x_coord in x_coords:
            tile = maptiles.get_image((x_coord, y_coord, zoom))
            tile_img = tile[0]
            row.append(tile)
            row_imgs.append(tile_img)
        tiles.append(row)
        tile_imgs.append(row_imgs)
        
    return (tiles, tile_imgs)
    


def stitch_image(tile_imgs):
    """ 
    stitches tiles downloaded by download_tiles together into 1 big map
    """
    # stitch image together
    rows = len(tile_imgs)
    columns = len(tile_imgs[0])
    tile_size = tile_imgs[0][0].size[0]

    # create new, empty, image
    big_img = Image.new("RGB", (columns*tile_size, rows*tile_size))

    # fill in the image
    for row in range(rows):
        for column in range(columns):
            big_img.paste(tile_imgs[row][column], box=(column*tile_size, row*tile_size))
    
    return big_img



# fetches all tiles for a given domain

# copied from sourcefile because the image stitcher stitches tiles together wrongly, mapboxtiles have NO overlap

def images_for_domain(maptilesinstance, target_domain, target_z):
    tiles = []

    def fetch_tile(tile):
        try:
            img, extent, origin = maptilesinstance.get_image(tile)
        except OSError:
            # Some services 404 for tiles that aren't supposed to be
            # there (e.g. out of range).
            raise
        img = np.array(img)
        #x = np.linspace(extent[0], extent[1], img.shape[1])
        #y = np.linspace(extent[2], extent[3], img.shape[0])
        return img, tile, origin

    with concurrent.futures.ThreadPoolExecutor(
            max_workers=maptilesinstance._MAX_THREADS) as executor:
        futures = []
        for tile in maptilesinstance.find_images(target_domain, target_z):
            futures.append(executor.submit(fetch_tile, tile))
        for future in concurrent.futures.as_completed(futures):
            try:
                img, tile, origin = future.result()
                tiles.append([img, tile, origin])
            except OSError:
                pass

    return tiles



# function to merge the tiles recieved from images_for_domain
def merge_tiles(tiles, image_pixel_size, area_tile_coords, TILESIZE):
    img = np.zeros((image_pixel_size[1], image_pixel_size[0], 3),np.uint8)
    origin = (area_tile_coords[0][0], area_tile_coords[1][0])

    for tile in tiles:
        # get the tilecoords in the map:
        mtc = (tile[1][0]-origin[0], tile[1][1]-origin[1])
        # get the pixel coords of the upperleft pixel
        ptc = (mtc[0]*TILESIZE, mtc[1]*TILESIZE)
        #print(ptc)
        # fill in the image, array is transpose of xy-coords for maptiles!
        img[ptc[1]:ptc[1]+TILESIZE,ptc[0]:ptc[0]+TILESIZE] = tile[0]
    
    return img



def gen_map_metadata(REGION_FILE, REGION_NAME, TOKEN, USER, DPI, TILESIZE=256, SEP=","):
    # collect data from csv
    regionfile = pd.read_csv(REGION_FILE, sep=SEP)
    regionfile = regionfile.set_index("Region")

    MAP_PARAMS = regionfile.loc[REGION_NAME]
    MAP_PARAMS = regionfile.loc[REGION_NAME]
    ZOOM = int(MAP_PARAMS["ZoomLevel"])
    LON = (MAP_PARAMS["LonMin"], MAP_PARAMS["LonMax"])
    LAT = (MAP_PARAMS["LatMin"], MAP_PARAMS["LatMax"])
    STYLE_NAME = MAP_PARAMS["StyleName"]
    STYLE_ID = MAP_PARAMS["StyleID"]
    CUSTOM_TILESET = MAP_PARAMS["CustomTileset"]
    FOLDER = MAP_PARAMS["Folder"]
    MAPNAME = "%s_%s"%(REGION_NAME, ZOOM)
    
    # if not custom map, overwrite tilesize to 512
    if not CUSTOM_TILESET:
        TILESIZE = 512

    # instance to download tiles
    if CUSTOM_TILESET:
        maptiles = cartotiles.MapboxStyleTiles(TOKEN, USER, STYLE_ID, cache=False)
    else:
        maptiles = cartotiles.MapboxTiles(TOKEN, STYLE_ID)

    # you shouldn't need to change these ;)
    CM_PER_INCH = 2.54
    EARTH_CIRCUMFERENCE = 40_075_016.686   # meter

    ###################################
    # Calculate all required metadata #
    ###################################

    # coordinate transform instance
    map_coords = gmaptiles.GlobalMercator(tileSize=TILESIZE)
    # constant
    PIXEL_SCALE_EQUATOR = EARTH_CIRCUMFERENCE/(2**ZOOM*TILESIZE)  # meter/pixel

    # construct polygon of area
    corners = [(LON[0], LAT[0]), (LON[0], LAT[1]), (LON[1], LAT[1]), (LON[1], LAT[0])]
    area = Polygon([ (np.mean(maptiles.tile_bbox(*tile_coords(*corner, ZOOM))[0]), np.mean(maptiles.tile_bbox(*tile_coords(*corner, ZOOM))[1]) ) for corner in corners])
    # find all tiles intersecting the polygon
    area_tile_coords = area_tiles(LON[0], LAT[0], LON[1], LAT[1], ZOOM) 

    # get extent and origin of map in different coordinate systems
    
    # extent in tile-ids
    extent_tile_id = [area_tile_coords[0][0], area_tile_coords[0][-1],
                    area_tile_coords[1][0], area_tile_coords[1][-1]]
    # extent in pseudo-mercator coordinates
    extent1 = maptiles.tile_bbox(*tile_coords(LON[0],LAT[0],ZOOM))
    extent2 = maptiles.tile_bbox(*tile_coords(LON[1],LAT[1],ZOOM))
    meterx = np.concatenate([extent1[0], extent2[0]])
    metery = np.concatenate([extent1[1], extent2[1]])
    extent_pseudo_mercator = (min(meterx), max(meterx),
                            min(metery), max(metery))
    del extent1, extent2, meterx, metery
    # extent in lat-lon degrees
    lat_min, lon_min = map_coords.MetersToLatLon(extent_pseudo_mercator[0], extent_pseudo_mercator[2])
    lat_max, lon_max = map_coords.MetersToLatLon(extent_pseudo_mercator[1], extent_pseudo_mercator[3])
    extent_lat_lon = (lat_min, lat_max, lon_min, lon_max)
    del lat_min, lon_min, lat_max, lon_max
    # extent in pixels
    min_point = (map_coords.MetersToPixels(extent_pseudo_mercator[0], extent_pseudo_mercator[2], ZOOM))
    max_point = (map_coords.MetersToPixels(extent_pseudo_mercator[1], extent_pseudo_mercator[3], ZOOM))
    map_pixel_extent = (int(min_point[0]), int(max_point[0]), int(min_point[1]), int(max_point[1]))
    map_pixel_origin = (int(min_point[0]), int(min_point[1]))
    # cleanup
    del min_point, max_point

    #####################################
    # Calculate extra info to show user #
    #####################################

    # precalculate the size of the resulting image in pixels
    image_pixel_size = (len(area_tile_coords[0])*TILESIZE,
                        len(area_tile_coords[1])*TILESIZE)

    # calculate the scale of the image (top and bottom, in meter/pixel)
    pixel_scale_bottom = np.cos(LAT[0]*np.pi/180)*PIXEL_SCALE_EQUATOR
    pixel_scale_top = np.cos(LAT[1]*np.pi/180)*PIXEL_SCALE_EQUATOR

    # calculate the physical size of the image
    image_phys_size = (image_pixel_size[0]/DPI*CM_PER_INCH,image_pixel_size[1]/DPI*CM_PER_INCH)

    # calculate how much larger the actual geographical area displayed is
    point0 = map_coords.LatLonToMeters(LAT[0],LON[0])
    point1 = map_coords.LatLonToMeters(LAT[1],LON[1])
    # area of region of interest in square meters (map coordinates)
    meter_area_of_interest = (point1[0]-point0[0])*(point1[1]-point0[1])
    # actual area covered by map (square meters, map coordinates)
    meter_area_mapped = (extent_pseudo_mercator[1]-extent_pseudo_mercator[0])*\
                        (extent_pseudo_mercator[3]-extent_pseudo_mercator[2])
    # calculate ratio
    surplus_area = 1-meter_area_of_interest/meter_area_mapped
    del point0, point1

    # create data array
    data = {"MapName": MAPNAME,
            "Region" : REGION_NAME, 
            "StyleName" : STYLE_NAME, 
            "StyleID" : STYLE_ID, 
            "IsPrivateStyle" : CUSTOM_TILESET, 
            "user" : USER, 
            "LatMin" : LAT[0], 
            "LatMax" : LAT[1], 
            "LonMin" : LON[0], 
            "LonMax" : LON[1], 
            "Zoom" : ZOOM, 
            "TileSize" : TILESIZE, 
            "PixWidth" : image_pixel_size[0], 
            "PixHeight" : image_pixel_size[1], 
            "PixScaleTop" : pixel_scale_top, 
            "PixScaleBottom" : pixel_scale_bottom,
            "MeterXMin" : extent_pseudo_mercator[0], 
            "MeterXMax" : extent_pseudo_mercator[1], 
            "MeterYMin" : extent_pseudo_mercator[2], 
            "MeterYMax" : extent_pseudo_mercator[3],
            "MapLatMin" : extent_lat_lon[0], 
            "MapLatMax" : extent_lat_lon[1], 
            "MapLonMin" : extent_lat_lon[2], 
            "MapLonMax" : extent_lat_lon[3],
            "MapPixXMin" : map_pixel_extent[0], 
            "MapPixXMax" : map_pixel_extent[1], 
            "MapPixYMin" : map_pixel_extent[2], 
            "MapPixYMax" : map_pixel_extent[3],
            "SurplusArea" : surplus_area}

    # create string with map information
    description = "------ SIZE ------\n\n"
    description += """Size of resulting image will be (in pixels): 
    width: %d
    height: %d
    yielding a %.2f Megapixel image\n\n""" %(image_pixel_size[0],image_pixel_size[1],
                                        image_pixel_size[0]*image_pixel_size[1]/10**6)

    description += """When printed in a resolution of %.1f DPI, the physical size of this print in cm is:
    width: %.2f
    height: %.2f
    yielding a %.3f m^2 print\n\n""" %(DPI,image_phys_size[0],image_phys_size[1],
                                    image_phys_size[0]*image_phys_size[1]/10**4)
    description += """The scale of the map in meter/pixel is:
    top: %.3f 
    bottom: %.3f\n\n""" %(pixel_scale_top, pixel_scale_bottom)

    description += "-----COVERAGE-------\n\n"
    description += """The actual geographical area covered ranges is between:
    lat: %.6f° - %.6f°
    lon: %.6f° - %.6f°\n""" %extent_lat_lon
    description += "This leads to a map that is %.2f%s bigger than the region of interest" %(100*surplus_area,"%")
    return maptiles, area, image_pixel_size, area_tile_coords, map_coords, data, description



def write_map_csv(MAP_FOLDER, data):
    # try to open the csv in case it already exists
    try: 
        df = pd.read_csv("../maps/maps_metadata.csv", sep=",")
    except:
        # csv doesn't exist, create new dataframe
        df = pd.DataFrame(columns=data.keys())
    df = df.set_index("MapName")

    # check if this map is already in the csv
    if data["MapName"] in df.index.values:
        # drop the row
        df = df.drop(data["MapName"])

    # append data array
    #df = df.append(data, ignore_index=True)
    df.loc[data["MapName"]] = [data[key] for key in list(data.keys())[1:]]

    # write out csv
    df.to_csv("../maps/maps_metadata.csv", sep=",")