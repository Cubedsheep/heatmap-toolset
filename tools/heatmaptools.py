import globalmaptiles as gmaptiles
import geopy.distance
from PIL import Image
import numpy as np
import pandas as pd
import gpxpy
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap



def LatLonToMapPixel(lat, lon, zoom, map_pixel_origin, TILESIZE):
    map_coords = gmaptiles.GlobalMercator(tileSize=TILESIZE)

    # convert lon-lat to meters
    xm, ym = map_coords.LatLonToMeters(lat, lon)
    # convert meters to pixel
    px, py = map_coords.MetersToPixels(xm, ym, zoom)
    mpx = px-map_pixel_origin[0]
    #mpy = map_pixel_origin[1]-py   # reverse because lat coordinates are in reverse wrt pseudo mercator
                                   # might have to be fixed for southern hemisphere, idk
    mpy = py-map_pixel_origin[1]
    
    return int(np.floor(mpx)), int(np.floor(mpy))



def generate_pixel_counts(GPX_SOURCE, MAP_DATA):
    # extract required data
    LAT = (MAP_DATA["MapLatMin"], MAP_DATA["MapLatMax"])
    LON = (MAP_DATA["MapLonMin"], MAP_DATA["MapLonMax"])
    IMG_SHAPE = (MAP_DATA["PixHeight"], MAP_DATA["PixWidth"])
    PIX_ORIGIN = (MAP_DATA["MapPixXMin"], MAP_DATA["MapPixYMin"])
    TILESIZE = MAP_DATA["TileSize"]
    ZOOM = MAP_DATA["Zoom"]

    MAP_SCALE = MAP_DATA["PixScaleTop"]

    print("reading in gpx file")
    # read in gpx file to create heatmap from
    with open(GPX_SOURCE, "r") as gpxfile:
        gpx = gpxpy.parse(gpxfile)
        
    # interval to use, use 1/4 of a pixel to miss as little pixels as possible
    INTERVAL = MAP_SCALE/4
    
    # count how many tracks go through each pixel

    tracknum = 0
    tracks_uniform_density = []
    pixel_sets = []

    for track in gpx.tracks:

        # print out some progress info
        tracknum += 1
        if tracknum % 100 == 0:
            print(tracknum)

        # keep the last 3 pixels to prevent "jumping" around and double
        # counting from pauses
        pixel_min_one = None
        pixel_min_two = None
        pixel_min_three = None
        distance = 0.0

        # iterate over the track segments and construct a uniform density version
        for segment in track.segments:
            # arrays to keep a track with uniformly spaced points
            # and the pixels it crosses
            unidens_track = []
            pixels = []

            # iterate over the points in this segment, skipping the first
            for i in range(1,len(segment.points)):
                # extract coordinates and delta
                coord1 = (segment.points[i-1].latitude, segment.points[i-1].longitude)
                coord2 = (segment.points[i].latitude, segment.points[i].longitude)
                d_lat = coord2[0]-coord1[0]
                d_lon = coord2[1]-coord1[1]
                # segment length in meter
                segment_length = geopy.distance.geodesic(coord1, coord2).m

                # add interpolated points to the track as long as distance smaller than segment_length
                while distance < segment_length:
                    factor = distance/segment_length
                    point_coords = [coord1[0]+factor*d_lat, coord1[1]+factor*d_lon]
                    unidens_track.append(point_coords)
                    # check if point in bbox, then add
                    if LAT[0] < point_coords[0] and point_coords[0] < LAT[1] and  LON[0] < point_coords[1] and point_coords[1] < LON[1]:
                        pixel = LatLonToMapPixel(point_coords[0], point_coords[1], ZOOM, PIX_ORIGIN, TILESIZE)
                        # check if new pixel was not one of the last pixels and add to list of pixels
                        if not (pixel in [pixel_min_one, pixel_min_two, pixel_min_three]):
                            pixels.append(pixel)
                            pixel_min_three = pixel_min_two
                            pixel_min_two = pixel_min_one
                            pixel_min_one = pixel
                    distance += INTERVAL
                # end of interpolation step

                # update distance for next piece of track
                distance -= segment_length
            # end iteration over segment

            # add uniform density segment and list of pixels to the lists
            tracks_uniform_density.append(unidens_track)
            pixel_sets.append(pixels)
        # end iteration over track
    # end iteration over list of tracks

    return pixel_sets



def generate_point_density(pixel_sets, MAP_DATA, DENOISE, DENOISE_RADIUS, WIDEN, WIDEN_RADIUS):
    # extract required data
    LAT = (MAP_DATA["MapLatMin"], MAP_DATA["MapLatMax"])
    LON = (MAP_DATA["MapLonMin"], MAP_DATA["MapLonMax"])
    IMG_SHAPE = (MAP_DATA["PixHeight"], MAP_DATA["PixWidth"])
    PIX_ORIGIN = (MAP_DATA["MapPixXMin"], MAP_DATA["MapPixYMin"])
    TILESIZE = MAP_DATA["TileSize"]
    ZOOM = MAP_DATA["Zoom"]

    MAP_SCALE = MAP_DATA["PixScaleTop"]

    point_dense = np.zeros(IMG_SHAPE, np.uint16)
    has_data = np.zeros(IMG_SHAPE, np.uint8)

    for pixel_set in pixel_sets:
        for pixel in pixel_set:
            point_dense[pixel[1],pixel[0]] += 1
            has_data[pixel[1],pixel[0]] = 1

    denoised = np.zeros_like(point_dense)

    # get a set of all non-zero pixels and iterate over it
    if DENOISE:
        all_pixels = set()
        all_pixels = all_pixels.union(*pixel_sets)
        for pixel in all_pixels:
            # get the x-size of the mask, be careful for edges
            if pixel[0] < DENOISE_RADIUS:
                minx = 0
                maxx = pixel[0]+DENOISE_RADIUS+1  # because ranges are exclusive
            elif pixel[0] > IMG_SHAPE[1]-DENOISE_RADIUS:
                minx = pixel[0]-DENOISE_RADIUS
                maxx = IMG_SHAPE[1]
            else:
                minx = pixel[0]-DENOISE_RADIUS
                maxx = pixel[0]+DENOISE_RADIUS+1

            # get the y-size of the mask, be careful for edges
            if pixel[1] < DENOISE_RADIUS:
                miny = 0
                maxy = pixel[1]+DENOISE_RADIUS+1  # because ranges are exclusive
            elif pixel[1] > IMG_SHAPE[0]-DENOISE_RADIUS:
                miny = pixel[1]-DENOISE_RADIUS
                maxy = IMG_SHAPE[0]
            else:
                miny = pixel[1]-DENOISE_RADIUS
                maxy = pixel[1]+DENOISE_RADIUS+1

            # set value of pixel to max of neighbourhood
            if not WIDEN:
                denoised[pixel[1],pixel[0]] = np.max(point_dense[miny:maxy,minx:maxx])
            else:
                # get mask to widen
                if pixel[0] < WIDEN_RADIUS:
                    wminx = 0
                    wmaxx = pixel[0]+WIDEN_RADIUS+1  # because ranges are exclusive
                elif pixel[0] > IMG_SHAPE[1]-WIDEN_RADIUS:
                    wminx = pixel[0]-WIDEN_RADIUS
                    wmaxx = IMG_SHAPE[1]
                else:
                    wminx = pixel[0]-WIDEN_RADIUS
                    wmaxx = pixel[0]+WIDEN_RADIUS+1

                # get the y-size of the mask, be careful for edges
                if pixel[1] < WIDEN_RADIUS:
                    wminy = 0
                    wmaxy = pixel[1]+WIDEN_RADIUS+1  # because ranges are exclusive
                elif pixel[1] > IMG_SHAPE[0]-WIDEN_RADIUS:
                    wminy = pixel[1]-WIDEN_RADIUS
                    wmaxy = IMG_SHAPE[0]
                else:
                    wminy = pixel[1]-WIDEN_RADIUS
                    wmaxy = pixel[1]+WIDEN_RADIUS+1
                    
                denoised[wminy:wmaxy,wminx:wmaxx] = np.max(point_dense[miny:maxy,minx:maxx])
            # end if
        # end for
        point_dense = denoised
        del denoised

        return point_dense



def generate_heat_layer(point_dense, LINEAR_MAX, LINEAR_WEIGHT, base_cmap, ALPHA):
    # create histogram
    max_count = np.max(point_dense)
    num_pixels, passages = np.histogram(point_dense, bins=max_count)

    # create inverse cdf based heatmap
    cum_pixels = np.zeros(len(num_pixels[LINEAR_MAX:]),np.uint32)
    cummulative_count = 0
    for i in range(LINEAR_MAX,len(num_pixels)):
        cummulative_count += num_pixels[i]
        cum_pixels[i-LINEAR_MAX] = cummulative_count

    # 1 color for each count
    colors = np.zeros((max_count,4))

    # create the linear segment of the colormap
    # the colors for [0,1/3] are linearly mapped to the first LINEAR_MAX bins
    colors[:LINEAR_MAX,:] = base_cmap(np.linspace(0,LINEAR_WEIGHT,LINEAR_MAX))

    # Cummulative part
    colors[LINEAR_MAX:,:] = base_cmap(LINEAR_WEIGHT + (1-LINEAR_WEIGHT)*cum_pixels/cum_pixels[-1])

    # set alpha
    colors[:,3] = ALPHA

    inv_cdf_color = ListedColormap(colors)
    inv_cdf_color = inv_cdf_color.with_extremes(under=np.array([0.0,0.0,0.0,0.0]))

    heat_layer = inv_cdf_color((point_dense.astype(np.float32)-1)/(max_count-1))

    return heat_layer