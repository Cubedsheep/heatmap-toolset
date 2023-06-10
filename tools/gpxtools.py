import os
import gpxpy
from gpxpy import gpxxml as gpxml
import pandas as pd

def clean_extensions(g: gpxpy.gpx.GPX):
    g.extensions = []
    for rte in g.routes:
        rte.extensions = []
        for pt in rte.points:
            pt.extensions = []
    for t in g.tracks:
        t.extensions = []
        for s in t.segments:
            s.extensions = []
            for p in s.points:
                p.extensions = []

def clean_waypoints(g: gpxpy.gpx.GPX):
    g.waypoints = []
    for rte in g.routes:
        rte.waypoints = []
        for pt in rte.points:
            pt.waypoints = []

# TODO: implement removing speed
def clean_and_merge(filenames, GPX_FOLDER, DEST_FOLDER, SIMPLIFY, REMOVE_TIME, REMOVE_ELEVATION, REMOVE_EXTENSIONS,
    REMOVE_WAYPOINTS, REMOVE_SPEED, MERGE_NAME, steps):
    total_files = len(filenames)
    if steps > 0:
        step = total_files // steps
        step = max((1,step))

    # open all gpxfiles in the GPX_FOLDER, read them, simplify them and save the xml string
    gpx_strings = []

    filenum = 0

    for filename in filenames:
        filenum += 1
        if filenum % step == 0 and steps > 0:
            print("Finnished file %i of %i, %2.0f percent of total files finnished" %(filenum, total_files, 100*filenum/total_files))
        
        if type(filename) is str:
            with open(GPX_FOLDER + filename, "r") as gpx_file:
                gpx = gpxpy.parse(gpx_file)
                if SIMPLIFY:
                    gpx.simplify()
                if REMOVE_TIME:
                    gpx.remove_time()
                if REMOVE_ELEVATION:
                    gpx.remove_elevation()
                if REMOVE_EXTENSIONS:
                    clean_extensions(gpx)
                if REMOVE_WAYPOINTS:
                    clean_waypoints(gpx)
                xml_string = gpx.to_xml()

                gpx_strings.append(xml_string)
            
    # merge simplified gpxs and write to file
    with open(DEST_FOLDER + MERGE_NAME, "w") as merged_gpxfile:
        print(gpxml.join_gpxs(gpx_strings), file=merged_gpxfile)


def convert_to_gpx(SOURCE_FOLDER, DEST_FOLDER, CONVERT_METADATA):

    EXTENSIONS = [".gpx", ".GPX", ".tcx", ".fit"]

    # unzip all zips in the SOURCE_FOLDER
    command = "gunzip " + SOURCE_FOLDER + "*.gz"
    os.system(command)

    # iterate over all files in SOURCE_DIR, convert the relevant
    # files and move them to DEST_FOLDER
    for name in os.listdir(SOURCE_FOLDER):
        extension = name[-4:]
        # test if its a file we are interested in
        if extension in EXTENSIONS:
            # just move gpxfile
            if extension == ".gpx":
                command = "mv %s%s %s%s" %(SOURCE_FOLDER, name, DEST_FOLDER, name)
                os.system(command)
            elif extension == ".GPX":
                command = "mv %s%s %s%s.gpx" %(SOURCE_FOLDER, name, DEST_FOLDER, name[:-4])
                os.system(command)
            elif extension == ".tcx":
                # cleanup fucking mess in first line of xml
                with open("%s%s"%(SOURCE_FOLDER,name), "r") as tcxin, open("%s%s"%(DEST_FOLDER,name), "w") as tcxout:
                    content = tcxin.read()
                    content = content.replace("          ", "")
                    tcxout.write(content)

                # convert xml to gpx and remove xml file
                command = """gpsbabel -i gtrnctr -f "%s%s" -o gpx -F "%s%s.gpx" """ %(DEST_FOLDER, name, DEST_FOLDER, name[:-4])
                os.system(command)
                os.system("rm %s%s"%(SOURCE_FOLDER,name))
                os.system("rm %s%s"%(DEST_FOLDER,name))
            elif extension == ".fit":
                command = """gpsbabel -i garmin_fit -f "%s%s" -o gpx -F "%s%s.gpx" """ %(SOURCE_FOLDER, name, DEST_FOLDER, name[:-4])
                os.system(command)
                remove = "rm %s%s"%(SOURCE_FOLDER,name)
                os.system(remove)
            else: print("Handling of file extension '%s' not implemented"%extension)

    if CONVERT_METADATA:
        # read in the strava metadata
        df = pd.read_csv(SOURCE_FOLDER + "metadata.csv", index_col="Unnamed: 0")
        # iterate over dataframe and update filenames
        for i in range(len(df)):
            file_name = df.iloc[i]["Filename"]
            if type(file_name) is str:
                df.iloc[i]["Filename"] = file_name.split("/")[1].split(".")[0] + ".gpx"

    df.to_csv(DEST_FOLDER + "metadata.csv")
