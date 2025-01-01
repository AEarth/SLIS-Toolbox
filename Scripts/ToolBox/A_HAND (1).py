import os, arcpy, re
from arcpy.sa import *


DamPoints = arcpy.GetParameterAsText(0)
s = arcpy.GetParameterAsText(1)
rootfolder = arcpy.GetParameterAsText(2)
DEM = arcpy.GetParameterAsText(3)
burnlen = arcpy.GetParameterAsText(4)


vainvnum = "VA" + s
arcpy.AddMessage(f"#-----Begin HAND process on {vainvnum}-----#")

arcpy.env.overwriteOutput = True

folder = os.path.join(rootfolder, s)
gdb = os.path.join(folder, s+".gdb")
scratchdir = os.path.join(folder, "scratch.gdb")


aprx = arcpy.mp.ArcGISProject("CURRENT")
arcpy.AddMessage(aprx.filePath)
m = aprx.activeMap

def insert_layer_data_toSymbol(TemplateName, LayerName):
    lyr = m.listLayers(TemplateName)[0]
    cp = lyr.connectionProperties
    cp['connection_info']['database'] = gdb
    cp['dataset'] = LayerName
    lyr.updateConnectionProperties(lyr.connectionProperties, cp)

def gdbname(suffix):
    return os.path.join(gdb, vainvnum+"_"+suffix)

def scratchname(suffix):
    return os.path.join(scratchdir, suffix)

os.makedirs(folder, exist_ok=True)

if os.path.exists(gdb) == False:
    arcpy.management.CreateFileGDB(folder, s+".gdb")

if os.path.exists(scratchdir) == False:
    arcpy.management.CreateFileGDB(folder, "scratch.gdb")

aprx = arcpy.mp.ArcGISProject("CURRENT")
m = aprx.activeMap
arcpy.AddMessage(f"Current Map Name: {m.name}")

arcpy.env.workspace = gdb
arcpy.env.scratchWorkspace = scratchdir
aprx.defaultGeodatabase = gdb
arcpy.AddMessage(f'Default geodatabase is set to:  {aprx.defaultGeodatabase}')
aprx.homeFolder = folder

arcpy.AddMessage('Home folder is set to: ' + aprx.homeFolder)

keywords = ["temp", "CulvPoints","startpoint","dpoint","_NHDTRACE"]
for lyr in m.listLayers():
    lyrname = str(lyr)
    if re.search(r"VA\d\d\d\d\d\d", lyrname):
        arcpy.AddMessage("Removing Layer: " + lyrname)
        m.removeLayer(lyr)
    for word in keywords:
        if re.search(word, lyrname):
            arcpy.AddMessage("Removing Layer: "+ lyrname)
            m.removeLayer(lyr)
    for lyr in m.listTables():
        lyrname = str(lyr)
        if re.search(r"VA\d\d\d\d\d\d", lyrname):
            arcpy.AddMessage("Removing Table: " + lyrname)
            m.removeTable(lyr)
        if re.search("MinTable", lyrname):
            m.removeTable(lyr)


#m.spatialReference = arcpy.Describe(DEM).spatialReference

# global DEM vars #
CRS = arcpy.Describe(DEM).spatialReference
griddesc = arcpy.management.GetRasterProperties(DEM, "CELLSIZEX")
gridsize = float(griddesc.getOutput(0))
gridsize2 = 4*gridsize
arcpy.AddMessage(f"Resample Gridsize: {gridsize2} ft" )

arcpy.env.workspace = scratchdir
arcpy.env.overwriteOutput = 'True'


#ADJUST DAM POINT IF NEEDED TO SAMPLE MORE CONSERVATIVE ELEVATION FOR SLOPE CALC

## PROCEDURE BEGINS HERE##
def run(burnlen):
    #for s in ToDoList:
    dam(s)
    nhdtrace(s)
    domainextent(s)
    #extentlist = domainextent.extentlist
    burnculverts(s, burnlen)
    if bool(burnculverts.DEMburn):
            DEMclip = burnculverts.DEMburn
    else:
        DEMclip = domainextent.DEMclip          
    fillflowd8(s, DEMclip)
    flowpathtrace(s, fillflowd8.DEMfill)
    handflood(s, dam.Height)
    # symbology(s)
    # cleanup(s)

#FUNCTION DEFINITIONS BEGIN HERE
def dam(invnum):
    arcpy.AddMessage("#------Checking Data on Dam from Dam Point Table------#")
    with arcpy.da.SearchCursor(DamPoints,['IdNumber', 'MaxV', 'TopH', 'Name'] ) as cursor:
        for row in cursor:
            if row[0] == invnum:
                dam.Id = row[0]
                dam.Volume = row[1]
                dam.Height = row[2]
                dam.Name = row[3]
                #dam.County = row[4]
                arcpy.AddMessage(f"{dam.Id} - {dam.Name} || Dam Height: {round(dam.Height,2)} || Dam Volume: {round(dam.Volume,2)}")


#Create & manually move dam point for trace
def dampoint(invnum):
    arcpy.AddMessage("#------Creating individual dam point feature------#")
    arcpy.management.SelectLayerByAttribute(DamPoints, "NEW_SELECTION", f"IdNumber = '{invnum}'")
    arcpy.conversion.FeatureClassToFeatureClass(DamPoints,cartogdb, vaid+"_DamPoint")



def nhdtrace(invnum):
     ## ARCGIS PRO TRACENETWORK STEPS ##
    arcpy.AddMessage("#------Performing NHD Downstream Trace------#")
    selectphrase = "IdNumber = '{}'".format(invnum)
    arcpy.management.SelectLayerByAttribute(DamPoints, "NEW_SELECTION", selectphrase, None)
    arcpy.conversion.FeatureClassToFeatureClass(DamPoints, scratchdir, "dpoint")
    nhdtrace.dpoint = os.path.join(scratchdir, "dpoint")
    CRSnhd = arcpy.Describe(r"NHD\NHDFlowline").spatialReference
    #arcpy.management.Project("dpoint", "startpoint" CRSnhd, "NO_PRESERVE_SHAPE", None, "NO_VERTICAL")

    #flow length regression equation
    #arcpy.AddMessage(f"Dam Height: {dam.Height} Dam Volume: {dam.Volume}")
    flowlen = round((dam.Volume*0.0007 + dam.Height*0.2818 -1.06),2)
    arcpy.AddMessage(f"predicted study length: {flowlen} miles")

    # start trace point on nearest flowline  ##
    arcpy.analysis.Near(nhdtrace.dpoint, r"NHD\NHDFlowline", None, "LOCATION", "NO_ANGLE", "PLANAR", "NEAR_FID NEAR_FID;NEAR_DIST NEAR_DIST;NEAR_X NEAR_X;NEAR_Y NEAR_Y")

    arcpy.management.MakeXYEventLayer(nhdtrace.dpoint, "NEAR_X", "NEAR_Y", "startpoint", CRSnhd, None)

    ## downstream trace by distance##

    length = str(flowlen) + " miles"

    #NHD is in decimal degrees therefor conversion from miles to DD required (look into reprojecting NHD another time)
    DDdist = flowlen * (2*math.pi*3959*math.cos(37.5)/360)**-1
    clause = "ADD 'Shape length' IS_GREATER_THAN {d} true".format(d = DDdist)

    EmptyBarrierLayer = r"Templates\EmptyBarrierLayer"

    arcpy.tn.Trace("NHD_TraceNetwork Trace Network", "DOWNSTREAM", "startpoint", EmptyBarrierLayer, "NO_DIRECTION", '', "EXCLUDE_BARRIERS", "DO_NOT_VALIDATE_CONSISTENCY", "DO_NOT_IGNORE_BARRIERS_AT_STARTING_POINTS", "IGNORE_INDETERMINATE_FLOW", None, clause, "BOTH_JUNCTIONS_AND_EDGES", None, None, "SELECTION", "NEW_SELECTION", "CLEAR_ALL_PREVIOUS_TRACE_RESULTS", '', "Trace_Results_Aggregated_Points", "Trace_Results_Aggregated_Lines", None, "DO_NOT_USE_TRACE_CONFIGURATION", '', None)
        ## copy / save trace ##

    outname = "t" + invnum + "_NHDTRACE"
    nhdtrace.outtrace = os.path.join(gdb, outname)

    CRS = arcpy.Describe(DEM).spatialReference
    ## remove z value export trace ##
    with arcpy.EnvManager(outputZFlag="Disabled", outputMFlag="Disabled", outputCoordinateSystem = CRS):
          arcpy.management.Dissolve(r"NHD\NHDFlowline", nhdtrace.outtrace, None, None, "MULTI_PART", "DISSOLVE_LINES")
    arcpy.management.SelectLayerByAttribute(r"NHD\NHDFlowline", "CLEAR_SELECTION", '', None)
    m.addDataFromPath(nhdtrace.outtrace)


def domainextent(invnum, buffer="1 Miles"):
    arcpy.AddMessage("#------Setting computational / DEM domain------#")
    vainvnum = "VA" + invnum
    trace = nhdtrace.outtrace

    ## bounding box on downstream trace ##
    arcpy.env.addOutputsToMap = False
    firstbound = os.path.join(scratchdir, "boundbox1")
    with arcpy.EnvManager(outputCoordinateSystem = CRS):
        arcpy.management.MinimumBoundingGeometry(trace, firstbound, "ENVELOPE", "NONE", None, "NO_MBG_FIELDS")

    ## buffer bounding box by 1 mi ##
    lastbound = os.path.join(scratchdir, "biggerbox")
    arcpy.analysis.Buffer(firstbound, lastbound, "1 Miles", "FULL", "FLAT", "ALL", None, "PLANAR")

    ## gridsize2 & clip dem ##
    arcpy.env.addOutputsToMap = True
    DEMclip = os.path.join(gdb, vainvnum + "_DEMclip")
    arcpy.AddMessage(DEMclip)
    desc = arcpy.Describe(lastbound)
    xmin = desc.extent.XMin
    xmax = desc.extent.XMax
    ymin = desc.extent.YMin
    ymax = desc.extent.YMax
    domainextent.extentlist = "{} {} {} {}".format(xmin, ymin, xmax, ymax)
    arcpy.AddMessage("env extent:" + domainextent.extentlist)

    with arcpy.EnvManager(extent=domainextent.extentlist):
        arcpy.management.Resample(DEM, DEMclip, gridsize2, "BILINEAR")
    domainextent.DEMclip = DEMclip


def burnculverts(invnum, burnlen):
    arcpy.AddMessage("#------Hydroenforcing DEM with culverts------#")
    #AUTOMATED BURN MIN TABLE FROM CULVERT TRANSECTS
    vainvnum = "VA"+invnum
    DEMclip = domainextent.DEMclip
    arcpy.env.workspace = scratchdir
    trace = nhdtrace.outtrace

    arcpy.management.SelectLayerByLocation(r"Civil\Roads_Rails", "INTERSECT", trace, "0 Feet", "NEW_SELECTION", "NOT_INVERT")
    # check if any selected
    desc = arcpy.Describe(r"Civil\Roads_Rails")
    string = desc.FIDSet
    #IF TRUE BURN CULVERTS
    if bool(string) == True:
        idlist = string.split(";")
        count = len(idlist)
        arcpy.AddMessage(f"Culvert count to burn: {count}")
        #BURN CULVERTS ON ROADS & RAILROADS
        trace = os.path.basename(nhdtrace.outtrace)
        arcpy.analysis.Intersect(f"{trace} #;Civil\\Roads_Rails #", "CulvPoints", "ONLY_FID", None, "POINT")
        #buffer point
        arcpy.analysis.Buffer("CulvPoints", "CulvBuf", f"{burnlen} Feet", "FULL", "ROUND", "NONE", None, "PLANAR")
        #burnlines lines
        arcpy.analysis.Intersect(f"{trace} #;CulvBuf #", "scratchlines", "ONLY_FID", None, "LINE")
        #dissolve to solve overlapping (two lane roads)
        arcpy.management.Dissolve("scratchlines", "burnlines", None, None, "SINGLE_PART", "DISSOLVE_LINES")
        #for sampling
        arcpy.analysis.Buffer("BurnLines", "samplebuf", "25 Feet", "FULL", "ROUND", "NONE", None, "PLANAR")
        #25% zonal min 
        arcpy.ia.ZonalStatisticsAsTable("samplebuf", "ORIG_FID", domainextent.DEMclip, f"VA{s}_MinTable", "DATA", "PERCENTILE", "CURRENT_SLICE", [25], "AUTO_DETECT")
        arcpy.management.JoinField("BurnLines", "OBJECTID", f"VA{s}_MinTable", "ORIG_FID", "PCT25")
        BurnLineRast = vainvnum+"_BurnLineRast"
        with arcpy.EnvManager(snapRaster=domainextent.DEMclip):
            arcpy.conversion.FeatureToRaster("BurnLines", "PCT25", BurnLineRast)
        BurnLineRastPath = os.path.join(scratchdir, BurnLineRast)
        #m.addDataFromPath(BurnLineRastPath)
        arcpy.env.workspace = gdb
        DEMclip = os.path.basename(domainextent.DEMclip)
        DEMburn = vainvnum+"_DEMburn"
        arcpy.management.MosaicToNewRaster(f"{BurnLineRastPath}; {DEMclip}", gdb, DEMburn , None, "32_BIT_FLOAT", 10, 1, "FIRST", "FIRST")
        burnculverts.DEMburn = DEMburn
        #m.addDataFromPath(os.path.join(gdb, DEMburn))
        arcpy.management.SelectLayerByAttribute(r"Civil\Roads_Rails", "CLEAR_SELECTION", '', None)

    else:
        arcpy.AddMessage("No Culverts To Burn, Will Perform Fill/D8 on Original DEM Clip")
        burnculverts.DEMburn = None

def fillflowd8(invnum, DEMclip):
    arcpy.AddMessage("#------Calculating D8 and Fill------#")
    #fill to create hydro dem#
    DEMfill = os.path.join(gdb, "VA"+invnum + "_DEMfill")
    Filled = arcpy.sa.Fill(DEMclip, None) 
    Filled.save(DEMfill)
    #D8 Flow Dir#
    FlowD8 = arcpy.sa.FlowDirection(DEMfill, "NORMAL", None, "D8")
    FlowD8.save(os.path.join(scratchdir,"FlowD8"))
    fillflowd8.DEMfill = DEMfill
    #m.addDataFromPath(DEMfill)

def flowpathtrace(invnum, DEMfill):
    arcpy.AddMessage("#------Creating new Flow Path based on Hydrocorrected DEM------#")
    #ArcHydro Flow Path
    import flowpathtracing
    flowd8 = os.path.join(scratchdir, "FlowD8")
    p_fpt = flowpathtracing.FlowPathTracing()
    runflow = p_fpt.executeAH(nhdtrace.dpoint, flowd8)
    runflow


    #save final streamline layer
    FlowPath = scratchname("FlowPath")
    streamline = "VA"+ invnum + "_streamline"
    arcpy.conversion.FeatureClassToFeatureClass(FlowPath, gdb, streamline)

    #convert line to raster (maintain extent and grid size/ snap )
    StreamRaster = os.path.join(gdb, "VA"+invnum+ "_StreamRaster")
    with arcpy.EnvManager(snapRaster=DEMfill, extent=domainextent.extentlist):
        arcpy.conversion.FeatureToRaster(streamline, "OBJECTID", StreamRaster, gridsize2)
    flowpathtrace.StreamRaster = StreamRaster

def handflood(invnum, Height):
    arcpy.AddMessage("#----Calculating Relative Height Flood Extent and Depth Raster----#")
    vainvnum = "VA" + invnum
    DEMfill = fillflowd8.DEMfill
    StreamRaster = flowpathtrace.StreamRaster
    flowd8 = os.path.join(scratchdir, "FlowD8")
    #Create HAND / built-in Vertical Flow distance function
    arcpy.CheckOutExtension("Spatial")
    HANDraster = os.path.join(gdb, vainvnum+"_HAND")
    HAND = arcpy.sa.FlowDistance(StreamRaster, DEMfill, flowd8, "VERTICAL", "D8", "MINIMUM")
    HAND.save(HANDraster)
    # Pull dam height from feature and calc max depth
    depth = Height*.6
    roundepth = round(depth,2)
    arcpy.AddMessage(f"Dam Height: {round(Height,2)}; Eval. Depth: {round(depth,2)}")
    #slice handraster to eval height and calc max depths
    floodraster = os.path.join(gdb, vainvnum + "_depthraster_HAND")
    inRas = Raster(HANDraster)
    FloodDepths = (Con(inRas < depth, inRas, None) - depth) * -1
    FloodDepths.save(floodraster)
    m.addDataFromPath(floodraster)
    #output extents
    rawextent = os.path.join(gdb, vainvnum + "_fldextntraw_HAND")
    with arcpy.EnvManager(outputZFlag="Disabled", outputMFlag="Disabled",overwriteOutput="True"):
            arcpy.ddd.RasterDomain(floodraster, rawextent, "POLYGON")
    #smooth output
    tolerance = gridsize2*2
    CRS = arcpy.Describe(DEMfill).spatialReference
    unit = CRS.linearUnitName
    if unit == '':
      unit = CRS.angularUnitName
    if unit == "Foot_US":
        unit = "Feet"
    simplifyphrase = "{} {}".format(tolerance, unit)
    arcpy.AddMessage(f"Simplify Output Extent by {tolerance} {unit}")
    arcpy.AddMessage(simplifyphrase)
    floodextent = os.path.join(gdb, vainvnum + "_floodextent_HAND")
    with arcpy.EnvManager(outputZFlag="Disabled", outputMFlag="Disabled", overwriteOutput="True"):
            arcpy.cartography.SimplifyPolygon(rawextent, floodextent, "BEND_SIMPLIFY", simplifyphrase, "0 SquareMeters", "RESOLVE_ERRORS", "NO_KEEP", None)
    m.addDataFromPath(floodextent)
    handflood.floodextent = floodextent

def cleanup(invnum):
    arcpy.AddMessage("#----Cleaning up layers----#")
    keywords = ["MinTable","NHDTRACE", "Filled", "HAND", "FlowPath", invnum, "FlowD8", "FloodDepths", "dpoint", "startpoint", "biggerbox", "HAND", "Burn", "Culv_", "CulvBuf", "CulvPoints", "pnt", "buf", "burnlines", "scratchlines"]
    #keepwords = ["floodextent", "depthraster"]
    for lyr in m.listLayers():
        lyrname = str(lyr)
        for word in keywords:
            if re.search(word, lyrname) and not re.search("floodextent",lyrname): #and not re.search("depthraster",lyrname):
                arcpy.AddMessage(f"Removing Layer:  {lyrname}")
                m.removeLayer(lyr)


# def symbology(invnum):
#     arcpy.AddMessage("#----Setting flood extent symbology----#")
#     floodextent = f"VA{invnum}_floodextent_HAND"
#     arcpy.management.ApplySymbologyFromLayer(floodextent, r"Templates\Flood Extent", None, "DEFAULT")


run(burnlen)

# arcpy.AddMessage("#------Setting flood extent symbology------#")
# floodextent = f"VA{invnum}_floodextent_HAND"

# insert_layer_data_toSymbol("Flood Extent", floodextent)
#arcpy.management.ApplySymbologyFromLayer(floodextent, r"Templates\Flood Extent", None, "DEFAULT")

# cleanup(s)
