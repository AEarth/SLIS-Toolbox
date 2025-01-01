import os, re, arcpy


def ScriptTool(DamPoints, GridIndex, DamPoint, StreamlineRTE, RoadX, Structures, PanelGridIndex, PanelFloodExtent):
    aprx = arcpy.mp.ArcGISProject("CURRENT")
    rootfolder = aprx.homeFolder

    arcpy.AddMessage(f"Project: {aprx.filePath}, HomeFolder: {rootfolder}")

    m = aprx.activeMap
    arcpy.AddMessage(f"Current Map Name: {m.name}")

    gdb = os.path.join(aprx.defaultGeodatabase)

    scratchgdb = os.path.join(rootfolder, "scratch.gdb")

    rootpath = aprxfolder.split('GIS Model')[0]

    MainModelPath = os.path.join(rootpath, 'GIS Model')

    SymbologyFolder = os.path.join(MainModelPath, 'Symbology Templates')

    arcpy.AddMessage("Symbology Dir:" + SymbologyFolder)

    a = re.search('\d\d\d\d\d\d', rootfolder)
    s = a.group(0)

    vaid = "VA" + s
    arcpy.AddMessage(f"Current Dam: {vaid}")

    def gdbname(suffix):
        return os.path.join(gdb, vaid + "_" + suffix)

    def tempgdbname(filename):
        return os.path.join(scratchgdb, filename)

    # INPUT PARAMS
    DamPoints = arcpy.GetParameterAsText(0)

    # OUTPUT PARAMS
    OutParamList = ["UnusedSpacer", "GridIndex", "DamPoint", "StreamlineRTE", "RoadX", "Structures", "PanelGridIndex",
                    "PanelFloodExtent"]

    for i in range(1, len(OutParamList)):
        globals()[OutParamList[i] + "Sym"] = os.path.join(SymbologyFolder, OutParamList[i] + ".lyrx")

        if "Panel" in OutParamList[i]:
            ReuseName = OutParamList[i].replace("Panel", "")
            globals()[OutParamList[i] + "File"] = gdbname(ReuseName)
            globals()[OutParamList[i]] = arcpy.SetParameter(i, gdbname(ReuseName))
        else:
            globals()[OutParamList[i] + "File"] = gdbname(OutParamList[i])
            globals()[OutParamList[i]] = arcpy.SetParameter(i, gdbname(OutParamList[i]))

        # arcpy.AddMessage(globals()[OutParamList[i]+"File"])
        # arcpy.AddMessage(globals()[OutParamList[i]])
        # arcpy.AddMessage(globals()[OutParamList[i]+"Sym"])
        i += 1

    GridIndexFile = gdbname("GridIndex")
    DamPointFile = gdbname("DamPoint")
    StreamlineRTEFile = gdbname("StreamlineRTE")
    RoadXFile = gdbname("RoadX")
    StructuresFile = gdbname("Structures")
    PanelFloodExtentFile = gdbname("FloodExtent")

    arcpy.AddMessage(f"#-----Setting Symbol Definitions for Out Parameters-----#:")

    params = arcpy.GetParameterInfo()

    for param in params:
        arcpy.AddMessage("Name: {}, Type: {}, Value: {}".format(param.name, param.parameterType, param.value))

    for i in range(1, len(params)):
        arcpy.AddMessage(f"{i} - {params[i].name}")
        params[i].symbology = globals()[OutParamList[i] + "Sym"]

    # Set defualt workspace
    arcpy.env.workspace = gdb
    arcpy.env.overwriteOutput = True

    streamline = gdbname("streamline")
    DepthRaster = vaid + "_depthraster"
    FloodExtent = vaid + "_floodextent"

    # Create dam point for map
    arcpy.conversion.FeatureClassToFeatureClass(DamPoints, gdb, vaid + "_DamPoint", where_clause=f"IdNumber = '{s}'")
    # Dam Point Symbology
    dampointlyr = f"VA{s}_DamPoint"
    # arcpy.management.ApplySymbologyFromLayer(dampointlyr, r"Templates\DamPoint", None, "DEFAULT")
    # insert_layer_data_toSymbol("DamPoint", dampointlyr)

    # GRID INDEX
    arcpy.cartography.StripMapIndexFeatures(streamline, GridIndexFile, "USEPAGEUNIT", 4000, "10.478 Inches",
                                            "7.6296 Inches", "HORIZONTAL", 10, 1, "WE_NS")
    m.addDataFromPath(GridIndexFile)

    # arcpy.management.ApplySymbologyFromLayer(gridfeatures, r"Templates\Grid Index Symbology", None, "DEFAULT")
    # insert_layer_data_toSymbol("Grid Index", gridfeatures)

    # LAYOUT MANUAL STEP: update PanelLocatorFrame with display option constraint: "linked map frame center", "(main) Map Frame", "GridIndex Layer"

    # aprx = arcpy.mp.ArcGISProject("CURRENT")
    PanelMap = aprx.listMaps("Panel Locator")[0]
    PanelMap.addDataFromPath(PanelGridIndexFile)
    PanelMap.addDataFromPath(PanelFloodExtentFile)

    # arcpy.AddMessage(f"Project: {aprx.filePath}, HomeFolder: {rootfolder}")

    # def insertpaneldata():
    #     arcpy.AddMessage("#-------Inserting Data to 'Panel Locator' Map Templates------#")
    #     insert_layer_data_toSymbol_MapPanel("Flood Extent Symbol", FloodExtent)
    #     insert_layer_data_toSymbol_MapPanel("Grid Index Symbol", gridfeatures)

    # insertpaneldata()

    # UPDATE DEPTH GRADIENT SYMBOL
    arcpy.AddMessage("#-------Updating Depth Raster Range Labels-------#")

    DepthRasterName = vaid + "_depthraster"
    DepthRasterPath = gdbname("depthraster")

    # arcpy.management.ApplySymbologyFromLayer(DepthRaster, r"Templates\Depth Symbology", None, "DEFAULT")
    # Update Depth Gradient Labels

    # Query max depth from Depth Raster Data
    MaxDepthResult = arcpy.GetRasterProperties_management(DepthRasterPath, "MAXIMUM")
    # Get the elevation standard deviation value from geoprocessing result object
    MaxDepth = MaxDepthResult.getOutput(0)

    # Insert depthraster data into template and update labels
    # insert_layer_data_toSymbol("Depth Raster", DepthRasterName)
    lyrobj = m.listLayers(DepthRasterName)[0]
    sym = lyrobj.symbology
    sym.colorizer.minLabel = " 0 ft"
    sym.colorizer.maxLabel = str(round(float(MaxDepth), 1)) + " ft"
    lyrobj.symbology = sym

    # CREATE MEASURES FROM STREAMLINE
    arcpy.AddMessage("#------Creating Measures from Streamline------#")
    Streamline = gdbname("streamline")

    arcpy.management.AddField(Streamline, "MileFrom", "SHORT")
    arcpy.management.AddField(Streamline, "MileEnd", "FLOAT")
    arcpy.management.CalculateField(Streamline, "MileFrom", "0", "PYTHON3", '', "SHORT")
    arcpy.management.CalculateField(Streamline, "MileEnd", "!Shape_Length!", "PYTHON3", '', "FLOAT")
    arcpy.lr.CreateRoutes(Streamline, "DestID", StreamlineRTEFile, "TWO_FIELDS", "MileFrom", "MileEnd", "LOWER_RIGHT",
                          0.00018939394, 0, "IGNORE", "INDEX")

    # m.addDataFromPath(StreamlineRTEFile)

    # Insert Streamline Data into Stream RTE template
    # arcpy.management.ApplySymbologyFromLayer(StreamlineRTE, r"Templates\Stream Measure", None, "DEFAULT")
    # insert_layer_data_toSymbol("Stream RTE", StreamlineRTE)

    # STRUCTURES Depths with lyrx symbology
    arcpy.AddMessage("#------Structure Impacts Procedure Start Here------#")
    DepthRaster = vaid + "_depthraster"
    FloodExtent = vaid + "_floodextent"
    buildingfc = vaid + "_Structures"

    StructureRef = r"Civil\MSFT_Buildings"  # r"Civil\Structures"

    arcpy.AddMessage("1. Selecting structures in flood extent..")
    arcpy.management.SelectLayerByLocation(StructureRef, "INTERSECT", FloodExtent, "0 feet", "NEW_SELECTION",
                                           "NOT_INVERT")

    # check if any selection present
    result = arcpy.management.GetCount(r"Civil\MSFT_Buildings")
    struct_ct = int(result[0])


    if struct_ct < 1:
        arcpy.AddMessage("2. No Structure Impacts Detected!")
    else:
        arcpy.conversion.FeatureClassToFeatureClass(StructureRef, gdb, "temp_Bldgs")

        ## Initial spatial join Building and Address (radius 4.5 m = 15 ft) ##
        arcpy.AddMessage("2. Bringing In Address Data..")
        arcpy.analysis.SpatialJoin("temp_Bldgs", r"civil\AddressPoints", buildingfc, "JOIN_ONE_TO_ONE", "KEEP_ALL",
                                   match_option="INTERSECT", search_radius="15 Feet")

        # Depth Calculation
        arcpy.AddMessage("3. Calculating Depth..")
        BldgDepthTable = os.path.join(gdb, vaid + "_BldgDepthTable")
        arcpy.ia.ZonalStatisticsAsTable(buildingfc, "OBJECTID", DepthRaster, BldgDepthTable, "DATA", "MEAN",
                                        "CURRENT_SLICE")
        arcpy.management.JoinField(buildingfc, "OBJECTID", BldgDepthTable, "OBJECTID_1", "AREA;MEAN")

        # Distance from Dam (LRS)
        arcpy.AddMessage("4. Determining Distance from Dam (mi) and from Thalweg (ft)...")
        BldgDistTable = os.path.join(gdb, vaid + "_BldgDistTable")
        arcpy.management.FeatureToPoint(buildingfc, "temp_BldgPoints", "CENTROID")
        arcpy.lr.LocateFeaturesAlongRoutes("temp_BldgPoints", StreamlineRTEFile, "DestID", "5000 Feet", BldgDistTable,
                                           "RID; Point; Dist2Dam", "FIRST", "DISTANCE", "ZERO", "FIELDS", "NO_M_DIRECTION")
        arcpy.management.JoinField(buildingfc, "OBJECTID", BldgDistTable, "ORIG_FID", "Dist2Dam;Distance")
        arcpy.management.CalculateField(buildingfc, "Distance", "abs(!Distance!)", "PYTHON3", '', "FLOAT")

        # arcpy.management.AlterField(buildingfc, "Join_Count", "class", "class")

        # --Retrieve XY Coordinates in DD for structures without address points
        arcpy.AddMessage("5. Adding WGS 1984 Lat Lon DD positions to Structure Attributes...")
        temp_bldgpnts_prj = tempgdbname("temp_BldgPoints_prj")
        WGS_sr = arcpy.SpatialReference("WGS 1984")
        arcpy.management.Project("temp_BldgPoints", temp_bldgpnts_prj, WGS_sr)
        arcpy.management.AddXY(temp_bldgpnts_prj)
        arcpy.management.JoinField(buildingfc, "OBJECTID", temp_bldgpnts_prj, "ORIG_FID", "POINT_Y;POINT_X")

        # layrx symbology update connection
        arcpy.AddMessage("6. Inserting Data into 'Impacted Structures' Layer Symbol Def...")
        # aprx = arcpy.mp.ArcGISProject("CURRENT")
        # mapx = aprx.activeMap
        # lyr = m.listLayers('Impacted Structures')[0]
        # cp = lyr.connectionProperties
        # cp['connection_info']['database'] = gdb
        # cp['dataset'] = buildingfc
        # lyr.updateConnectionProperties(lyr.connectionProperties, cp)

        ## label structure types manually after add field
        arcpy.management.AddField(buildingfc, "Type", "TEXT")

    # m.addDataFromPath(StructuresFile)

    ## ROAD CROSSING STEPS

    arcpy.AddMessage("#------Road Impacts Procedure Start Here------#")
    # Option 1) exclude depth samples from Bridge decks
    tempintroad = "temp_introad"
    # roadfc = vaid+"_RoadX"
    streamline = vaid + "_streamline"
    RoadDepthTable = vaid + "_RoadDepthTble"
    BridgeCulvert = r"Civil\Bridges & Culverts"

    phrase = f"'Civil\Roads_Rails' #; {FloodExtent} #"  # used to be r"Roads etc\RoadsTypes"
    arcpy.analysis.Intersect(phrase, tempintroad, "ALL", None, "LINE")

    # check if any features resulted
    result = arcpy.management.GetCount(tempintroad)
    road_ct = int(result[0])

    if road_ct < 1:
        arcpy.AddMessage("2. No Structure Impacts Detected!")

    else:
        tempbridgeculv = tempgdbname("temp_bridgeculv")

        arcpy.management.SelectLayerByLocation(BridgeCulvert, "INTERSECT", FloodExtent, "", "NEW_SELECTION", "NOT_INVERT")
        arcpy.conversion.FeatureClassToFeatureClass(BridgeCulvert, scratchgdb, "temp_bridgeculv")

        # subtract from roads
        bridgebuffer = tempgdbname("road_erase_buff")
        arcpy.analysis.Buffer(tempbridgeculv, bridgebuffer, "5 Feet", "FULL", "FLAT", "NONE", None, "PLANAR")
        arcpy.analysis.Erase(tempintroad, tempbridgeculv, RoadXFile, None)

        # road symbology
        # arcpy.management.ApplySymbologyFromLayer(roadfc, r"Templates\Impacted Roads", None, "DEFAULT")
        # insert_layer_data_toSymbol("Impacted Roads", roadfc)

        # sample depths (excludes bridges)
        arcpy.ia.ZonalStatisticsAsTable(RoadXFile, "OBJECTID", DepthRaster, RoadDepthTable, "DATA", "MEAN", "CURRENT_SLICE",
                                        90, "AUTO_DETECT")
        arcpy.management.JoinField(RoadXFile, "OBJECTID", RoadDepthTable, "OBJECTID_1", "MEAN")

        # ROAD LINEAR REFERENCING STEPS
        roaddist_table = vaid + "_roaddist_tble"

        temproadpoints = tempgdbname("temp_roadpoints")

        arcpy.management.FeatureToPoint(RoadXFile, temproadpoints, "INSIDE")
        arcpy.lr.LocateFeaturesAlongRoutes(temproadpoints, StreamlineRTEFile, "DestID", "3000 Feet", roaddist_table,
                                           "RID; Point; Dist2Dam", "FIRST", "DISTANCE", "ZERO", "NO_FIELDS", "M_DIRECTON")
        arcpy.management.CalculateField(roaddist_table, "Distance", "abs(!Distance!)", "PYTHON3", '', "FLOAT")
        arcpy.management.JoinField(RoadXFile, "OBJECTID", roaddist_table, "INPUTOID", "Dist2Dam;Distance")

if __name__ == '__main__':
    # Tool parameter accessed with GetParameter or GetParameterAsText
    aprx = arcpy.mp.ArcGISProject("CURRENT")
    aprxfolder = aprx.homeFolder
    m = aprx.activeMap

    a = re.search('\d\d\d\d\d\d', aprxfolder)
    s = a.group(0)

    gdb = os.path.join(aprx.defaultGeodatabase)

    vaid = "VA" + s


    def gdbname(suffix):
        return os.path.join(gdb, vaid + "_" + suffix)


    GridIndexFile = gdbname("GridIndex")
    DamPointFile = gdbname("DamPoint")
    StreamlineRTEFile = gdbname("StreamlineRTE")
    RoadXFile = gdbname("RoadX")
    StructuresFile = gdbname("Structures")
    PanelFloodExtentFile = gdbname("FloodExtent")

    DamPoints = arcpy.GetParameterAsText(0)
    GridIndex = arcpy.SetParameter(1, GridIndexFile)
    DamPoint = arcpy.SetParameter(2, DamPointFile)
    StreamlineRTE = arcpy.SetParameter(3, StreamlineRTEFile)
    RoadX = arcpy.SetParameter(4, RoadXFile)
    Structures = arcpy.SetParameter(5, StructuresFile)
    PanelGridIndex = arcpy.SetParameter(6, GridIndexFile)
    PanelFloodExtent = arcpy.SetParameter(7, PanelFloodExtentFile)

    ScriptTool(DamPoints, GridIndex, DamPoint, StreamlineRTE, RoadX, Structures, PanelGridIndex, PanelFloodExtent)