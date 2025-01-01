import arcpy, re
from arcpy import *
from arcpy.sa import Con
import arcpy.sa as SA
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from statistics import mean
#import seaborn as sns

import warnings


def ScriptTool(DamPoints, LiDAR, FloodExtent, DepthRaster):

    #Initializing Stuff
    warnings.filterwarnings('ignore')
    rootdir = aprxfolder.split('GIS Model')[0]
    MainModelPath = os.path.join(rootdir,'GIS Model')
    SymbologyFolder = os.path.join(MainModelPath, 'Symbology Templates')

    aprx = arcpy.mp.ArcGISProject("CURRENT")
    arcpy.AddMessage(aprx.filePath)
    m = aprx.activeMap

    arcpy.AddMessage(f"Current Map: {m.name}")

    gdb = os.path.join(aprx.defaultGeodatabase)

    arcpy.AddMessage(f'Default geodatabase is set to:  {gdb}')

    homefolder = aprx.homeFolder

    a = re.search('\d\d\d\d\d\d', homefolder)

    s = a.group(0)

    vainvnum = "VA"+s

    arcpy.AddMessage(f'Dam Inventory:  {vainvnum}')


    scratchdir = os.path.join(homefolder, "scratch.gdb")

    arcpy.env.workspace = gdb

    def gdbname(suffix):
        return os.path.join(gdb, f"{vainvnum}_{suffix}")

    DamPoints = arcpy.GetParameterAsText(0)
    LiDAR = arcpy.GetParameterAsText(1)

    FloodExtentFile =  gdbname("floodextent")
    DepthRasterFile = gdbname("depthraster")

    FloodExtent = arcpy.SetParameter(2, FloodExtentFile)
    DepthRaster = arcpy.SetParameter(3, DepthRasterFile)


    params = arcpy.GetParameterInfo()
    
    arcpy.AddMessage("#------Script Tool Params-------#")
    for param in params:
        arcpy.AddMessage("Name: {}, Type: {}, Value: {}".format(param.name, param.parameterType, param.value))

    FloodExtentSym = os.path.join(SymbologyFolder, "FloodExtent.lyrx")
    DepthRasterSym = os.path.join(SymbologyFolder, "DepthRaster.lyrx")

    arcpy.AddMessage(f"Symbol Def Files: \n {FloodExtentSym} \n {DepthRasterSym}")

    params[2].symbology = FloodExtentSym
    params[3].symbology = DepthRasterSym

    # DamPointFile = gdbname("DamPoint")
    # DamPoint = arcpy.SetParameter(3, DamPointFile)
    # DamPointSym = os.path.join(SymbologyFolder, "Dam Point.lyrx")

    def scratchname(suffix):
        return os.path.join(scratchdir, suffix)

    def insert_layer_data_toSymbol(TemplateName, LayerName):
        lyr = m.listLayers(TemplateName)[0]
        cp = lyr.connectionProperties
        cp['connection_info']['database'] = gdb
        cp['dataset'] = LayerName
        lyr.updateConnectionProperties(lyr.connectionProperties, cp)

    def dam(invnum):
            with arcpy.da.SearchCursor(DamPoints,['IdNumber', 'MaxV', 'TopH', 'Name', 'DA', 'CN', 'TopA','PMP_06'] ) as cursor:
                for row in cursor:
                    if row[0] == invnum:
                        dam.Id = row[0]
                        dam.V = row[1]
                        dam.H = row[2]
                        dam.Name = row[3]
                        dam.DA = row[4]
                        dam.CN = row[5]
                        dam.SA = row[6]
                        dam.PMP = row[7]
                        #dam.County = row[4]
                        arcpy.AddMessage(f"{dam.Id} - {dam.Name} || Dam Height: {round(dam.H,2)} || Dam Volume: {round(dam.V,2)}")

    #assigning some global file path / layer name variables
    DEM = gdbname("DEMclip")
    CRS = arcpy.Describe(DEM).spatialReference

    #new code: subset stream for basins
    streamline = gdbname("streamline")
    demfill = gdbname("DEMfill")
    #FlowD8 = gdbname("FlowD8")
    FlowD8 = os.path.join(scratchdir, "FlowD8")
    basinpoints = gdbname("basinpoints")
    accumraster = gdbname("flowaccum")
    snappour = gdbname("snappour")
    basins = gdbname("basins")


    arcpy.AddMessage("#---------Hydraulic Properties Module Begin---------#")

    dam(s)
    RiverLength = int(0)

    arcpy.AddMessage("#---------Determining how to slice reach segments---------#")

    with arcpy.da.SearchCursor(streamline, 'SHAPE@LENGTH') as cursor:
        for row in cursor:
            RiverLength = RiverLength + row[0]

    RiverLengthmi = round(RiverLength/5280,2)
    arcpy.AddMessage("Measured Reach Length: " + str(RiverLengthmi) + " mi")

    #convert miles to percent
    segmentcount = round((RiverLength/5280),0)
    PercentPoints = 100/segmentcount
    seglength = RiverLength/segmentcount
    arcpy.AddMessage(f"Percent: {round(PercentPoints,2)}% Seg Count: {segmentcount} || seglength: {round(seglength,1)} ft")


    arcpy.AddMessage("#---------Generating SubBasins---------#")
    arcpy.management.GeneratePointsAlongLines(streamline, basinpoints, "PERCENTAGE", None, PercentPoints-0.01) #"END_POINTS" not using
    basinpointslyr = "VA"+s+"_basinpoints"

    #Re-use temp D8 Raster from HAND Process
    FlowD8 = os.path.join(scratchdir, "FlowD8")
    #FlowD8 = arcpy.sa.FlowDirection(demfill, "NORMAL", None, "D8"); 
    #FlowD8.save(FlowD8)

    accum_raster = arcpy.sa.FlowAccumulation(FlowD8, None, "FLOAT", "D8"); accum_raster.save(accumraster)
    snap_raster = arcpy.sa.SnapPourPoint(basinpointslyr, accumraster, 30, "OBJECTID"); snap_raster.save(snappour)
    basin_raster = arcpy.sa.Watershed(FlowD8, snappour, "Value"); basin_raster.save(basins)



    arcpy.AddMessage("#---------Applying Ineffective Flow Areas---------#")
    InefFlowPoly = gdbname("InefPoly")

    #convert to raster after drawing
    Temp_InfctRast = "temp_InefctvRast"
    with arcpy.EnvManager(snapRaster=basins):
        arcpy.conversion.PolygonToRaster(InefFlowPoly, "Value", Temp_InfctRast, "CELL_CENTER", "NONE", 10, "BUILD")
    InefFlowRast = "VA" + s + "_basins_inef"
    arcpy.management.MosaicToNewRaster(f"{Temp_InfctRast};{basins}", gdb, InefFlowRast, CRS, "8_BIT_UNSIGNED", 10, number_of_bands=1, mosaic_method="FIRST")



    #ENSURES USE OF INEFFECTIVE FLOW AREAS FOR HYDRAPROPS
    basins = InefFlowRast


    arcpy.AddMessage("#---------Calculating Average Slope---------#")
    with arcpy.EnvManager(scratchWorkspace=gdb, workspace=gdb):
        arcpy.management.GeneratePointsAlongLines(streamline, "temp_elvpts", "PERCENTAGE", None, 100, "END_POINTS")
    ELVpts = gdbname("elevpts") 

    #ADJUST ELEVATION SAMPLE POINT IF NEEDED

    arcpy.sa.ExtractValuesToPoints("temp_elvpts", LiDAR, ELVpts,"NONE", "VALUE_ONLY")
    elev_list = []
    cursor = arcpy.da.SearchCursor(ELVpts, ['OBJECTID','RASTERVALU'])
    for row in cursor:
        elev_list.append(row[1])

    #strseglen = RiverLength*(PercentPoints/100) #seglength already calced
    riseft = [a-b for a,b in zip(elev_list, elev_list[1:])]
    meanslope = riseft[0]/RiverLength
    arcpy.AddMessage("Average Slope: " + str(round((meanslope*100),2)) + "%")

    ###### END OF CODE TO FIND RIVER SLOPE

    #Insert Hand Model and Stage level
    HAND = gdbname("HAND")
    #HAND = HAND_Raw #SA.ExtractByMask(HAND_Raw, inZoneData)
    maxH = dam.H*0.8
    MAXwaterLevel = maxH #input


    arcpy.AddMessage("#---------Creating Slope Raster for Hydra. Prop. Calcs---------#")
    #OUTSIDE LOOP CALCULATE TOTAL DEM SLOPE RASTER:
    DEMclip = gdbname("DEMclip")
    DEM = Raster(DEMclip) #input

    sloperaster = gdbname("sloperaster")
    arcpy.ddd.Slope(DEMclip, sloperaster, "PERCENT_RISE", 1, "PLANAR", "METER")


    arcpy.AddMessage("#---------Creating Mannings N Value Raster---------#")
    #OUTSIDE LOOP CALCULATE MAIN MANNINGS N RASTER:
    domainextentpoly = scratchname("biggerbox")
    desc = arcpy.Describe(domainextentpoly)
    xmin = desc.extent.XMin
    xmax = desc.extent.XMax
    ymin = desc.extent.YMin
    ymax = desc.extent.YMax
    extentphrase = "{} {} {} {}".format(xmin, ymin, xmax, ymax)

    LCraster = gdbname("LCrast")
    arcpy.management.Clip("LandCover_VA", extentphrase, LCraster, domainextentpoly, "255", "NONE" , "NO_MAINTAIN_EXTENT") #NONE instead of maintain_clipping extent

    nvalraster = gdbname("nvals")
    arcpy.ddd.Reclassify(LCraster, "Value", "11 4;21 1;22 1;31 3;41 14;42 8;51 12;61 4;71 3;81 4;82 4;91 10;92 7;254 254", nvalraster, "DATA")

    #describe raster cell grid size
    descRast = arcpy.Describe(DEM)
    x_cell = descRast.meanCellWidth
    y_cell = descRast.meanCellHeight
    cellArea = x_cell * y_cell


    arcpy.AddMessage("#---------Begin Analyzing Hydraulic Properties---------#")
    #THIS WHERE LOOP BEGINS
    depth = maxH*0.2 #set starting depth
    Depth_step = maxH*0.2 #set depth slice iteration step

    loop = 1
    while loop <= 5:

        arcpy.AddMessage("Analyzing Water Level: " + str(round(depth,1)) + "ft")

        # MY REPLACEMENT VERSION
        intdepth = int(round(depth,0))
        depthshort = f"_{intdepth}"

        floodraster = gdbname("depthraster" + depthshort)
        inRas = Raster(gdbname("HAND"))
        FloodDepths = (Con(inRas < depth, inRas, None) - depth) * -1
        FloodDepths.save(floodraster)

        #SURFACE AREA TABLE: Calculate pixel count and Surface Area using Zonal Statistics as Table
        FloodAreaTable = gdbname("areavol"+depthshort)
        outZSaT = SA.ZonalStatisticsAsTable(basins, "VALUE", floodraster, FloodAreaTable, "DATA", "SUM")
        #multiply sum of depths by cell area / 100 sqft to get volume
        arcpy.management.AddField(FloodAreaTable, "Vol", "FLOAT")
        arcpy.management.CalculateField(FloodAreaTable, "Vol", "!SUM!*100", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")
        #sum of depths no longer needed
        #arcpy.management.DeleteField(FloodAreaTable, "SUM", "DELETE_FIELDS")

        #Calculating channel bed area of inundation zone; eqn 5 from Zheng et al. (2018)

        #janky way to clip surface DEM to depthraster (rast)
        #DEM_clip = (DEM + floodraster) - floodraster
        sloperasterslice = gdbname("sloperaster"+depthshort)
        SlopeRasterExtract = arcpy.sa.ExtractByMask(sloperaster, floodraster)
        #SlopeRasterExtract.save(sloperasterslice)

        #create raster for calculating bed area (account for vertical slopes)
        channelbedraster = gdbname("ChannelBedRast"+depthshort)
        ChannelBedRast = cellArea * SA.SquareRoot(1 + SA.Square((SlopeRasterExtract/100)))
        #ChannelBedRast.save(channelbedraster)

        BedAreaTable = gdbname("ChanBedArea"+depthshort)
        SA.ZonalStatisticsAsTable(basins, "VALUE", ChannelBedRast, BedAreaTable, "DATA", "SUM")
        arcpy.management.AlterField(BedAreaTable, "SUM", "BedArea", "BedArea", "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")

        #mannings n slice calculations
        ManNTable = gdbname("nvals"+depthshort)
        ManNRastSlice = gdbname("ManNRast"+depthshort)
        ManNRastlvl = arcpy.sa.ExtractByMask(nvalraster, floodraster)
        #ManNRastlvl.save(ManNRastSlice)

        arcpy.ia.ZonalStatisticsAsTable(basins, "VALUE", ManNRastlvl, ManNTable, "DATA", "MEAN", "CURRENT_SLICE")
        arcpy.management.AlterField(ManNTable, "MEAN", "ManNint", "ManNint", "DOUBLE", 8, "NULLABLE", "DO_NOT_CLEAR")

        #NEW: JOIN TABLES UP
        arcpy.management.JoinField(FloodAreaTable, "VALUE", BedAreaTable, "VALUE", "BedArea")
        arcpy.management.JoinField(FloodAreaTable, "VALUE", ManNTable, "VALUE", "ManNint")

        FloodStatsTable = FloodAreaTable

        #Add depth specific field (other field calcs after all merged for optimization)
        arcpy.AddField_management(FloodStatsTable, "WaterLvl_H", "DOUBLE", 10, 7)
        
        arcpy.management.CalculateField(FloodStatsTable, "WaterLvl_H", depth, "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")

        # Keep the loop going
        loop += 1
        Push_lvl = Depth_step * (loop)
        depth = round(Push_lvl, 2)
    #LOOP ENDS HERE

    arcpy.AddMessage("Completed Hydraulic Property Analysis!")

    ##append all loop tables into one table

    arcpy.AddMessage("#---------Aggregating Hydra. Properties & Calculating Rating Curves---------#")

    TableList = arcpy.ListTables(vainvnum+"*areavol*")
    arcpy.AddMessage(TableList)
    OrigTable = os.path.join(gdb,TableList[0])

    TableDirs = []
    for table in TableList[1:]:
        TableDirs.append(os.path.join(gdb,table))

    separator = ";"
    AppendPhrase = separator.join(TableDirs)

    arcpy.management.Append(AppendPhrase, OrigTable, "NO_TEST")

    HydraProp = vainvnum+"_HydraProp"
    arcpy.conversion.TableToTable(OrigTable, gdb, HydraProp, where_clause = "VALUE <> 99")

    arcpy.management.DeleteField(HydraProp, "SUM", "DELETE_FIELDS")



    ###Calculate for eqn 7-10 in Zheng et al. (2018)
    arcpy.AddField_management(HydraProp, "Length", "SHORT", 10, 7)
    arcpy.AddField_management(HydraProp, "ChanWidth", "DOUBLE", 10, 7)
    arcpy.AddField_management(HydraProp, "CrossArea", "DOUBLE", 10, 7)
    arcpy.AddField_management(HydraProp, "WetPerm", "DOUBLE", 10, 7)
    arcpy.AddField_management(HydraProp, "HydroRad", "DOUBLE", 10, 7)
    arcpy.AddField_management(HydraProp, "Slope", "DOUBLE", 10, 7)
    arcpy.AddField_management(HydraProp, "Flow_Q", "DOUBLE", 10, 7)

    arcpy.management.CalculateField(HydraProp, "Length", seglength, "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")
    arcpy.management.CalculateField(HydraProp, "ChanWidth", "!AREA!/!Length!", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")
    arcpy.management.CalculateField(HydraProp, "CrossArea", "!Vol!/!Length!", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")
    arcpy.management.CalculateField(HydraProp, "WetPerm", "!BedArea!/!Length!", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")
    arcpy.management.CalculateField(HydraProp, "HydroRad", "!CrossArea!/!WetPerm!", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")
    arcpy.management.CalculateField(HydraProp, "Slope", meanslope, "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")

    #calc flow
    arcpy.management.CalculateField(HydraProp, "Flow_Q", "(1.49 / (!ManNint!/100)) * !CrossArea! * (!HydroRad! ** 0.6667)  * (!Slope! ** 0.5)", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")

    #calc velocities
    arcpy.AddField_management(HydraProp, "Vel", "FLOAT", 10, 7)
    arcpy.management.CalculateField(HydraProp, "Vel", "!Flow_Q!/!CrossArea!", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")


    #Synthetic Rating Curve Charting
    FloodStatsTable = HydraProp

    ###Code to calculate Stage-Discharge Curves for Reach Segments
    arcpy.AddMessage("Creating Figures....")
    GID = []
    Q = [] # Discharge
    H = [] # Stage Heights

    with arcpy.da.SearchCursor(HydraProp, ['VALUE','WaterLvl_H','Flow_Q']) as cursor:
        for row in cursor:
            GID.append(row[0])
            H.append(row[1])
            Q.append(row[2])

    giduniq = np.unique(GID)

    data = {'gid': GID,  'Height': H, 'Flow': Q}
    df = pd.DataFrame(data)

    #initiate matplotlib fig
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    #empty lists for pwr eq constants a, b
    alist = []
    blist = []

    #slice the dataframe by GID
    for gid in giduniq:
        df_gid = df.loc[df['gid'] == gid]
        x = df_gid['Flow'].tolist()
        y = df_gid['Height'].tolist()
        x.sort()
        y.sort()

        #power regression equation
        f1 = lambda xvar,a,b: a*xvar**b
        
        const, pcov = curve_fit(f1, x,  y)
        a1 = const[0] 
        b1 = const[1]
        alist.append(a1)
        blist.append(b1)

        astr = str(round(a1,4))
        bstr = str(round(b1,4))
        yout = a1*x**b1
        StageEq = f"GID{gid}: H = {astr} * Q^({bstr})"
        #FlowEq = f"GID{gid}: Q = (H / {astr}) ^ (1/{bstr})"
        arcpy.AddMessage(StageEq)
        #arcpy.AddMessage(FlowEq)
        # #CALC R^2's
        # residuals = y - yout
        # ss_res = numpy.sum(residuals**2)
        # ss_tot = numpy.sum((y-np.mean(y))**2)
        # r_square = 1 - (ss_res / ss_tot)
        # arcpy.AddMessage("R^2 =" + str(round(r_square,4)))

        plt.scatter(x, y, s =7)
        
        #make curve equations graphically smoother
        xcurve = np.arange(min(x), max(x), 1000)
        ycurve = const[0]*xcurve**const[1]
        plt.plot(xcurve, ycurve, linestyle='solid', label = StageEq)

    #plt.axes.set_xticks(xcurve)
    plt.grid(which='both', axis="both")
    plt.title("River Segment Rating Curves (power fit)")
    plt.xlabel("Discharge / Q (cfs)")
    plt.ylabel("Water Level / H (ft)")

    plt.legend(loc='lower right')
    savedir = os.path.join(homefolder,s+"_RatingCurves.png")
    fig.savefig(savedir)
    plt.show()

    arcpy.AddMessage(f"Figure Saved here: {savedir}")
    #END SETTING CHART UP


    # Froehlich Breach Parameter Equations
    arcpy.AddMessage("#---------Calculating Dam Breach Parameters---------#")

    def ftime08(Vol,Hbr):
        Vol = Vol * 43560 #acft to ft^3
        ftime_s = 63.2 * (Vol / (32.17*(Hbr**2)))**0.5
        ftime_min = round(ftime_s / 60,2)
        #arcpy.AddMessage("Btime: " + str(ftime_min) + " min (Frhlch, 2008)")
        return(ftime_min)

    def Qp_NWS(bavg, Hw, Tf, SA):
        Tf = Tf/60 #convert back to hours
        Kq = 23.4*(SA/bavg) # instantaenous flow reduction factor
        Qp = 3.1 * bavg * Hw**1.5 * (Kq/(Kq + Tf*Hw**0.5))**3
        Qp = round(Qp,1)
        #arcpy.AddMessage(str(Qp) + " cfs (NWS SMPDBK, 1991)")
        return Qp

    def bavgwidth08(Vol, Hbr):
        Vol = Vol * 43560 #convert ac-ft to ft^3
        Km = 1.3 #1.0 for piping
        bavg = 0.27 * Km * Vol **(1/3)
        bavg = round(bavg,1)
        #arcpy.AddMessage("Bavg: " + str(bavg) + " ft" + " (Frhlch, 2008)")
        return round(bavg,2)

    def PMF_Dewberry(PMP, CN, A):
        PMF = ((((80*PMP**2))/(0.75*PMP+(600/CN)-6))+(700*PMP/CN))*(A-0.0017*A**2)
        PMF = round(PMF, 1)
        #arcpy.AddMessage(str(PMF) + " cfs (Dewberry, 2008)")
        return PMF  

    arcpy.AddMessage("----Failure Time---")
    Ft08 = ftime08(dam.V,dam.H)
    arcpy.AddMessage("Btime: " + str(Ft08) + " min (Frhlch, 2008)")

    arcpy.AddMessage("----Breach Average Width---")
    Bavg08 = bavgwidth08(dam.V, dam.H)
    arcpy.AddMessage("Bavg: " + str(Bavg08) + " ft" + " (Frhlch, 2008)")

    #NWS SMPDBK Peak Discharge Equation
    arcpy.AddMessage("----DamBrch Peak Discharge Estimate---")
    QpNWS = Qp_NWS(Bavg08, dam.H, Ft08, dam.SA)
    arcpy.AddMessage(str(round(QpNWS,0)) + " cfs")
    arcpy.AddMessage("----PMF Discharge Estimate---")
    PMF1 = PMF_Dewberry(dam.PMP, dam.CN, dam.DA)
    arcpy.AddMessage(str(round(PMF1,0)) + " cfs")

    arcpy.AddMessage("----Total Peak Discharge (Brch + PMP)---")
    Qptot = QpNWS + PMF1
    arcpy.AddMessage(str(round(Qptot,0)) + " cfs")


    #START CALS FOR H-TABLE
    arcpy.AddMessage("#---------Solving for Normal Depths based on Subreach Rating Curves & Simulated Attenuation---------#")

    def Atten(dist, Qp):
        Qpatt = Qp * ((0.0013*(dist**2)) - (0.0625*dist) + 1)
        #Qpatt = Qp* (-0.0424*dist + 0.9976)
        return round(Qpatt,1)

    def Dist(flow):
        for i in range(len(alist)):
            segmi = seglength/5280
            distmi = (segmi*(i+1)) - (segmi/2)
            distmi = round(distmi,2)
            flowatten = Atten(distmi, flow)
            H = alist[i] * flowatten**blist[i]
            H = round(H,2)
            Distrows.append((i+1, distmi, flowatten, H))
        arcpy.AddMessage(Distrows)


    Distrows = []
    Dist(Qptot)

    #make archydro H lookup table

    Htable = vainvnum+"_Htable"
    arcpy.management.CreateTable(gdb, Htable, None, '', '')

    #ADD FIELDS
    arcpy.AddField_management(Htable, "COMID", "SHORT", 10, 7)
    arcpy.AddField_management(Htable, "distmi", "FLOAT", 10, 7)
    arcpy.AddField_management(Htable, "Q", "FLOAT", 10, 7)
    arcpy.AddField_management(Htable, "H", "FLOAT", 10, 7)

    # Open an InsertCursor
    arcpy.AddMessage("#---------Setting up Normal Depths Table---------#")
    cursor = arcpy.da.InsertCursor(Htable,['COMID', 'distmi', 'Q', 'H'])
    for row in Distrows:
        cursor.insertRow(row)
    del cursor

    #ARCHYDRO DEFINE HAND FLOOD FROM TABLE 
    arcpy.AddMessage("#---------Drawing Results using HAND---------#")

    import defineHANDbasedflooddepthandextentfromtable

    arcpy.env.workspace = gdb

    #outdepthrast = vainvnum+"_depthraster"
    #outfloodextent = vainvnum+"_floodextent"
    basins = gdbname("basins")
    temp_floodextent = scratchname("temp_floodextent")

    pProcessor = defineHANDbasedflooddepthandextentfromtable.ClassOp()
    RunFloodTable = pProcessor.executeAH(HAND,basins,Htable, DepthRasterFile, temp_floodextent, scratch_workspace = "", scratch_folder="", leading_space="")
    RunFloodTable

    #dslv smooth symbolize
    arcpy.management.Dissolve(temp_floodextent, "temp_dslv", None, None, "MULTI_PART", "DISSOLVE_LINES")
    arcpy.cartography.SimplifyPolygon("temp_dslv", FloodExtentFile, "BEND_SIMPLIFY", "20 Feet", "0 SquareFeet", "RESOLVE_ERRORS", "NO_KEEP", None)

    #arcpy.management.ApplySymbologyFromLayer(outfloodextent, r"Templates\Flood Extent", None, "MAINTAIN")
    #insert_layer_data_toSymbol("Flood Extent", outfloodextent)
    #depth raster symbology
    #arcpy.management.ApplySymbologyFromLayer(outdepthrast, r"Templates\Depth Symbology", None, "DEFAULT")
    #insert_layer_data_toSymbol("Depth Raster", outdepthrast)


    #remove temp files

    # keywords = ["temp_dslv", "temp_floodextent", "flood_poly", "sloperaster","LCrast","nvals","elvpts","basin_raster","snap_raster","basins_inef","accum_raster","StreamRaster", "temp"]
    # for lyr in m.listLayers():
    #     lyrname = str(lyr)
    #     for word in keywords:
    #         if re.search(word, lyrname):
    #             arcpy.AddMessage("Removing Layer: ", lyrname)
    #             m.removeLayer(lyr)



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
        return os.path.join(gdb, vaid+"_"+suffix)

    DamPoints = arcpy.GetParameterAsText(0)
    LiDAR = arcpy.GetParameter(1)
    FloodExtent = arcpy.SetParameter(2, gdbname("floodextent"))
    DepthRaster = arcpy.SetParameter(3, gdbname("depthraster"))


    ScriptTool(DamPoints, LiDAR, FloodExtent, DepthRaster)
