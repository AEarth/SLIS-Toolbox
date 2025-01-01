import arcpy
import os, re


import warnings
warnings.filterwarnings('ignore')

aprx = arcpy.mp.ArcGISProject("CURRENT")
arcpy.AddMessage(aprx.filePath)
m = aprx.activeMap
arcpy.AddMessage(f"Current Map Name: {m.name}")

gdb = os.path.join(aprx.defaultGeodatabase)

arcpy.AddMessage(f'Default geodatabase is set to:  {gdb}')

rootfolder = aprx.homeFolder

a = re.search('\d\d\d\d\d\d', rootfolder)

s = a.group(0)



vainvnum = "VA"+s

scratchdir = os.path.join(rootfolder, "scratch.gdb")


arcpy.env.workspace = gdb

def gdbname(suffix):
    return os.path.join(gdb, vainvnum+"_"+suffix)

def scratchname(suffix):
    return os.path.join(scratchdir, suffix)


DEM = gdbname("DEMclip")
CRS = arcpy.Describe(DEM).spatialReference


#CREATE EMPTY INEF FLOW POLYGON
ineflayername = vainvnum+"_InefPoly"
InefFlowPoly = gdbname("InefPoly")
templatepoly = r"Templates\IneffectiveFlow_poly"
arcpy.management.CreateFeatureclass(gdb, ineflayername,"POLYGON",templatepoly, "DISABLED", "DISABLED", CRS)
arcpy.management.AddField(InefFlowPoly, "VALUE", "SHORT")
m.addDataFromPath(InefFlowPoly)
arcpy.management.CalculateField(InefFlowPoly, "Value", "99", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")

#arcpy.management.ApplySymbologyFromLayer(ineflayername, templatepoly, None, "MAINTAIN")


arcpy.AddMessage("#------USER MUST DRAW INEFFECTIVE FLOW AREAS USING InefPoly------#")
arcpy.AddMessage("#---------Then Run Step C - Hydraulic Properties---------#")
