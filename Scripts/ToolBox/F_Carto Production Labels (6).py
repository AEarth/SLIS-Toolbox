#--------BEFORE CONTINUING RUN VELOCITY AND TTP SCRIPT----------#
import re, os


## LAYOUT CODE (UPDATE SURROUND)
DamPoints = arcpy.GetParameterAsText(0)

aprx = arcpy.mp.ArcGISProject("CURRENT")
arcpy.AddMessage(aprx.filePath)
m = aprx.activeMap
arcpy.AddMessage(f"Current Map Name: {m.name}")
gdb = os.path.join(aprx.defaultGeodatabase)
arcpy.AddMessage(f'Default geodatabase is set to: {gdb}')

def gdbname(suffix):
    return os.path.join(gdb, vaid+"_"+suffix)

rootfolder = aprx.homeFolder

a = re.search('\d\d\d\d\d\d', rootfolder)
s = a.group(0)
vaid = "VA" + s

def dam(invnum):
        with arcpy.da.SearchCursor(DamPoints,['IdNumber', 'Name', "County"] ) as cursor:
            for row in cursor:
                if row[0] == invnum:
                    dam.Id = row[0]
                    dam.Name = row[1]
                    dam.County = row[2]
                    arcpy.AddMessage(f"{dam.Id} - {dam.Name}|| County: {dam.County}")


dam(s)

#Update Map Surround Text
arcpy.AddMessage("#------Updating Map Surround Text--------#")
p = arcpy.mp.ArcGISProject("Current")
for lyt in p.listLayouts():
    if lyt == p.listLayouts('SIMS LAYOUT')[0]:
        for elm in lyt.listElements("TEXT_ELEMENT"): 
                if elm.name == "Main Name Title":
                    nl = '\n'
                    newtitle = f"DCR Dam Break Inundation Study{nl}{dam.Name} (VA-{dam.Id})"
                    arcpy.AddMessage("new title:"+ nl + newtitle)
                    elm.text = newtitle
                elif elm.name == "County Title":
                    nl = '\n'
                    newtitle = f"{dam.County}, VA"
                    arcpy.AddMessage("new county:"+ nl + newtitle.upper())
                    elm.text = newtitle.upper() 
                elif elm.name == "Inventory Text":
                        newtitle = f"{dam.Id}"
                        arcpy.AddMessage("new id:" + nl + newtitle)
                        elm.text = newtitle

#Update Depth Gradient Legend Label Text
arcpy.AddMessage("#------Updating Depth Gradient Legend Labels--------#")
aprx = arcpy.mp.ArcGISProject("CURRENT")
SIMSmap = "SIMS MAP"
mapx = aprx.listMaps(SIMSmap)[0]
arcpy.AddMessage(m.name)
DepthRaster = f"VA{s}_depthraster"
lyrobj = m.listLayers(DepthRaster)[0]
sym = lyrobj.symbology
#maxdepthtext = str(round(float(sym.colorizer.maxLabel),2))
maxdepthtext = sym.colorizer.maxLabel
p = arcpy.mp.ArcGISProject("Current")
for lyt in p.listLayouts():
    for elm in lyt.listElements("TEXT_ELEMENT"): 
        if elm.name == "maxdepth":
            nl = '\n'
            arcpy.AddMessage("old max depth label:" + nl + elm.text)
            arcpy.AddMessage("new max depth label:" + nl + maxdepthtext)
            elm.text = maxdepthtext


#DROP ROAD IMPACT FIELDS BEFORE EXPORT
fc = gdbname("RoadX")
fields = arcpy.ListFields(fc) 
 
keepFields = ["OBJECTID", "Shape", "STREET_NAME_FULL", "VDOT_RTE_NM", "VDO_RTE_NUMBER", "VDOT_RTE_CATEGORY_CD", "VDOT_TRAFFIC_AADT_NBR","AADT_Text","Route_Name","MEAN", "Dist2Dam", "Distance","TTP", "Shape_Length"]

dropFields = [x.name for x in fields if x.name not in keepFields]
# delete fields
arcpy.DeleteField_management(fc, dropFields)

#Formatting etc on roads fields
arcpy.AddMessage("#------Cleaning up Fields/Formatting on Road Impacts Data for Table Export--------#")
arcpy.management.AddField(fc, "FloodLenft", "SHORT", 6, 2, None, '', "NULLABLE", "NON_REQUIRED", '')
arcpy.management.CalculateField(fc, "FloodLenft", "round(!Shape_Length!*3.28084,-2)", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")
arcpy.management.CalculateField(fc, "Distance", "round(!Distance!,0)", "PYTHON3", '', "FLOAT")
arcpy.management.CalculateField(fc, "Dist2Dam", "round(!Dist2Dam!,2)", "PYTHON3", '', "FLOAT")
arcpy.management.CalculateField(fc, "TTP", "round(!TTP!,2)", "PYTHON3", '', "FLOAT")
arcpy.management.CalculateField(fc, "Mean", "round(!Mean!,1)", "PYTHON3", '', "FLOAT")
arcpy.management.CalculateField(fc, "AADT_Text", '"{:,}".format(!VDOT_TRAFFIC_AADT_NBR!)', "PYTHON3", '', "FLOAT")


#EXPORT IMPACT TABLES AFTER REVIEW (DEV)
outputroadexcel = os.path.join(rootfolder, f"VA{s}_RoadX.xls")
arcpy.conversion.TableToExcel("VA"+s+"_RoadX", outputroadexcel, "NAME", "CODE")


#DROP BUILDING IMPACT FIELDS BEFORE EXPORT
arcpy.AddMessage("#------Cleaning up Fields/Format on Structures Impacts Data for Table Export--------#")
fc = gdbname("Structures")
fields = arcpy.ListFields(fc) 
 
keepFields = ["OBJECTID", "Type", "Shape", "Join_Count", "FULLADDR", "PO_NAME", "ZIP_5", "AREA","MEAN", "Dist2Dam", "Distance", "POINT_X", "POINT_Y", "TTP", "Shape_Length", "Shape_Area"]

dropFields = [x.name for x in fields if x.name not in keepFields]
# delete fields
arcpy.DeleteField_management(fc, dropFields)

#Formatting etc on structures fields
arcpy.management.CalculateField(fc, "Distance", "round(!Distance!,0)", "PYTHON3", '', "FLOAT")
arcpy.management.CalculateField(fc, "Dist2Dam", "round(!Dist2Dam!,2)", "PYTHON3", '', "FLOAT")
arcpy.management.CalculateField(fc, "POINT_X", "round(!POINT_X!,5)", "PYTHON3", '', "FLOAT")
arcpy.management.CalculateField(fc, "POINT_Y", "round(!POINT_Y!,5)", "PYTHON3", '', "FLOAT")
arcpy.management.CalculateField(fc, "TTP", "round(!TTP!,2)", "PYTHON3", '', "FLOAT")
arcpy.management.CalculateField(fc, "Mean", "round(!Mean!,1)", "PYTHON3", '', "FLOAT")
arcpy.management.CalculateField(fc, "FULLADDR", "!FULLADDR!.title()", "PYTHON3", '', "TEXT")
arcpy.management.CalculateField(fc, "PO_NAME", "!PO_NAME!.title()", "PYTHON3", '', "TEXT")

arcpy.management.AddField(fc, "Footprint", "LONG", 6, 2, None, '', "NULLABLE", "NON_REQUIRED", '')
arcpy.management.CalculateField(fc, "Footprint", "round(!Shape_Area!*10.76,-2)", "PYTHON3", '', "TEXT", "NO_ENFORCE_DOMAINS")

outputstructexcel = os.path.join(rootfolder, f"VA{s}_Structures.xls")

arcpy.conversion.TableToExcel(fc, outputstructexcel, "NAME", "CODE")