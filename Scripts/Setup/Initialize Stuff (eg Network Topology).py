# this script may need to be run once the project has been moved to a new computer, but in practice was not needed

import os

aprx = arcpy.mp.ArcGISProject("CURRENT")
arcpy.AddMessage(aprx.filePath)
m = aprx.activeMap
arcpy.AddMessage(f"Current Map Name: {m.name}")
gdb = os.path.join(aprx.defaultGeodatabase)
arcpy.AddMessage(f'Default geodatabase is set to: {gdb}')

rootfolder = aprx.filepath


TraceNetworkRelDir = r"Data\TraceNetwork\NHD_H_Virginia_State_GDB\NHD_H_Virginia_State_GDB.gdb\Hydrography\NHD_TraceNetwork"
TraceNetworkFullDir = os.path.join(rootfolder,TraceNetworkRelDir)


arcpy.tn.EnableNetworkTopology(TraceNetworkFullDir, 10000, "ENABLE_TOPO")


