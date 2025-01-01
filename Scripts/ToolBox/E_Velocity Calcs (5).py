import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from statistics import mean
#import seaborn as plt
import os, re

import warnings


def ScriptTool(DamPoints, Roadcallout):

    warnings.filterwarnings('ignore')

    aprx = arcpy.mp.ArcGISProject("CURRENT")
    arcpy.AddMessage(aprx.filePath)
    m = aprx.activeMap
    arcpy.AddMessage(f"Current Map Name: {m.name}")
    gdb = os.path.join(aprx.defaultGeodatabase)
    arcpy.AddMessage(f'Default geodatabase is set to: {gdb}')
    rootfolder = aprx.homeFolder

    scratchdir = os.path.join(rootfolder, "scratch.gdb")

    a = re.search('\d\d\d\d\d\d', rootfolder)
    s = a.group(0)
    vainvnum = "VA" + s

    rootpath = aprxfolder.split('GIS Model')[0]
    MainModelPath = os.path.join(rootpath,'GIS Model')
    SymbologyFolder = os.path.join(MainModelPath, 'Symbology Templates')

    def gdbname(suffix):
        return os.path.join(gdb, vainvnum+"_"+suffix)

    def scratchname(suffix):
        return os.path.join(scratchdir, suffix)

    RoadcalloutFile = gdbname("Roadcallout")

    DamPoints = arcpy.GetParameterAsText(0)
    Roadcallout = arcpy.SetParameter(1, RoadcalloutFile)

    params = arcpy.GetParameterInfo()

    for param in params:
        arcpy.AddMessage("Name: {}, Type: {}, Value: {}".format(param.name, param.parameterType, param.value))


    RoadcalloutSym = os.path.join(SymbologyFolder, "RoadCallout.lyrx")

    params[1].symbology = RoadcalloutSym


    HydraProp = gdbname("HydraProp")
    streamline = gdbname("streamline")

    #River Length (duplicate code)
    RiverLength = int(0)

    with arcpy.da.SearchCursor(streamline, 'SHAPE@LENGTH') as cursor:
        for row in cursor:
            RiverLength = RiverLength + row[0]

    RiverLengthmi = round(RiverLength/5280,2)
    print("Predicted Reach Length: " + str(RiverLengthmi) + " mi")

    #convert miles to percent
    segmentcount = round((RiverLength/5280),0)
    PercentPoints = 100/segmentcount
    seglength = RiverLength/segmentcount

    #Duplicate code for calculating ftime, qptot, etc
    #START CALS FOR H-TABLE
    
    arcpy.AddMessage("#---------Calculating Reach-Average Velocity Curves---------#")

    ###Code to calculate Stage-Discharge Curves for Reach Segments
    print("Creating Figures....")
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

    # #initiate matplotlib fig
    # fig = plt.figure()
    # ax1 = fig.add_subplot(111)

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
        print(StageEq)
        #print(FlowEq)
        # #CALC R^2's
        # residuals = y - yout
        # ss_res = numpy.sum(residuals**2)
        # ss_tot = numpy.sum((y-np.mean(y))**2)
        # r_square = 1 - (ss_res / ss_tot)
        # print("R^2 =" + str(round(r_square,4)))

        # plt.scatter(x, y, s =7)
        
        # #make curve equations graphically smoother
        # xcurve = np.arange(min(x), max(x), 1000)
        # ycurve = const[0]*xcurve**const[1]
        # plt.plot(xcurve, ycurve, linestyle='solid', label = StageEq)

    # #plt.axes.set_xticks(xcurve)
    # plt.grid(which='both', axis="both")
    # plt.title("River Segment Rating Curves (power fit)")
    # plt.xlabel("Discharge / Q (cfs)")
    # plt.ylabel("Water Level / H (ft)")

    # plt.legend(loc='lower right')




    def Atten(dist, Qp):
        Qpatt = Qp * ((0.0013*(dist**2)) - (0.0625*dist) + 1)
        #Qpatt = Qp* (-0.0424*dist + 0.9976)
        return round(Qpatt,1)

    seglength

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
                        arcpy.AddMessage(f"{dam.Id} - {dam.Name} || Dam Height: {dam.H} || Dam Volume: {dam.V}")
    dam(s)
    # Froehlich Breach Parameter Equations

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

    #arcpy.AddMessage("----Failure Time---")
    Ft08 = ftime08(dam.V,dam.H)
    arcpy.AddMessage(str(round(Ft08,0)) + " min")

    #arcpy.AddMessage("----Breach Average Width---")
    Bavg08 = bavgwidth08(dam.V, dam.H)
    Bavg08

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

    Distrows = []
    Dist(Qptot)


    #START SETTING CHART UP
    FloodStatsTable = HydraProp

    ###Code to calculate Stage-Discharge Curves for Reach Segments
    print("Creating Figures....")
    GID = []
    Q = [] # Discharge
    H = [] # Stage Heights
    BA = []
    Vol = []
    W = []
    XA = []
    P = []
    R = []
    S = []
    Vel = []

    with arcpy.da.SearchCursor(HydraProp, ['VALUE','WaterLvl_H','Flow_Q', 'BedArea','Vol','ChanWidth','CrossArea','WetPerm','HydroRad','Slope', 'Vel']) as cursor:
        for row in cursor:
            GID.append(row[0])
            H.append(row[1])
            Q.append(row[2])
            BA.append(row[3])
            Vol.append(row[4])
            W.append(row[5])
            XA.append(row[6])
            P.append(row[7])
            R.append(row[8])
            S.append(row[9])
            Vel.append(row[10])

    # attempt to extract distances. Distance should be calculated with Hydraulic properties most likely
    # giduniq = np.unique(GID)
    # gidrange = np.ptp(giduniq,axis=0)

    # #convert gid's to miles for slope chart
    # distances = [(strseglen*giduniq[i])/5280 for i in range(gidrange)]

    Qp1 = Qptot

    #COMMENTED OUT DUE TO NOT BEING ABLE TO CLONE PYTHON ENVIRONMENT AND INSTALL SNS
    data = {'gid': GID,  'Height (ft)': H, 'Flow (cfs)': Q, 'BedArea (ft^2)': BA, 'Vol (ft^3)': Vol, 'ChanWidth (ft)': W, 'XS-Area (ft^2)': XA, 'WetPerim (ft)': P, 'HydraRad (ft)': R, 'Slope (%)': [i * 100 for i in S], 'Vel (ft/s)': Vel}
    df = pd.DataFrame(data)

    # fig, axes = plt.subplots(3, 3, sharey=False, figsize=(16,10))
    # sns.set_theme(style="whitegrid", palette="Set2", font_scale = 0.625)
    # j = sns.lineplot(ax=axes[0,0], data=data, x="Flow (cfs)", y="Height (ft)", hue = 'gid', palette="Set2", legend=False, markers=True, dashes=True)
    # sns.lineplot(ax=axes[0,1], data=data, x="BedArea (ft^2)", y="Height (ft)", hue = 'gid', palette="Set2", legend=False, markers=True, dashes=True)
    # sns.lineplot(ax=axes[0,2], data=data, x="ChanWidth (ft)", y="Height (ft)", hue = 'gid', palette="Set2", legend=False, markers=True, dashes=True)
    # sns.lineplot(ax=axes[1,0], data=data, x="XS-Area (ft^2)", y="Height (ft)", hue = 'gid', palette="Set2", legend=False, markers=True, dashes=True)
    # sns.lineplot(ax=axes[1,1], data=data, x="WetPerim (ft)", y="Height (ft)", hue = 'gid', palette="Set2", legend=False, markers=True, dashes=True)
    # sns.lineplot(ax=axes[1,2], data=data, x="Height (ft)", y="HydraRad (ft)", hue = 'gid', palette="Set2", legend=False, markers=True, dashes=True)
    # sns.lineplot(ax=axes[2,0], data=data, x="gid", y="Slope (%)")
    # sns.lineplot(ax=axes[2,1], data=data, x="Height (ft)", y="Vel (ft/s)", hue = 'gid', palette="Set2", legend=False, markers=True, dashes=True)
    # g = sns.lineplot(ax=axes[2,2], data=data, x="Flow (cfs)", y="Vel (ft/s)", hue = 'gid', palette="Set2", legend=True)
    # g.set_xlim(0,Qp1)
    # j.set_xlim(0,Qp1)
    # xlabels = ['{:,.0f}'.format(x) + 'K' for x in g.get_xticks()/1000]
    # xlabels = ['{:,.0f}'.format(x) + 'K' for x in j.get_xticks()/1000]
    # g.set_xticklabels(xlabels)
    # j.set_xticklabels(xlabels)

    # fig.savefig(os.path.join(rootfolder,s+"_HydraProp.png"), dpi=150)
    # plt.show()

    #VELOCITY CURVES
    #initiate matplotlib fig
    fig = plt.figure(figsize = (8,5))
    ax1 = fig.add_subplot(111)

    #empty lists for pwr eq constants a, b
    alist = []
    blist = []

    giduniq = np.unique(GID)

    #slice the dataframe by GID
    for gid in giduniq:
        df_gid = df.loc[df['gid'] == gid]
        x = df_gid['Flow (cfs)'].tolist()
        y = df_gid['Vel (ft/s)'].tolist()
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
        VelocEq = f"GID{gid}: H = {astr} * Q^({bstr})"
        #FlowEq = f"GID{gid}: Q = (H / {astr}) ^ (1/{bstr})"
        print(VelocEq)
        # #CALC R^2's
        residuals = y - yout
        ss_res = numpy.sum(residuals**2)
        ss_tot = numpy.sum((y-np.mean(y))**2)
        r_square = 1 - (ss_res / ss_tot)
        print("R^2 =" + str(round(r_square,4)))

        plt.scatter(x, y, s =7)
        
        #make curve equations graphically smoother
        xcurve = np.arange(min(x), max(x), 1000)
        ycurve = const[0]*xcurve**const[1]
        plt.xlim(0,Qp1)
        plt.plot(xcurve, ycurve, linestyle='solid', label = VelocEq)

    #plt.axes.set_xticks(xcurve)
    plt.grid(which='both', axis="both")
    plt.title("Velocity - Discharge Functions (power fit)")
    plt.xlabel("Discharge / Q (cfs)")
    plt.ylabel("Velocity (ft/s)")

    plt.legend(loc='upper left')
    fig.savefig(os.path.join(rootfolder,s+"_Vel.png"), dpi=150)
    plt.show()

    #END SETTING CHART UP

    def Atten(dist, Qp):
        Qpatt = Qp * ((0.0013*(dist**2)) - (0.0625*dist) + 1)
        #Qpatt = Qp* (-0.0424*dist + 0.9976)
        return round(Qpatt,1)

    segmi = seglength/5280

    Vrows = []
    def Vfunct(flow):
        for i in range(len(alist)):
            distmi = (segmi*(i+1)) - (segmi/2)
            distmi = round(distmi,2)
            flowatten = Atten(distmi, flow)
            V1 = alist[i] * flowatten**blist[i]
            V1 = round(V1,2)
            Vrows.append((i+1, distmi, flowatten, V1))
        return(Vrows)

    Vfunct(Qptot) #run function and input flow

    def TTPfunct(structdist):
        Varray = numpy.array(Vrows)
        Velocities = Varray[:,3] #get array of velocities
        weightedgid = structdist / segmi
        gidmath = math.trunc(weightedgid) #gives whole number gid
        if gidmath == 0:
            timechunk = structdist*5280/Velocities[0] #convert dist to feet div by first vel
            timechunk = timechunk/3600 # convert seconds to hrs
            TTP = timechunk + (ftimehrs/2)
            return TTP
        elif gidmath > 0:
            Vsum = sum(Velocities[0:gidmath]) #sum all complete gid seg paths
            lenavg = len(Velocities[0:gidmath])
            Vavg1 = Vsum/lenavg #avg
            gidremain = weightedgid - gidmath #remaining section decimal
            Velremain = Velocities[gidmath] # velocity gid rounded up (trunc rounded down)

            bigdist = (gidmath * seglength) # completed dist in ft 
            fractdist = (gidremain * seglength) # fractional dist in ft
            timefract = fractdist/Velremain/3600 # convert to hrs
            timechunk = bigdist/Vavg1/3600
            TTP = timefract + timechunk + (ftimehrs/2)
            return TTP

    #add failure time to TTP calc

    def dam(invnum):
            with arcpy.da.SearchCursor(DamPoints,['IdNumber', 'MaxV', 'TopH'] ) as cursor:
                for row in cursor:
                    if row[0] == invnum:
                        dam.Id = row[0]
                        dam.V = row[1]
                        dam.H = row[2]
                        print(f"{dam.Id} || Dam Height: {dam.H} || Dam Volume: {dam.V}")

    def ftime08(Vol, Hbr):
        Vol = Vol * 43560 #acft to ft^3
        ftime_s = 63.2 * (Vol / (32.17*(Hbr**2)))**0.5
        ftime_min = round(ftime_s / 60,2)
        print("Btime: " + str(ftime_min) + " min (Frhlch, 2008)")
        return(ftime_min)


    dam(s)

    ftimehrs = ftime08(dam.V, dam.H)/60
    # Open TTP for structures and Roads

    arcpy.AddMessage("#---------Calculating Structure Time to Peak---------#")

    buildingfc = gdbname("Structures")

    try:
        arcpy.AddField_management(buildingfc, "TTP", "FLOAT", 6, 2)
    except:
        pass

    cursor = arcpy.da.UpdateCursor(buildingfc,['Dist2Dam', 'TTP'])
    for row in cursor:
        if row[0] == None:
            row[1] == None
        else:
            row[1] = TTPfunct(row[0])
        cursor.updateRow(row)
    del cursor


    arcpy.AddMessage("#---------Calculating Road Time to Peak---------#")

    roadfc = gdbname("RoadX")

    try:
        arcpy.AddField_management(roadfc, "TTP", "FLOAT", 6, 2)
    except:
        pass

    cursor = arcpy.da.UpdateCursor(roadfc,['Dist2Dam', 'TTP'])
    for row in cursor:
        if row[0] == None:
            row[1] == None
        else:
            row[1] = TTPfunct(row[0])
        cursor.updateRow(row)
    del cursor

    #ROAD LABEL CALLOUT STEPS - if callout not working might need to run TTP script
    roadcalloutlyr = gdbname("Roadcallout")
    arcpy.env.addOutputsToMap = False
    arcpy.analysis.Buffer(roadfc, roadcalloutlyr, "15 Feet", "FULL", "ROUND", "NONE", None, "PLANAR")

    # #layrx symbology update connection
    # lyr = m.listLayers(r'Road Callout')[0]
    # cp = lyr.connectionProperties
    # cp['connection_info']['database'] = gdb
    # cp['dataset'] = os.path.basename(roadcalloutlyr)
    # lyr.updateConnectionProperties(lyr.connectionProperties, cp)
    # arcpy.env.addOutputsToMap = True


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

    RoadcalloutFile = gdbname("Roadcallout")

    DamPoints = arcpy.GetParameterAsText(0)
    Roadcallout = arcpy.SetParameter(1, RoadcalloutFile)

    ScriptTool(DamPoints, Roadcallout)