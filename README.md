## Screening Level Inundation Study Toolbox
### An ArcGIS Pro toolset for hydraulic modeling, impact modeling, and cartographic production

**Basic Model Summary**
1) Downstream reach is performed based on NHD flow lines and an emperically estimated "danger reach" that is a function of the dam's size
2) A healthy buffer around this danger reach defines the model domain
3) Various calculations are performed to determin segment-average hydraulic property tables based on channel morphology and landcover friction coefficients to develop synthetic stage-discharge rating curves.
4) The dam breach discharge is calculated via an NWS empirical equation the dam parameters
5) The estimated discharge is "routed" through the reach segments with an attenuation factor applied
6) Impacts analysis on roads and structures are performed including arrival times based on velocity calculations.
7) Results are provided both spatially with symbology applied and as excel tables 

To read more of the technical details please see the [User guide & Technical Model Documentation](https://github.com/AEarth/SLIS-Toolbox/blob/main/User%20Guide%20%26%20Technical%20Model%20Documentation.pdf)

## Project Structure
* **Data**: Reference data (not included in repo); land cover, nhd tracenetwork, address/buildings, spatial template files, etc
* **Runs**: Where each map projects,  geodatabase, and outputs is stored individually by folder ids
* **Scripts**: The primary multi-step geoprocessing toolbox located here
* **Symbology Templates**: All symbology layer templates that support the cartographic automation

## Resulting Data and Map Products

### Example Figure Exports 
<table border="0">
 <tr>
    <td><b style="font-size:30px">Rating Curves</b></td>
    <td><b style="font-size:30px">Velocity Curves</b></td>
 </tr>
 <tr>
    <td> <img src="https://github.com/user-attachments/assets/8d287bb1-1f72-4338-bdda-a9edbf4625e5" alt="alt text" width="400"></td>
    <td> <img src="https://github.com/user-attachments/assets/e86fe069-bc20-49b6-b3d2-3017c2f9336e" alt="alt text" width="400"></td>
 </tr>
   <tr>
    <td><b style="font-size:30px">Empirical Discharge Comps</b></td>
    <td><b style="font-size:30px">Elevation Profile</b></td>
 </tr>
 <tr>
    <td> <img src="https://github.com/user-attachments/assets/a5eb443d-3ed4-4f0f-9111-c6d034392a82" alt="alt text" width="400"></td>
    <td> <img src="https://github.com/user-attachments/assets/7a9a83bc-2624-45a5-9a55-3effbee14849" alt="alt text" width="400"></td>
 </tr>
</table>

### Drainage Area Map
![Drainage Area Map](https://github.com/user-attachments/assets/180b689c-8e9c-4241-830e-857d903e3f28)

### Sample Map Panel
![Map Sample](https://github.com/user-attachments/assets/50815b9f-e67d-4a69-9c5a-4ea5796baa27)

### Impact Tabular Exports
[VA023013_RoadX.xls](https://github.com/user-attachments/files/18287267/VA023013_RoadX.xls)

[VA023013_Structures.xls](https://github.com/user-attachments/files/18287272/VA023013_Structures.xls)

### Complete Report Sample
Complete reports were assembled via the a mailmerge template that intakes the various figures and data fields to dynamically fill tables and figures. 
Manual review and modifications were required to review and format impact tables in the word document. 

[075004_DCRSIM_20220511.pdf](https://github.com/user-attachments/files/18287265/075004_DCRSIM_20220511.pdf)

