# Getting all the dependencies

I suggest using Conda for environment management, as it helps ensure consistent functionality.

To install Conda, please refer to the following guide: [Miniconda Installation](https://docs.anaconda.com/miniconda/).

To create a new environment with all the dependencies, use the following command in the Github directory of the project:

```
conda env create -f environment.yml
```

This will generate a new environment named "geo" for this project.

# Generating Points on a Map

The scripts use a MultiLineString in a `.csv` file to generate faults. This allows the user to specify the exact latitude and longitude for each point on the fault.

The `.csv` file should follow this format:

| id | Name | to\_join | multilinestring |
| :---- | :---- | :---- | :---- |
| 923 | Calaveras (north) | a1 | MULTILINESTRING((-121.81311 37.45749000000001 0.0,\<x\> \<y\> \<z\>,\<x\> \<y\> \<z\>)) |
| 902 | Hayward (north) | a2 | MULTILINESTRING((-121.81311 37.45749000000001 0.0,\<x\> \<y\> \<z\>,\<x\> \<y\> \<z\>)) |

The `.csv` file can have more than four columns, but the final column needs to hold the MultiLineString. The second column should be the **Name** column, and the third column should be the **to\_join** column.

To visualize the points with a KML file, use the merge script to generate a KML file for visualization.

```
python merge_and_smooth.py <location_of_file>.csv --disable-smoothing --kml-output <output_location>.kml
```

[This website is very useful for visualizing the MultiLineStrings for each of the faults](https://editor.gislayer.com/)  
[Could also use google earth for visualization and generation](https://earth.google.com/)

## Merging Faults

The merge\_and\_smooth.py script can also merge two or more lines together. The **to\_join** column is used to specify which faults should be joined. The algorithm merges lines based on the values in this column, which should follow the pattern `<identifier><number>`. For example, use values like `a1`, `b1`, `a2`. The script merges all items with the same identifiers together in numerical order.

The to\_join column can also have empty values

## Smoothing the faults

the script smooths the faults using a [Centripetal Catmull–Rom](https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline) spline to interpolate between the points and create a smooth line. Smoothing can be disabled with the `--disable-smoothing` tag.

Can use the –resolution tag to change the distance between the points. As well as –alpha tag to change how the cat mull rom behaves (more information on the wiki above).  
![][image1]

The \`--plot\` tag can be used to visualize the lines with matplotlib. The blue line is unsmoothed and the red one is smoothed 

### Sample Command for generating smoothed line

```
python merge_and_smooth.py outputs\BayModel1\output.csv --csv-output outputs\BayModel1\merged_catmul_500m_a0.5.csv --resolution 500 --plot
```

## Extracting points from a .kmz file

I have also written a Jupyter notebook to help extract MultiLineStrings from different datasets, such as the fault and fold datasets. The notebook is named `KMLEdit.ipynb` and is located in the meshgeneration repository. This notebook is a bit chaotic, you would need to understand the code to use it, sorry :( 

* Change the KMZ file to load in the second cell of the notebook.  
* Use the third cell to filter based on the name in the dataset.  
* Use the final cell to write the data back to a CSV file if required.

# Generating Topography and Fault STLs

The bulk of the work is in the `meshgen` script. It takes in a CSV file with relaxed requirements: it only needs a **Name** column, and the final column should contain the MultiLineStrings.

The `meshgen` script retrieves topography from the USGS’s 3D elevation dataset and generates the topography as well as the faults in relation to this topography, so you won't need to align the two manually. It also generates a bounding box for both the faults and the topography, which can be used as the topography if you want to create a simpler mesh.

## Example Command

```
python meshgen.py outputs\BayModel1\merged_catmul_500m_a0.5.csv 
--fault-output outputs\BayModel1\fault.stl 
--fault-resolution 200 
--topography-output outputs\BayModel1\topo.stl 
--topography-step 10 
--surrounding-region 0.7 
--bounding-box-output outputs\BayModel1\bbox.stl 
--bb-depth-below-fault 14 
--bb-mesh-size 1000 
--meta-data-output outputs\BayModel1\meta 
--bb-height-above-topography 4 
--bb-distance-from-topography 4 
--fault-depth 14 
--fault-height 5
```

This command takes a bit of time to finish, about 10-15 minutes.  
The output should look something like this

```
Using Fast Path
Downloading topography for region: (np.float64(-123.69627), np.float64(35.72582999999996), np.float64(-120.27386), np.float64(39.462520000000026))
Num points for topography : (1657040, 3)
Generating faults for 2 sections
Generating Walls: 100%|██████████████████████████████████████████████████████████████████| 2/2 [00:00<00:00,  8.57it/s]
Generating bounding box
Info	: Meshing 1D...
Info	: [  0%] Meshing curve 1 (Line)
Info	: [ 10%] Meshing curve 2 (Line)
Info	: [ 20%] Meshing curve 3 (Line)
Info	: [ 30%] Meshing curve 4 (Line)
Info	: [ 40%] Meshing curve 6 (Line)
Info	: [ 50%] Meshing curve 7 (Line)
Info	: [ 50%] Meshing curve 8 (Line)
Info	: [ 60%] Meshing curve 9 (Line)
Info	: [ 70%] Meshing curve 11 (Line)
Info	: [ 80%] Meshing curve 12 (Line)
Info	: [ 90%] Meshing curve 16 (Line)
Info	: [100%] Meshing curve 20 (Line)
Info	: Done meshing 1D (Wall 0.0062671s, CPU 0s)
Info	: Meshing 2D...
Info	: [  0%] Meshing surface 1 (Plane, Frontal-Delaunay)
Info	: [ 20%] Meshing surface 13 (Surface, Frontal-Delaunay)
Info	: [ 40%] Meshing surface 17 (Surface, Frontal-Delaunay)
Info	: [ 50%] Meshing surface 21 (Surface, Frontal-Delaunay)
Info	: [ 70%] Meshing surface 25 (Surface, Frontal-Delaunay)
Info	: [ 90%] Meshing surface 26 (Plane, Frontal-Delaunay)
Info	: Done meshing 2D (Wall 15.6381s, CPU 13.0781s)
Info	: 333871 nodes 670712 elements
Saving bounding box mesh outputs\BayModel1\bbox.stl
Info	: Writing 'outputs\BayModel1\bbox.stl'...
Info	: Done writing 'outputs\BayModel1\bbox.stl'
Saving topography : outputs\BayModel1\topo.stl
Saving faults : outputs\BayModel1\fault.stl
Center: [-8.17344854e-08 -6.13386789e-08  6.36873824e+06]
Rotational Matrix: [[ 0.88988707 -0.1763185   0.42072889]
 [-0.1763185   0.71766974  0.67369276]
 [-0.42072889 -0.67369276  0.6075568 ]]
Generating meta file at outputs\BayModel1\meta.h5
```

## Changing the topography resolution

The script can be used to determine the available resolutions for a CSV file. Use the `--just-check-res` tag to generate a bounding box for all the different points in the CSV. Then, expand that bounding box with the specified expansion term `--surrounding-region <number>`, which defaults to 0.01 (in latitude and longitude). The script will then output the different resolutions that can be used for that region.

The output will look something like this:

```
{'1m': True, '3m': True, '5m': False, '10m': True, '30m': True, '60m': False, 'topobathy': True}
```

## Fault Options

### **`--fault-output TEXT`**

**Description:** Specifies the filepath where the fault output will be saved. This output file contains the stl. If not provided, the output will not be saved. 

Make sure that the output also specifies the file format as it is used to infer the format.   
These are the supported formats  
[Abaqus](http://abaqus.software.polimi.it/v6.14/index.html) (`.inp`), ANSYS msh (`.msh`), [AVS-UCD](https://lanl.github.io/LaGriT/pages/docs/read_avs.html) (`.avs`), [CGNS](https://cgns.github.io/) (`.cgns`), [DOLFIN XML](https://manpages.ubuntu.com/manpages/jammy/en/man1/dolfin-convert.1.html) (`.xml`), [Exodus](https://nschloe.github.io/meshio/exodus.pdf) (`.e`, `.exo`), [FLAC3D](https://www.itascacg.com/software/flac3d) (`.f3grid`), [H5M](https://www.mcs.anl.gov/~fathom/moab-docs/h5mmain.html) (`.h5m`), [Kratos/MDPA](https://github.com/KratosMultiphysics/Kratos/wiki/Input-data) (`.mdpa`), [Medit](https://people.sc.fsu.edu/~jburkardt/data/medit/medit.html) (`.mesh`, `.meshb`), [MED/Salome](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/med-file.html) (`.med`), [Nastran](https://help.autodesk.com/view/NSTRN/2019/ENU/?guid=GUID-42B54ACB-FBE3-47CA-B8FE-475E7AD91A00) (bulk data, `.bdf`, `.fem`, `.nas`), [Netgen](https://github.com/ngsolve/netgen) (`.vol`, `.vol.gz`), [Neuroglancer precomputed format](https://github.com/google/neuroglancer/tree/master/src/neuroglancer/datasource/precomputed#mesh-representation-of-segmented-object-surfaces), [Gmsh](https://gmsh.info/doc/texinfo/gmsh.html#File-formats) (format versions 2.2, 4.0, and 4.1, `.msh`), [OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) (`.obj`), [OFF](https://segeval.cs.princeton.edu/public/off_format.html) (`.off`), [PERMAS](https://www.intes.de) (`.post`, `.post.gz`, `.dato`, `.dato.gz`), [PLY](https://en.wikipedia.org/wiki/PLY_\(file_format\)) (`.ply`), [STL](https://en.wikipedia.org/wiki/STL_\(file_format\)) (`.stl`), [Tecplot .dat](http://paulbourke.net/dataformats/tp/), [TetGen .node/.ele](https://wias-berlin.de/software/tetgen/fformats.html), [SVG](https://www.w3.org/TR/SVG/) (2D output only) (`.svg`), [SU2](https://su2code.github.io/docs_v7/Mesh-File/) (`.su2`), [UGRID](https://www.simcenter.msstate.edu/software/documentation/ug_io/3d_grid_file_type_ugrid.html) (`.ugrid`), [VTK](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf) (`.vtk`), [VTU](https://vtk.org/Wiki/VTK_XML_Formats) (`.vtu`), [WKT](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) ([TIN](https://en.wikipedia.org/wiki/Triangulated_irregular_network)) (`.wkt`), [XDMF](https://xdmf.org/index.php/XDMF_Model_and_Format) (`.xdmf`, `.xmf`).

**Default:** None

### **`--fault-height INTEGER`**

**Description:** Defines the vertical height of the fault relative to the topography. This value, in kilometers, specifies how high the fault extends above the surface.

**Default:** 2 km

### **`--fault-depth INTEGER`**

**Description:** Sets the depth of the fault below the topography, measured in kilometers. 

**Default:** 4 km

### **`--fault-resolution FLOAT`**

**Description:** Determines the resolution of the fault model by specifying the size of the triangles used to represent the fault surface, measured in meters. A smaller value results in higher resolution with finer detail, whereas a larger value results in a coarser representation. This setting influences the accuracy of the fault representation and computational load, so choose a resolution that balances detail with performance needs.

**Default:** 50 meters

## Topography Options

#### **`--topography-output TEXT`**

**Description:** Defines the filepath where the topography output will be saved. This file will be the stl file for the topography mesh. If not provided, the output will not be saved.

Make sure that the output also specifies the file format as it is used to infer the format.   
These are the supported formats  
[Abaqus](http://abaqus.software.polimi.it/v6.14/index.html) (`.inp`), ANSYS msh (`.msh`), [AVS-UCD](https://lanl.github.io/LaGriT/pages/docs/read_avs.html) (`.avs`), [CGNS](https://cgns.github.io/) (`.cgns`), [DOLFIN XML](https://manpages.ubuntu.com/manpages/jammy/en/man1/dolfin-convert.1.html) (`.xml`), [Exodus](https://nschloe.github.io/meshio/exodus.pdf) (`.e`, `.exo`), [FLAC3D](https://www.itascacg.com/software/flac3d) (`.f3grid`), [H5M](https://www.mcs.anl.gov/~fathom/moab-docs/h5mmain.html) (`.h5m`), [Kratos/MDPA](https://github.com/KratosMultiphysics/Kratos/wiki/Input-data) (`.mdpa`), [Medit](https://people.sc.fsu.edu/~jburkardt/data/medit/medit.html) (`.mesh`, `.meshb`), [MED/Salome](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/med-file.html) (`.med`), [Nastran](https://help.autodesk.com/view/NSTRN/2019/ENU/?guid=GUID-42B54ACB-FBE3-47CA-B8FE-475E7AD91A00) (bulk data, `.bdf`, `.fem`, `.nas`), [Netgen](https://github.com/ngsolve/netgen) (`.vol`, `.vol.gz`), [Neuroglancer precomputed format](https://github.com/google/neuroglancer/tree/master/src/neuroglancer/datasource/precomputed#mesh-representation-of-segmented-object-surfaces), [Gmsh](https://gmsh.info/doc/texinfo/gmsh.html#File-formats) (format versions 2.2, 4.0, and 4.1, `.msh`), [OBJ](https://en.wikipedia.org/wiki/Wavefront_.obj_file) (`.obj`), [OFF](https://segeval.cs.princeton.edu/public/off_format.html) (`.off`), [PERMAS](https://www.intes.de) (`.post`, `.post.gz`, `.dato`, `.dato.gz`), [PLY](https://en.wikipedia.org/wiki/PLY_\(file_format\)) (`.ply`), [STL](https://en.wikipedia.org/wiki/STL_\(file_format\)) (`.stl`), [Tecplot .dat](http://paulbourke.net/dataformats/tp/), [TetGen .node/.ele](https://wias-berlin.de/software/tetgen/fformats.html), [SVG](https://www.w3.org/TR/SVG/) (2D output only) (`.svg`), [SU2](https://su2code.github.io/docs_v7/Mesh-File/) (`.su2`), [UGRID](https://www.simcenter.msstate.edu/software/documentation/ug_io/3d_grid_file_type_ugrid.html) (`.ugrid`), [VTK](https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf) (`.vtk`), [VTU](https://vtk.org/Wiki/VTK_XML_Formats) (`.vtu`), [WKT](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) ([TIN](https://en.wikipedia.org/wiki/Triangulated_irregular_network)) (`.wkt`), [XDMF](https://xdmf.org/index.php/XDMF_Model_and_Format) (`.xdmf`, `.xmf`).

**Default:** None

#### **`--topography-resolution INTEGER`**

**Description:** Sets the resolution for the topography data in meters. This value needs to be one of the values that are available in the 3dep dataset. Use the `--just-check-res` to find the different options

**Default:** 30 meters

#### **`--topography-step INTEGER`**

**Description:** Defines the stride or sampling step for processing the topography data. This parameter controls how frequently data points are sampled from the topography, affecting both the resolution and the processing time. A step value of 1 indicates full resolution, while larger values will reduce the data density.

**Default:** 1

#### **`--surrounding-region FLOAT`**

**Description:** Specifies the size of the bounding box around the topographic region in latitude and longitude units. This parameter allows you to extend the area of interest away from the faults.

**Default:** 0.01

#### **`--topo-solver [vtk|scipy|custom]`**

**Description:** Selects the solver to be used for processing the point cloud data of the topography. Options include:

* `vtk`: Suitable for smaller point clouds but may crash with larger datasets.  
* `scipy`: Alternative solver for handling larger data.  
* `custom`: A specialized solver, though it does not support chunked downloads. (Also the fastest to compute)

**Default:** custom

#### **`--compare-solver / --no-compare-solver`**

**Description:** When enabled with `--compare-solver`, the tool will compare the generated topography with the results produced by VTK's Delaunay implementation. This comparison helps validate the accuracy of the topography generation process. If not enabled (`--no-compare-solver`), no comparison will be performed.

**Default:** \--no-compare-solver

#### **`--compare-topo-resample / --no-compare-topo-resample`**

**Description:** When enabled with `--compare-topo-resample`, a plot will be displayed comparing the topography before and after convolution. This visual comparison helps assess the impact of resampling on the topographic data.

**Default:** \--no-compare-topo-resample

## Extrusion Options 

In most cases not required, unless you want to duplicate the topography on the bottom of the mesh. This is disabled when the extrusion depth is 0 km

#### **`--extrusion-solver [custom|pyvista]`**

**Description:** Specifies the solver to be used for extruding the topography mesh. Options include:

* `custom`: A specialized solver for extrusion tasks.  
* `pyvista`: An alternative solver that is slower.

**Default:** custom

#### **`--extrude-surface-to-depth FLOAT`**

**Description:** Defines how deep the topography mesh should be extruded, measured in kilometers. This parameter allows you to extend the topographic surface into the subsurface, which can be useful for simulations requiring an extended depth profile.

**Default:** 0.0 km

## Bounding Box Options

The bounding box is generated with gmsh and will be in the same coordinate system as the topography mesh as well as the fault mesh.

#### **`--bounding-box-output TEXT`**

**Description:** Specifies the filepath for storing the bounding box mesh. If not provided, the bounding box will not be generated. 

Make sure that the output also specifies the file format as it is used to infer the format. 

This supports the following datatypes

* VTK (.vtk) \- Suitable for triangle meshes and visualizing in tools like ParaView.  
* STL (.stl) \- Can be used for triangle meshes in 3D printing and CAD applications.  
* GMSH Mesh Format (.msh) \- Supports triangle elements in 2D meshes.

**Default:** None

#### **`--bb-distance-from-topography INTEGER`**

**Description:** This setting defines the corners of the bounding box and how it will be positioned relative to the topography. Setting a value of 0 means the corners of the bounding box will align with the corners of the topography. To ensure the bounding box intersects slightly with the topography, adjust this value to create a thin strip that can be deleted in SimModeler.

If the mesh is not watertight in SimModeler, increase this value in increments of 1 and test if the newly generated mesh becomes watertight.

**Default:** 1 

#### **`--bb-mesh-size FLOAT`**

**Description:** Defines the size of the bounding box mesh, measured in meters. This parameter controls the granularity of the mesh, affecting both the detail and the computational load. A larger size results in a coarser mesh, while a smaller size increases the resolution.

**Default:** 500 meters

#### **`--bb-depth-below-fault FLOAT`**

**Description:** Specifies how deep the bounding box should extend below the topography, measured in kilometers. Make sure this is more than the depth of the fault

**Default:** 2 km

#### **`--bb-height-above-topography FLOAT`**

**Description:** Defines the height of the bounding box above the topography, measured in kilometers. This parameter sets how high the bounding box extends above the topographic surface

**Default:** 0.5 km

#### **`--plot-bb / --no-plot-bb`**

**Description:** When enabled with `--plot-bb`, the bounding box will be displayed in the GMSH user interface before saving. This option allows you to visualize and verify the bounding box before finalizing the output. If not enabled (`--no-plot-bb`), the bounding box will be processed without visualization.

**Default:** \--no-plot-bb

## Miscellaneous Options

#### **`--plot / --no-plot`**

**Description:** When enabled with `--plot`, the tool will display the fault and topography mesh using pyvista. This visualization helps in reviewing the generated meshes and ensuring that they meet the expected criteria. If not enabled (`--no-plot`), the meshes will not be displayed.  
Enabling this makes the script take significantly longer to compute. 

**Default:** \--no-plot

#### **`--verbose / --no-verbose`**

**Description:** When enabled with `--verbose`, the tool will provide detailed output during processing, which can be useful for debugging and in-depth analysis. If not enabled (`--no-verbose`), the output will be less detailed.

**Default:** \--no-verbose

#### **`--fast-path-disabled / --no-fast-path-disabled`**

**Description:** When enabled with `--fast-path-disabled`, the fast path for processing will be disabled, and only the meshio library will be used. This option can be useful for troubleshooting or ensuring compatibility. If not enabled (`--no-fast-path-disabled`), the fast path will be used if available.

**Default:** \--no-fast-path-disabled

#### **`--meta-data-output TEXT`**

**Description:** Specifies the filename for the HDF5 file that will store the center and rotational matrix metadata. Providing this filepath ensures that the metadata is saved and available for further use, such as for data analysis or integration with other tools.

**The metadata also records the command used to generate the different meshes. This section is very important. Ensure that you generate the metadata if you plan to create `.nc` files or use the script to find coordinates for forced rupture, as those scripts rely on the metadata to generate the necessary information.**

**Default:** None

# Generating Tetrahedral Mesh (with SimModeler)

## Preparing the mesh

To generate the mesh with SimModeler:

1. Import the bounding box, topography, and faults as discrete meshes.  
2. In the Discrete tab, use the **Union** tool to combine the bounding box with the topography.  
3. Use the **Delete** tool to remove the thin section outside the bounding box. You should now have two regions in SimModeler.  
4. If you don't see two regions, go back to generating the meshes and use the `-bb-distance-from-topography` tag to increase the intersection distance between the bounding box and the topography.  
5. Delete the region above the topography.  
6. Use the **Union** tool to combine the remaining region with the fault walls. This should separate the faults into those above and below the topography.  
7. Delete all faults floating above the topography.

## Generating the Tetrahedrons

1. **Mesh Case Selection**  
   * Ensure **Mesh case 1** is already present.  
2. **Mesh Size Settings**  
   * **Volume Mesh Size**: Set to **Absolute**, value **5000**.  
     * *Example*: For Sumatra, 250 km not on fault.  
   * **Fault Mesh Size**: Set to **Absolute**, value **1000**.  
     * *Example*: For Sumatra, 2.5-5 km on fault.  
3. **Gradation Rate**  
   * Unselect all faces.  
   * Set to **0.15**.  
     * *Example*: For Alice, Gambit, use \<0.1 to avoid reflections.  
     * *Example*: For Sumatra, 0.15 to 0.2 (Thomas used 0.3).  
4. **Surface/Volume Shape Metric**  
   * Select **Aspect Ratio** (value typically around **8**).  
5. **Surface/Volume Meshing**  
   * Set **Smoothing Level** to **4**.  
   * Choose **Type** as appropriate.  
6. **Mesh Size Propagation**  
   * Set **Distance** to **1000**.  
   * Set **Factor** to **2**.

After setting all these settings in the meshing tab. Press Generate mesh.

*This information is from Prof. Madden’s SeisSol notes* 

[I’ve uploaded the SimModeler documentation to the shared drive. If you want to understand how each of these variables affects the mesh, the documentation contains various examples for reference.](https://drive.google.com/open?id=1sUJ_0c8F7CD2adi5XMQY9e8REd-9CGI6)

## Specifying the boundary conditions

In the analysis tab

1. **Create a New Case**  
   * Right-click in the blank area of the cases box.  
   * Select **New** from the context menu.  
   * Enter the name **“Analysis case 1”** (this name will be used when meshing with PumGen).  
2. **Select Faces**  
   * Select the appropriate faces in your model.  
   * Click the **\+** tab and choose **BC** (Boundary Conditions).  
3. **Set Boundary Conditions**  
   * **Dynamic Rupture BC**: Apply this to the fault faces.  
   * **Absorbing BC**: Apply this to the side faces and bottom.  
   * **Free Surface BC**: Apply this to the surface/topography faces.

*This information is from Prof. Madden’s SeisSol notes*

After completing the steps, go to the File tab and select the Export Analysis Case button. Specify the filename and format, ensuring you use the .neu format, as this is the only format available for export. The export process uses a Lisp script located in the v001/ directory under case/seiSol in the SimModeler directory (gambit\_seisSol.sxp). This process may take a considerable amount of time (20-30 minutes).

If you need to export in a different file format later, the SDK for SimModeler should have been added to our license based on discussions with SimMertix representatives. You should be able to download it from their website. The SDK will allow you to write a custom C++ export plugin to handle tetrahedrons and boundary conditions in other file formats.

For very large meshes, the export might fail if there are more than 2^9 tetrahedrons in the mesh. In such cases, open the Lisp script and replace all occurrences of `i9` with `i15`. Alternatively, you can use the script stored in the Mesh-Generation repository.

### About gambit files

About \- [https://web.stanford.edu/class/me469b/handouts/gambit\_write.pdf](https://web.stanford.edu/class/me469b/handouts/gambit_write.pdf)

## Converting mesh for use with Seissol

You can use PumGen to convert the mesh into the H5/XDMF format required by SeisSol. Documentation for compiling PumGen can be found in the [SeisSol compilation document](https://docs.google.com/document/u/2/d/1cQN3C2WDc-FRLzx5l7yBwXLy4xcXKQJh-JS2688VcRY/edit).

```

./pumgen <location_of_simmodeler_export>.neu <filename_for_h5_file>
```

# Finding Coordinates for Forced Rupture

You can use the `get_point.py` script to find the XYZ coordinates of a latitude and longitude within the mesh to set the rupture position. Note that the script requires the metafile to be generated for the mesh.

#### **`--lat TEXT`**

**Description:** Specifies the latitude coordinate to convert into a point on the mesh. This latitude value will be used to determine the corresponding point in the mesh model.

**Default:** 37.341206616740266

#### **`--long TEXT`**

**Description:** Defines the longitude coordinate to convert into a point on the mesh. Similar to latitude, this longitude value will determine the corresponding point in the mesh model. 

**Default:** \-121.88075569799896

#### **`--depth FLOAT`**

**Description:** Sets the depth of the point relative to sea level, measured in meters. A negative value indicates a depth below sea level, while a positive value would represent an elevation above sea level.

**Default:** \-1000.0 meters

#### **`--round-to-closest-point / --no-round-to-closest-point`**

**Description:** When enabled with `--round-to-closest-point`, the point will be adjusted to the nearest point on the fault line. This rounding ensures that the point is aligned with the fault geometry. If not enabled (`--no-round-to-closest-point`), the point will not be rounded and will retain its original coordinates.

**Default:** \--round-to-closest-point

# Generating NC files for Initialization

The get\_point.py script can also be used to generate various NetCDF (.nc) files for initializing the model in SeisSol. This requires the metafile for the mesh.

#### **`--point-field-resolution FLOAT`**

**Description:** Defines the resolution of the netCDF files, specified in meters. This parameter sets the granularity of the data points in the output files, with a smaller resolution providing more detailed data. With smaller values this scales exponentially so make sure the number is only as big as required. 

**Default:** 100 meters

#### **`--output-distance-from-faults TEXT`**

**Description:** Specifies the filepath for generating a netCDF file that contains distances from points to the fault lines. This output file will be created at the provided location.

**Default:** None

#### **`--output-distance-from-topo TEXT`**

**Description:** Defines the filepath for generating a netCDF file that contains distances from points to the topography. This file will be created at the specified location and will include distance measurements from each point to the nearest topographic feature, aiding in topographic analysis and modeling.

**Default:** None

#### **`--chunk-size INTEGER`**

**Description:** Sets the size of chunks used during the generation process. This is to make sure that the script does not crash due to insufficient RAM, because the nc file could easily take more RAM than available.

**Default:** 50

You can easily add additional functions for NetCDF files. The code for this is quite simple and is located at the end of the `get_point.py` script. If the function you want to generate cannot be created using the combinations of fault distance and topography distance NetCDFs in the EASI file, you can modify or extend the code accordingly.

[image1]: <data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAnAAAAE7CAYAAABdZLyRAAAbUklEQVR4Xu3de5TVZb3H8QcoNTrIyQ6nVWm14lgnzTx4O5WWCTNIICKmaTiKkk4IKGkXnEAFb2QqWgNjgoimghqB5p2LuMxzCkUpcLxgSi6Nm3JAFBxkz3yPv42znXkcYHjYl9+zP+/XWr/Fnu8Df8j+Dp+Pe8/e2xkAAACi4vwBAAAA0o0CBwAAEBkKHAAAQGQocAAAAJGhwAEAAESGAgcAABAZChwAAEBkKHAAAACRocABAABEhgIHAAAQGQocAABAZChwAAAAkaHAAQAARIYCBwAAEBkKHAAAQGQocAAAAJGhwAEAAESGAgcAABAZChwAAEBkKHAAAACRocABAABEhgIHAAAQGQocAABAZChwAAAAkaHAAQAARIYCBwAAEBkKHAAAQGQocAAAAJGhwAEAAESGAgcAABAZChwAAEBkKHAAAACRocABAABEhgIHAAAQGQocAABAZChwAAAAkaHABVq/fr29+eabXFxcXFxcXCW6kixWRYELlCwOAAAoHeUspsAFUl4aAADSQDmLKXCBlJcGAIA0UM5iClwg5aUBACANlLOYAhdIeWkAAEgD5SymwAVSXhoAANJAOYspcIGUlwYAgDRQzmIKXCDlpQEAIA2Us5gCF0h5aQAASAPlLKbABSrE0mQaM7Zg+QJ/DAAA2lCILI4FBS5QvpcmKW9urMte+0zYxz8GAACefGdxTChwgfK9NMkjb80FLrkAAMD25TuLY0JTCFSIpeEpVAAA2q8QWRwLClwg5aUBACANlLOYAhdIeWkAAEgD5SymwAVSXhoAANJAOYspcIHSsjQLFpjNdb3sjjv8EwAAyltasrgUKHCB0rY0mzebneuus+/st9o/AgCgLKUti4uJAhcozUszoHKjjXejbN06/wQAgPKR5iwuNApcoBiW5tFHze5x/W3cOP8EAID4xZDFhUKBCxTb0hz90fnWt8tj/hgAgGjFlsX5RIELFOvSfKvzIjth9z/6YwAAohNrFucDBS5Q7Evzi1+YPeD62MKF/gkAAHGIPYt3BQUuULkszerVZmPcJXZK//L47wEA6CiXLA5BgQtUjkvT44vrbbS71DIZ/wQAgPQpxyxuLwpcoHJemilTzB5yve3OO/0TAADSo5yzeEcocIFUluapp8zudCfa8uX+CQAApaWSxW2hwAVSW5rnnjOb7QbY5Mn+CQAApaGWxS1R4AKpLk3ykV3Jz8kNq97iHwEAUFSqWZygwAVSXppmAzrcY989aJU/BgCgKJSzmAIXSHlpfN/o9qL16fK4PwYAoKCUs5gCF0h5abblmJ4bbair88cAABSEchZT4AIpL82O/PxnTXb4Z172xwAA5JVyFlPgAikvTXv12m+FHdvpPn8MAEBeKGcxBS6Q8tLsrCGnvmuj3Hg+4QEAkFfKWUyBC6S8NKGuvdbsd67KVvHCVQBAHihnMQUukPLS7KoFC8zmuAp77DH/BACA9lPOYgpcIOWlyZeVK7d+TNf48f4JAAA7ppzFFLj31NXVWZcuXbK3q6qqrEePHrkz51z28ikvTb41NZkd6+6xqVP9EwAAtk05iz/cTERlMhmrr6+3fv362YwZM3Lzzp0726BBg3JfNzQ0ZBdGeWkKZcsWs0PcEzb6gHv8IwAAPkQ5iylwLVRXV9szzzzjj9ukvDTF0L272VxX4Y8BAMhRzmIKXCDlpSm2/fYzu9sdm32qFQCAZspZTIELpLw0pXLZZWYjP3mrPwYAiFLOYgpcIOWlKaXGRrPr3Y/8MQBAkHIWU+ACKS9NGnR3L9rGjf4UAKBEOYspcIGUlyYtvn5YY/YROQCAJuUspsAFUl6aNDn6U4vtrrv8KQBAgXIWU+ACKS9N2iQF7jsff8IfAwDKnHIWU+ACKS9NWk1zp9sVV/hTAEC5Us5iClwg5aVJsw0bzC5zv7C33/ZPAADlRjmLKXCBlJcmBr8a9YZ9seNyfwwAKCPKWUyBC6S8NDE5eY/ZdtoP3vXHAIAyoJzFFLhAyksTo1+6UVZZ6U8BADFTzmIKXCDlpYnZTZ3OtN//3p8CAGKknMUUuEDKSxO75K67051otbX+CQAgJspZTIELpLw05WTGLZvtYnexrVzpnwAA0k45iylwgZSXphzNn282x1XYiSc0+UcAgJRSzmIKXCDlpSl3Q05vtIddpf3yl/4JACBNlLOYAhdIeWmUHPrVTXa3O9YyGf8EAFBqyllMgQukvDSK7r576yc89D5krX8EACgR5SymwAVSXhp1X+r2f3ZB1zp/DAAoMuUspsAFUl4abNW9u9nk3Yb7YwBAkShnMQUukPLSoLUePcxu71jljwEABaacxRS4QMpLg7b17Gl2fYeh/hgAUCDKWUyBC6S8NNi+E04wG+9G+WMAQJ4pZzEFLpDy0qB9+ne63wb22+yPAQB5opzFFLhAykuDnXOdO9cG9G/0xwCAXaScxRS4QMpLgzA/HLzF7nX9/DEAIJByFlPgAikvDXbNtW6kPwIABFDOYgpcIOWlwa5LfjZuwQJ/CgDYGcpZTIELpLw0yI+BR6yxhQv9KQCgvZSzmAIXSHlpkD+nfG2JLVniTwEA7aGcxRS4QMpLg/wa+YW7bdkyfwoA2BHlLKbABVJeGuRfzaem2iuv+FMAwPYoZzEFLpDy0qAwxuw10R8BALZDOYspcIGUlwaFc5/r648AANugnMUUuEDKS4PCuqbDT+zBB/0pAMCnnMUUuEDKS4PCmzPH7Nsd/+SPAQAtKGcxBS6Q8tKgeK5x59kjj/hTAEBCOYspcIGUlwbFlbzZb5X7nT8GAHnKWUyBC6S8NCiNS7r9xnr18qcAoEs5iylwgZSXBqXT1GRW64bbq6/6JwCgRzmLKXDvqaursy5dumRvV1VVWY8ePXJnEydOtMMOOyz3dTPlpUHpzbyk3jq5Lf4YAKQoZzEF7n2ZTMbq6+utX79+NmPGjNx89uzZLX7XB5SXBukxr2OlbdjgTwFAg3IWU+BaGDx4sD377LOtZiNGjLDGxsZWs4Ty0iB95rujbN48fwoA5U05iylwgZSXBum0aJHZLe5Ufj4OgAzlLKbABVJeGqTb009v/Uiuq67yTwCgvChnMQUukPLSIA7JK1YvdzX231992z8CgLKgnMUUuEDKS4P4VB70hv3EXWVt/DgnAERLOYspcIGUlwbxmjjRbK7j3YABlAflLKbABVJeGsSv6pQmu9Bd4o8BICrKWUyBC6S8NCgf3/xGk12+21h/DABRUM5iClwg5aVB+fncPk02dq9f+2MASDXlLKbABVJeGpSvvTq/YzX73+OPU2PB8gWWacz4YwCilLOYAhdIeWlQ/v7drbaLjl3sj0tq7wl7mxvrshclDkBCOYspcIGUlwY6vu7+164ck44PW20ub8mVPBIHAMpZTIELpLw00PLuu2Yj3bWpeA+5pLztM2EffwxAlHIWU+ACKS8NNNXUmP3WVfvjklu1yuxSN9o2bfJPAJQ75SymwAVSXhpgvBuVys9aHTzwTTvbTbKNG/0TAOVIOYspcIGUlwZIvPWW2TzX01au9E9Kb+ZMs9vdD6zPgSv8IwBlRDmLKXCBlJcGaOnWW83Oc9f449SYPNlsjquwuXP9EwCxU85iClwg5aUB2tL/oNfs4M+m8OG49yU/K3ebG+SPAURMOYspcIGUlwbYnqGuzsam+NO5ZrsB/ghApJSzmAIXSHlpgB1J3npkpjvetmzxT9JhkjvbHwGIkHIWU+ACKS8N0F6XjcvYD9zt1tTkn5RejbvcHwGIjHIWU+ACKS8NsLN+fvpq298t9cclV/mtd2xoh9/6YwCRUM5iClwg5aUBQl3w7f+xvTqu88clN7DfZjv5a/X+GEDKKWcxBS6Q8tIAu+rC/5hun/iXd/1xSSUfFXaXOyGVT/cCaJtyFlPgAikvDZAvV3zyKvvsZ9LVmI7Yd5UtW+ZPAaSRchZT4AIpLw2Qb7V7/NQOOMCfls63911hL77oTwGkjXIWU+ACKS8NUCjXdjrf/vIXf1oald1fsr//3Z8CSBPlLKbABVJeGqCQnn/e7AHXx156yT8pvv5fWEKJA1JMOYspcIGUlwYolovGNFqtG25LS/gOJCtWmFV2mu+PAaSAchZT4AIpLw1QbK+8YvaQ620jRvgnxfOwq/RHAEpMOYspcIGUlwYopYEfe9B6f/xxf1wU81xP3mYESBHlLKbABVJeGiANzj7b7I/umKK/WvQ7e79oZ5yy2R8DKAHlLKbABVJeGiBNkhcZzHW97Omn/ZPC2bJl61O6mzb5JwCKSTmLKXCBlJcGSKOVK83OcFNt2NHFe9noueeaDe18sz8GUCTKWUyBC6S8NEDafWXfLTbRDbNMxj8pjOPcLHvgAX8KoNCUs5gCF0h5aYBYTJu29YUHyVuBFNrUqWbnfn+lPwZQQMpZTIELpLw0QGxmzizOK0jvv9/suF4b/DGAAlHOYgpcIOWlAWI1dKjZFHemP86rW24xu+MOfwqgEJSzmAIXSHlpgNjt/+UtNmOGP82fyi/+3Roa/CmAfFPOYgpcIOWlAcrBmjVm97m+2U95KISb3Wn+CECeKWcxBS6Q8tIA5eS++8ymu5OtsdE/2XWz3HH+CEAeKWcxBS6Q8tIA5egn52y2Q9yT/niX9Tr8HX8EIE+Us5gCF0h5aYByduFnbrTPfc6fhvvHP8y+secz/hhAHihnMQXufZlMxurr67O3nfvgr6Vr1642NHnp2vsaGhqyC6O8NICC5c9usm+6x23tWv8kzA3urIK/jQmgRjmLKXAtVFdXZ39dk/x08w4oLw2g5PdX/8MOdQvzUr5+2PEm27jRnwIIpZzFFLj3DBkyxLp165a9XVNTk5tPmzYtOz/rrLNys2bKSwMoOm3QFjvCPWbvP1AfrGKvp2z5cn8KIIRyFlPgAikvDaBs3Tqz2zpU2Sfd6/b88/5p+1QctNbmzvWnAHaWchZT4AIpLw2ArZKfjxvd8XI7yc2wZcv80+0784yMXXmlPwWwM5SzmAIXSHlpAHzY66+bHdLxKbvJnW4vveSftm3KFLMDdgt8GA+AdBZT4AIpLw2A7Vu92qxTh0Z7yPX2j9p0yUfG2vDh/hTAjihnMQUukPLSAGi/ujqz3VyDLexaYdv7uIfkVa7DPjrZDvvIIv8IwDYoZzEFLpDy0gAIU1ubvM9kkz3RtZd/1MqAw1+3Ce7H9uyz/gmAlpSzmAIXSHlpAOy6664z+7xbbuP3neof5axYYTbC/cb6fPll/wiAaWcxBS6Q8tIAyK/Pf3qz/eXGpf64lcN7bLRzdr/BHwPSlLOYAhdIeWkA5F/yliKL/rWnP/6QQ7+03ir2XOiPAUnKWUyBC6S8NAAKI3lfufnuKH/cpiM6P2Wnf2+DPwakKGcxBS6Q8tIAKKyD/u0Vf9Sm5JWrl7saW8gDchClnMUUuEDKSwOg8Prs/ogNOLbJH2/T+e5qW8Q7kECMchZT4AIpLw2A4rnS/czGjG5fkfvzn7c+BZs8MgcoUM5iClwg5aUBUFyZjNmP3PV2160N/lGbTjrJbGDnB/0xUHaUs5gCF0h5aQCURvJ5q2e6ybZ0++84knOeuyb7Z4BypZzFFLhAyksDoLSSn3U755vt+4G36mqzfke+5Y+BsqCcxRS4QMpLA6D0hgwxe+01f9q2224zu/pqfwrETzmLKXCBlJcGQDqc5m72R9t0dv9XbfFifwrETTmLKXCBlJcGQEo0NZnbiX/F5+/R195+258C8VLO4p341kdLyksDID0aG83muIpskevb1+ypp/zf0dol7kJ/BERLOYspcIGUlwZAejX99W/23MGDbPZHTrA+7oFssTvlFLOXX/7g9ySvZAXKgXIWU+ACKS8NgPhs2WI2dqxZVYfbrNI9bPNcT97wF9FTzmIKXCDlpQFQHmYNuNk6uMYd/lycG+ts7wl7+2Og5JSzmAIXSHlpAJSXv514qXV0mVZPszZLylvztWD5Av8YKCnlLKbABVJeGgDl6YW6efZf7mmbO/eDWfLIW3OByzRmPjgAUkA5iylwgZSXBkD5mzJyqR3nZtn69f4JkB7KWUyBC6S8NAB0zJpldr/7rs2b558ApaecxRS4QMpLA0DPhg1mQ12dTZninwClo5zFFLhAyksDQFfyxsFHdHjcBned7R8BRaecxRS4QMpLAwCJ886z7PvJ/elP/glQHMpZTIELpLw0ANDSunVmv3XVNn++fwIUlnIWU+ACKS8NALRl2TKz37kq459HFItyFlPgAikvDQBsT12d2RxX4Y+BvFPOYgpcIOWlAYD2OKZvo1V+jB+QQ+EoZzEFLpDy0gDAzvi1O8cfAXmhnMUUuEDKSwMAO+t8d7WtWuVPgV2jnMUUuEDKSwMAIX78Y7Mju7/qj4FgyllMgQukvDQAEKqpyewe198fA0GUs5gCF0h5aQBgV53rrrOGBn8K7BzlLKbABVJeGgDIh4O+ssmm3tjkj4F2U85iCpwl71lUZ126dMneds5lr2YXXnih9erVK/d1M+WlAYB8muYG+yOgXZSzmAL3vkwmY/X19dnbLQvc4sWLc7dbUl4aAMi3/9xjuT3zjD8Ftk85iylwLVRXV2d/XbNmTW7mF7iGhobswigvDQAUwsiRZoO/u9ofA9uknMUUuPcMGTLEunXrlr1dU1OTm0+bNo2nUAGgiFauNLu4w1h/DLRJOYspcIGUlwYACq2y60J/BHyIchZT4AIpLw0AFFryfnGT3Nn+GGhFOYspcIGUlwYAimWe6+mPgBzlLKbABVJeGgAoppvdaf4IyFLOYgpcIOWlAYBi+6GbYps3+1OoU85iClwg5aUBgFI46hOLbdkyfwplyllMgQukvDQAUCoDe663qVP9KVQpZzEFLpDy0gBAKdX+psm+dzyfoQrtLKbABVJeGgAotSefNOvRnX+H1SlnMQUukPLSAEAarF1rdugeS/wxhChnMQUukPLSAECa/NqdY5mMP4UC5SymwAVSXhoASJtBu8/0RxCgnMUUuEDKSwMAafRT9yu74gp/inKmnMUUuEDKSwMAafX222a3ulNsxQr/BOVIOYspcIGUlwYA0m7m5LX2afdPf4wyo5zFFLhAyksDALGYuvdFts8+/hTlQjmLKXCBlJcGAGIz0x1v48b5U8ROOYspcIGUlwYAYtTQYPazj06wr7h6e/FF/xQxUs5iClwg5aUBgNitWWN2kzvDevb0TxAT5SymwAVSXhoAKCfJU6sPuD5WdWKDf4SUU85iClwg5aUBgHI1YoTZH9xAO3/YO/4RUkg5iylwgZSXBgAUfH9Ag01yZ1tTk3+CtFDOYgpcIOWlAQAlo0aZzXEV9u67/glKTTmLKXCBlJcGABRNmmT2sKv0xygh5SymwAVSXhoAUJY8pVrjLreTvs9zq6WmnMUUuEDKSwMA2Gp4xzrredhb/hhFopzFFLhAyksDAGitX49/2nBXa6++6p+gkJSzmAIXSHlpAABte+EFs3tdP6voxdOrxaCcxRS4QMpLAwDYsdMqV9jF7mL761/9E+SLchZT4AIpLw0AoP1Wrza7z/W12lr/BLtKOYspcIGUlwYAsPOSV6+evvt0O7zTn62x0T9FCOUspsAFUl4aAMCuuflms7mul110kX+CnaGcxRS4QMpLAwDIn+M/8YgNcLOtocE/wY4oZzEFLpDy0gAA8u/ee80ecr1t+HD/BNuinMUUuEDKSwMAKKzKz79g39p/rT+GRzmLKXCBlJcGAFAcQwZn7B7X36ZP90+QUM5iClwg5aUBABRX8qrVQz72jF3kxtr69f6pLuUspsAFUl4aAEDpPPmk2Y/c9Va97yP+kRzlLKbABVJeGgBAOowaZTbdnWy/GrfJP5KgnMUUuEDKSwMASJ/kjYKHu1o78cAX/KOypZzFFLhAyksDAEi3M4c02lR3hmUy/kl5Uc5iClwLN954o40ZM6bV7MADD2z1dTPlpQEAxGHCBLMH3dH2xhv+SXlQzmIK3PsWLVpk9fX1VlNT4x/ZunXrcrcbGhqyC6O8NACAuDz8sNklbowteqK8PoRVOYspcO+54IILbMmSJTZnzhzbc889c4/C1dbW2qOPPmoVFRXen9BeGgBAnJYvNzvTTbZbp7zjH0VJOYspcIGUlwYAELdNm8z6uXtt7IjX/aOoKGcxBS6Q8tIAAMpD8srVGe4ke/xx/yQOyllMgQukvDQAgPLy2mtmtW643XCDf5JuyllMgQukvDQAgPK0ebPZKDfeRo6I4/1HlLOYAhdIeWkAAOVvkLvNjjlygz9OFeUspsAFUl4aAADSQDmLKXCBlJcGAIA0UM5iClwg5aUBACANlLOYAhdIeWkAAEgD5SymwAVSXhoAANJAOYspcIGUlwYAgDRQzmIKXCDlpQEAIA2Us5gCF0h5aQAAaK8/PPsHc2OddRjbwT/aZcpZTIELpLw0AAC0V1Lemq98U87i/P9tilBeGgAA2otH4AqDAhdIeWkAAEgD5SymwAVSXhoAANJAOYspcIGUlwYAgDRQzmIKXCDlpQEAIA2Us5gCF0h5aQAASAPlLKbABVJeGgAA0kA5iylwgZSXBgCANFDOYgpcIOWlAQAgDZSzmAIXSHlpAABIA+UspsAFcs5xcXFxcXFxlfhSpftfXiaU/+8jrZT/QUkr7pP04T5JH+6TuHBvAQAARIYCBwAAEBkKXOSOOuqo3O2qqqrc7c2bN9uoUaNyX6N4hg0blrs9bty43O3FixfnbqO4Wt4n++23X+423yel0/I+aYmn8UqnOU+OPPLIVvNVq1bZ9ddf32qG0uM7JeWScGn+Qc3nnnuu1dnSpUtzt5csWWLTp09vcWpWX1/f6mvkx/buk8Sll16a/fX222/3TsxOPfVUf4Q8aO99kpg7d26LE75PCqm99wmKY3vfJy3zxC9wiY0bN/ojlBgFLjLJozjNj+SsXLnS1qxZYxMnTrRly5bZ6NGjrba2NnuWFIXOnTu3/KMokJb3SfIoaHK/JJJf169fn7tPkn8wO3bsmPtzKJxt3SeTJk3K3U7wfVI827pPWs4TU6ZMsbq6utzXKJzmv/ukvDXnSeLggw/O/o/NrFmzsl/37t3bunXr1vKPIgUocAAAAJGhwAEAAESGAgcAABAZChwAAEBkKHAAAACRocABAABEhgIHAAAQGQocAABAZChwAAAAkaHAAQAARIYCBwAAEBkKHAAAQGQocAAAAJGhwAEAAESGAgcAABAZChwAAEBkKHAAAACRocABAABEhgIHAAAQGQocAABAZChwAAAAkaHAAQAARIYCBwAAEBkKHAAAQGQocAAAAJGhwAEAAESGAgcAABAZChwAAEBk/h/SOVm4/+ouaAAAAABJRU5ErkJggg==>
