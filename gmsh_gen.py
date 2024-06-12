# import gmsh
# import sys
#
# # Initialize the Gmsh API
# gmsh.initialize()
#
# # Create a new model
# gmsh.model.add("points")
#
# # Sample list of points with xyz coordinates
# points = [
#     (0, 0, 0),
#     (1, 0, 0),
#     (0, 1, 0),
#     (0, 0, 1),
#     (1, 1, 1),
# ]
#
# # Add points to the model
# for i, (x, y, z) in enumerate(points, start=1):
#     gmsh.model.geo.addPoint(x, y, z, 1.0, i)
#
# # Synchronize the model to make sure the points are processed
# gmsh.model.geo.synchronize()
#
# # Optionally, you can add visualization options like a physical group
# gmsh.model.addPhysicalGroup(0, [i for i in range(1, len(points) + 1)], 1)
# gmsh.model.setPhysicalName(0, 1, "MyPoints")
#
# # Generate the mesh (although not strictly necessary for points)
# gmsh.model.mesh.generate(0)
#
# # Save the model to a file (optional)
# # gmsh.write("points.msh")
#
# # Launch the GUI to visualize the points
# gmsh.fltk.run()
#
# # Finalize the Gmsh API
# gmsh.finalize()
import gmsh  # Download gmsh.py, and libgmsh files from gmsh-sdk

model = gmsh.model
factory = model.occ
mesh = model.mesh

gmsh.initialize()
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.2)  # max mesh size
# gmsh.option.setNumber("General.Terminal", 1)  # can be useful
model.add("Example")  # add a model

 # Create overlapping boxes (can be any shapes)
b1 = (3, factory.addBox(0,0,0, 1, 1, 1))  # 1x1x1 box with "bottom left" at 0,0,0
b2 = (3, factory.addBox(0.5, 0.5, 0.5, 1, 1, 1))

 # Calculate the intersection. With removeObject, removeTool as False, this
 # creates an object if there is an intesection
intersect = factory.intersect([b1], [b2], removeObject=False, removeTool=False)[0]
factory.synchronize()

if len(intersect):
    # if there is an intersection, do what you want to do.
    factory.remove(intersect, True)  # remove created intersection objects
    factory.synchronize()

mesh.generate(3)

 # Visualise
model.setVisibility(model.getEntities(3),0)  # turn volumes off
gmsh.fltk.run()  # start gmsh
gmsh.finalize()