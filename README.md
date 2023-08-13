# TimberMeshGen
more documented readme coming
Please contact Patrick Grant at p21.grant@qut.edu.au with any questions and clarification.

The files are run in the following order:
- TimberMaskRefined
- TimberSpectral
- Splinefitting
- BoxSegementing (which calls RingZoning and ExportGeo)

  Then mesh the geometry in Gmsh and output the mesh file as a .m file. This is then read by ImportMesh3D (make sure to not clear the workspace)
