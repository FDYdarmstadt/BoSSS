# Renders the Phi=0 isosurface (droplet interface) of the exported terminal
# timestep as a 3D contour and saves a screenshot.
#
# Usage: pvpython leidenfrost_contour.py <export_dir> <output_png>
#   <export_dir>  directory containing the BoSSS TecPlot export (state_*.plt)
#   <output_png>  path to write the screenshot to
#
# If no arguments are given, defaults to ./ParaViewLeidenfrost and
# ./ParaViewLeidenfrost/leidenfrost_contour.png (relative to the current
# working directory the script is invoked from).

import glob
import os
import sys

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

exportDir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(".", "ParaViewLeidenfrost")
outputPng = sys.argv[2] if len(sys.argv) > 2 else os.path.join(exportDir, "leidenfrost_contour.png")

pltFiles = sorted(glob.glob(os.path.join(exportDir, "state_*.plt")))
if not pltFiles:
    raise RuntimeError(f"No exported TecPlot file (state_*.plt) found in '{exportDir}'. "
                        f"Run the .Export().To('{exportDir}').Do() cell first.")
pltFile = pltFiles[-1]

reader = VisItTecplotBinaryReader(registrationName=os.path.basename(pltFile), FileName=pltFile)
reader.UpdatePipelineInformation()
# select every available mesh zone / point field -- avoids hardcoding zone/field names
reader.MeshStatus = list(reader.MeshStatus.Available)
reader.PointArrayStatus = list(reader.PointArrayStatus.Available)

renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.Background = [1.0, 1.0, 1.0]

contour1 = Contour(registrationName='PhiContour', Input=reader)
contour1.ContourBy = ['POINTS', 'Phi']
contour1.Isosurfaces = [0.0]
contour1.PointMergeMethod = 'Uniform Binning'

contourDisplay = Show(contour1, renderView1, 'GeometryRepresentation')
contourDisplay.Representation = 'Surface'
contourDisplay.DiffuseColor = [0.85, 0.1, 0.1]

renderView1.ResetCamera(True)
renderView1.Update()

SaveScreenshot(outputPng, renderView1, ImageResolution=[1024, 1024])
print(f"Wrote {outputPng}")
