# Renders the Phi=0 (droplet interface) and Phi2=0 (solid interface) isosurfaces
# of the exported terminal timestep as 3D contours, mirrors them across the
# x=0 and z=0 symmetry planes to reconstruct the full droplet (the simulated
# domain only covers x, z in [-L, 0]), and saves a screenshot.
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


def mirrored_contour(fieldName, registrationPrefix, color, opacity=1.0):
    """Contour(fieldName == 0), then mirror across x=0 and z=0 (the domain's own
    x/z maximum, since the simulated quarter-domain is x, z in [-L, 0]) so the
    full droplet is reconstructed from the simulated quarter."""
    contour = Contour(registrationName=f'{registrationPrefix}Contour', Input=reader)
    contour.ContourBy = ['POINTS', fieldName]
    contour.Isosurfaces = [0.0]
    contour.PointMergeMethod = 'Uniform Binning'

    reflectX = Reflect(registrationName=f'{registrationPrefix}ReflectX', Input=contour)
    reflectX.Plane = 'X Max'
    reflectX.CopyInput = 1

    reflectZ = Reflect(registrationName=f'{registrationPrefix}ReflectZ', Input=reflectX)
    reflectZ.Plane = 'Z Max'
    reflectZ.CopyInput = 1

    display = Show(reflectZ, renderView1, 'GeometryRepresentation')
    display.Representation = 'Surface'
    display.DiffuseColor = color
    display.Opacity = opacity
    return reflectZ


mirrored_contour('Phi', 'Phi', [0.65, 0.85, 1.0], opacity=0.5)     # droplet interface, light blue, translucent
mirrored_contour('Phi2', 'Phi2', [0.65, 0.9, 0.65])                # solid interface, light green

renderView1.ResetCamera(True)
renderView1.Update()

SaveScreenshot(outputPng, renderView1, ImageResolution=[1024, 1024])
print(f"Wrote {outputPng}")
