# Renders the Phi=0 (droplet interface) and Phi2=0 (solid interface) isosurfaces
# of the exported terminal timestep as 3D contours, mirrors them across the
# z=0 symmetry plane (and, if requested, also across x=0) to reconstruct the
# full droplet, and saves a screenshot.
#
# Usage: pvpython leidenfrost_contour.py <export_dir> <output_png> [mirror] [zoom]
#   <export_dir>  directory containing the BoSSS TecPlot export (state_*.plt)
#   <output_png>  path to write the screenshot to
#   [mirror]      'xz' (default) mirrors across both x=0 and z=0 (the two cases
#                 whose domain only covers x, z in [-L, 0]); 'z' mirrors only
#                 across z=0 (the inclined-gravity/sliding case, which is not
#                 symmetric across x=0). For 'z', the camera is framed on the
#                 droplet's own bounding box, padded by 1/zoom on each side
#                 (it may have slid off to one side of the much larger solid
#                 plate, so fitting the whole scene leaves it tiny/off-center).
#   [zoom]        for 'xz': extra camera zoom factor on top of the whole-scene
#                 fit (default 1.4). For 'z': divides the droplet's own extent
#                 to get the padding margin around it -- bigger = tighter crop
#                 (e.g. zoom=2 pads by half the droplet size on each side).
#
# If no arguments are given, defaults to ./ParaViewLeidenfrost,
# ./ParaViewLeidenfrost/leidenfrost_contour.png, 'xz' and 1.4 (relative to the
# current working directory the script is invoked from).

import glob
import os
import sys

from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

exportDir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(".", "ParaViewLeidenfrost")
outputPng = sys.argv[2] if len(sys.argv) > 2 else os.path.join(exportDir, "leidenfrost_contour.png")
mirror = sys.argv[3] if len(sys.argv) > 3 else "xz"
zoom = float(sys.argv[4]) if len(sys.argv) > 4 else 1.4

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
    """Contour(fieldName == 0), then mirror across z=0 (and, if `mirror` == 'xz',
    also across x=0 -- the domain's own x/z maximum, since that quarter/half-domain
    is simulated over x, z in [-L, 0]) so the full droplet is reconstructed from the
    simulated portion."""
    contour = Contour(registrationName=f'{registrationPrefix}Contour', Input=reader)
    contour.ContourBy = ['POINTS', fieldName]
    contour.Isosurfaces = [0.0]
    contour.PointMergeMethod = 'Uniform Binning'

    last = contour
    if "x" in mirror:
        reflectX = Reflect(registrationName=f'{registrationPrefix}ReflectX', Input=last)
        reflectX.Plane = 'X Max'
        reflectX.CopyInput = 1
        last = reflectX

    reflectZ = Reflect(registrationName=f'{registrationPrefix}ReflectZ', Input=last)
    reflectZ.Plane = 'Z Max'
    reflectZ.CopyInput = 1

    display = Show(reflectZ, renderView1, 'GeometryRepresentation')
    display.Representation = 'Surface'
    ColorBy(display, None)  # disable scalar coloring, otherwise DiffuseColor is ignored
    display.DiffuseColor = color
    display.AmbientColor = color
    display.Opacity = opacity
    return reflectZ


dropletContour = mirrored_contour('Phi', 'Phi', [0.65, 0.85, 1.0], opacity=0.5)    # droplet interface, light blue, translucent
mirrored_contour('Phi2', 'Phi2', [0.65, 0.9, 0.65])                                # solid interface, light green

if mirror == "z":
    # the droplet may have slid off to one side of the (much larger) solid plate,
    # so fitting the camera on the whole scene leaves it tiny and off-center;
    # frame directly on the droplet's own (padded) bounding box instead. Padding
    # the bounds themselves, rather than fitting tight and then calling Zoom(),
    # gives a predictable result independent of the camera's projection/clipping
    # behavior (a large Zoom() factor on an already-tight fit can push the near
    # clip plane through the geometry and make it disappear).
    dropletContour.UpdatePipeline()
    xmin, xmax, ymin, ymax, zmin, zmax = dropletContour.GetDataInformation().GetBounds()
    marginX = (xmax - xmin) / zoom
    marginY = (ymax - ymin) / zoom
    marginZ = (zmax - zmin) / zoom
    renderView1.ResetCamera(xmin - marginX, xmax + marginX, ymin - marginY, ymax + marginY, zmin - marginZ, zmax + marginZ)
else:
    renderView1.ResetCamera(True)
    renderView1.Update()
    GetActiveCamera().Zoom(zoom)

renderView1.Update()

SaveScreenshot(outputPng, renderView1, ImageResolution=[1024, 1024])
print(f"Wrote {outputPng}")
