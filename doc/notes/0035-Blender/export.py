from paraview.simple import *
startFrame = 0
endFrame = 494
for i in range(startFrame,endFrame+1):
	AnimationScene = GetAnimationScene()
	AnimationScene.AnimationTime = i
	exporters = servermanager.createModule("exporters")
	x3d = exporters.X3DExporter(FileName="Z:/Blender/Sphere/blender files/Sphere_" + str(i) + ".x3d")
	view = GetActiveView()
	x3d.SetView(view.SMProxy)
	x3d.Write()
	del x3d
