import bpy 
#---------------------------------Definitions---------------------------------------# 
path_tox3d='Z:/Blender/Sphere/blender files/' 
path_toimage='Z:/Blender/Sphere/blender files/' 
filename_ofx3d='Sphere_' 
filename_ofimage='Sphereimg_' 
startframe_ofx3d=0 
endframe_ofx3d=494
#-----------------------------------------------------------------------------------# 

bpy.ops.object.add() 
bpy.data.objects['Empty'].select=True 
bpy.data.objects['Empty'].name="Aux_Empty" 
bpy.context.active_object.location =(0.0, 0.0, 0.0) 
bpy.data.objects['Aux_Empty'].select=False 
  
from math import radians

camera = bpy.context.object
camera.location = (1.0, 0.0, 1.0)
camera.rotation_euler = (radians(45), 0.0, radians(30))

for currentframe in range(startframe_ofx3d,(endframe_ofx3d+1)):    
    bpy.context.scene.frame_set(currentframe) 
    bpy.data.objects['Aux_Empty'].select=True 
    bpy.context.scene.objects.active = bpy.data.objects['Aux_Empty'] 
    obj = bpy.context.object 
    obj.location[2] = 0.0 
    obj.keyframe_insert(data_path='location') 
    bpy.data.objects['Aux_Empty'].select=False 
    #---------------------------Import the x3d file--------------------------------# 
    bpy.ops.import_scene.x3d(filepath=path_tox3d+filename_ofx3d+str(currentframe)+'.x3d', filter_glob="*.x3d;*.wrl", axis_forward='Z', axis_up='Y') 
    #-------------------------------------------------------------------------------# 
	
    #------------rename the geometry to "imported geometry"----------------#
    #--Change the "ShapeIndexedFaceSet" according to your geometry--#
    obs = []
    for ob in bpy.data.objects:
        if ob.name.startswith("Shape_IndexedFaceSet"):
            ob.select = True
            bpy.context.scene.objects.active = ob
            bpy.context.active_object.material_slots[0].material.use_vertex_color_paint=True
            ob.select = False


	#-----------------------Set the rendering options-----------------------------#
    bpy.data.scenes['Scene'].render.resolution_percentage=100 

    bpy.ops.render.render(write_still=True) 

    bpy.data.images['Render Result'].file_format='PNG' 
    bpy.data.images['Render Result'].save_render(filepath=path_toimage+filename_ofimage+str(currentframe)+'.png')
    
    for ob in bpy.data.objects:
        if ob.name.startswith("Shape_IndexedFaceSet"):
            ob.select = True
    bpy.ops.object.delete()
    
    for ob in bpy.data.objects:
        if ob.name.startswith("DirectLight"):
            ob.select = True
    bpy.ops.object.delete()

	#--------------------------------End of the loop------------------------------#
bpy.data.objects['Aux_Empty'].select=True 
bpy.ops.object.delete()

