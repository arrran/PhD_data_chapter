# Here I download and process three datasets

1. Rema strips
2. icesat lines
3. cryosat lines


#How I loaded rema stripes 

On Qgis, i loaded the stripes (multipolygon) and the line.
Then used the function "select by location" to select all polygons which intersected the line. Then "save selected variables as"

from attribute table get name.


#DIdnt do the following:

#to make features a variable:

variable = iface.activeLayer()

#then

line_feature = channel_line
area_layer = rema_strips

output_strips = []
cands = area_layer.getFeatures(QgsFeatureRequest().setFilterRect(line_feature.geometry().boundingBox()))
for area_feature in cands:
    if line_feature.geometry().intersects(area_feature.geometry()):
        areas.append(area_feature.id())

area_layer.select(areas)
