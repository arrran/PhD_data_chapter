#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 09:41:23 2020

@author: arran
"""
import numpy as np
import os
import sys

import matplotlib.pyplot as plt

import pandas as pd








# attribute_table_stripes_over_channel
df = pd.read_csv('/Users/home/whitefar/DATA/REMA_2m_strips/KAMB_CHANNEL/attribute_table_stripes_over_channel.txt',delimiter='\t')

#acquisition date
df.acquisition
df.acquisit_1
# both are the same

#stripe name:
df.name


df = pd.read_csv('/Users/home/whitefar/DATA/REMA_2m_strips/KAMB_CHANNEL/attribute_table_stripes_over_channel.txt',delimiter='\t')


df = pd.read_csv('/home/arran/PHD/DATA/REMOTE_SENSING/attribute_table_stripes_over_channel.txt',delimiter='\t')

for stripe in df.name:
    crop_modus











#!/usr/bin/env python
# coding: utf-8

# # Crop REMA stripes by the polygon from 3 (around channel)
# 
# based on https://rasterio.readthedocs.io/en/stable/topics/masking-by-shapefile.html 
# but also using tarfile
# 
# This script unzips the raster, crops it, then saves a tiff

# In[1]:


import rasterio as rio
import rasterio.mask
import tarfile
import fiona

import os
import sys
#working_dir = '/Volumes/arc_02/REMOTE_SENSING/'
#working_dir = "/Users/home/whitefar/DATA/Channel/line_channel_traced_for_50kms/"
working_dir = '/Users/home/whitefar/DATA/FIELD_ANT_19/'

sys.path.append(os.path.abspath(working_dir))
sys.path.append(os.path.abspath('/Users/home/whitefar/DATA/code/'))
#sys.path.append(os.path.abspath('/Users/home/whitefar/DATA/Channel/line_channel_traced_for_50kms/'))
os.chdir(working_dir)


# In[2]:


#zipped_stripes_path = "/Users/home/whitefar/DATA/Channel/line_channel_traced_for_50kms/test_stripe/"
zipped_stripes_path = "/Volumes/arc_02/REMOTE_SENSING/REMA/REMA_stripes_over_kamb_channel/"

#cropped_stripes_path = zipped_stripes_path
#cropped_stripes_path = "/Users/home/whitefar/DATA/Channel/line_channel_traced_for_50kms/cropped_stripes/"
cropped_stripes_path = "/Volumes/arc_02/REMOTE_SENSING/REMA/REMA_stripes_1920_fieldsite/"


zipped_stripes_fnames = [f for f in os.listdir(zipped_stripes_path) if os.path.isfile(os.path.join(zipped_stripes_path, f))]

temp_tiffs_path = "/Users/home/whitefar/DATA/Channel/line_channel_traced_for_50kms/stripe_tiffs/" #save them temporararily

#zipped_stripes_fnames


# In[3]:


#these filenames came from looking at the polygon file of stripes in qgis

for_mosaic = ["SETSM_WV01_20161109_1020010058134D00_10200100576C9100_seg1_2m_v1.0",
              "SETSM_WV01_20161027_1020010056BA2E00_1020010058DFD800_seg1_2m_v1.0",
             "SETSM_WV02_20161220_1030010061866000_103001006007BA00_seg1_2m_v1.0",
             "SETSM_WV03_20161220_1040010026391800_1040010026480500_seg1_2m_v1.0"]

zipped_stripes_fnames = [name + '.tar.gz' for name in for_mosaic]
zipped_stripes_fnames


# In[4]:


with fiona.open("/Users/home/whitefar/DATA/FIELD_ANT_19/study_area_buffer.shp", "r") as shapefile:
    field_area = [feature["geometry"] for feature in shapefile]


# In[ ]:


with fiona.open("/Users/home/whitefar/DATA/REMA_2m_strips/polygons_of_tiffs/objectid_145041.gpkg", "r") as shapefile:
    field_area = [feature["geometry"] for feature in shapefile]


# In[5]:


def crop_stripe(crop_to_area,zipped_stripe_fname, zipped_stripes_path, cropped_stripes_path, temp_tiffs_path):
       
    #remove .tar.gz from the name
    stripe_name = zipped_stripe_fname[:-7]
    
    zipped_stripe_path =  zipped_stripes_path + zipped_stripe_fname 
    
    #path for output tiff
    tiff_stripe_fname = stripe_name + "_dem.tif"
    
    #extract the .tar.gz and save as .tiff
    with tarfile.open(zipped_stripe_path) as tar:
        stripe_tiff = tar.extract(member=tiff_stripe_fname, path=temp_tiffs_path)
        
    #crop the .tiff    
    with rio.open(temp_tiffs_path + tiff_stripe_fname) as src:
        #print(src.transform)
        out_image, out_transform = rasterio.mask.mask(src, crop_to_area,crop=True)
        data_type = out_image.dtype
        out_meta = src.meta.copy()
        
    #print(out_transform)
    
    out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform,
                 "dtype" : data_type})
    
    cropped_stripe_path = cropped_stripes_path + tiff_stripe_fname
    
    #write the cropped tiff to file
    with rio.open(cropped_stripe_path, "w", **out_meta) as dest:
        dest.write(out_image)
    
    #remove the uncropped tiff
    os.remove(temp_tiffs_path + tiff_stripe_fname)
        
    print(stripe_name + "cropped and written to "+ cropped_stripe_path)




# In[6]:


for i, zipped_stripe_fname in enumerate(zipped_stripes_fnames):
    
    crop_stripe(field_area, zipped_stripe_fname, zipped_stripes_path, cropped_stripes_path, temp_tiffs_path)
    
    print(f"{i+1}/{len(zipped_stripes_fnames)}")
       


# In[ ]:

#!/usr/bin/env python
# coding: utf-8

# # Merge tiffs
# 
# Merge the tiffs of cropped over the channel area / field site radar lines into one tiff.

# In[1]:


import numpy as np
import rasterio 
import rasterio.merge
import os
import sys
import fiona

#working_dir = '/Volumes/arc_02/REMOTE_SENSING/'
working_dir = "/Users/home/whitefar/DATA/FIELD_ANT_19/2m_REMA_mosaic_DEM/polygons_of_tiffs"


sys.path.append(os.path.abspath(working_dir))
sys.path.append(os.path.abspath('/Users/home/whitefar/DATA/code/'))
os.chdir(working_dir)


# In[2]:


polygons_path = "/Users/home/whitefar/DATA/Channel/line_channel_traced_for_50kms/cropped_stripes/"
cropped_REMA_stripes_path = "/Volumes/arc_02/REMOTE_SENSING/REMA/REMA_stripes_1920_fieldsite/"

output_path = "/Users/home/whitefar/DATA/FIELD_ANT_19/2m_REMA_mosaic_DEM/"

selected_stripes_fnames = [f for f in os.listdir(cropped_REMA_stripes_path) if os.path.isfile(os.path.join(cropped_REMA_stripes_path, f))]
selected_stripes_fnames


# In[3]:


# selected_stripes_polygon = "/Users/home/whitefar/DATA/FIELD_ANT_19/2m_REMA_mosaic_DEM/polygons_of_tiffs/polygons_for_merge.gpkg"
# with fiona.open(selected_stripes_polygon, "r") as shapefile:
#     selected_stripe_names = [feature["properties"]["name"] for feature in shapefile]
# selected_stripes_fnames = [names + "_dem.tif"  for names in selected_stripe_names]
# selected_stripes_fnames


# In[6]:


from contextlib import ExitStack

with ExitStack() as stack:
    files = [stack.enter_context(rasterio.open(cropped_REMA_stripes_path + fname)) for fname in selected_stripes_fnames]
    # All opened files will automatically be closed at the end of
    # the with statement, even if attempts to open files later
    # in the list raise an exception




# with rasterio.open(cropped_REMA_stripes_path + selected_stripes_fnames[0],"r") as src0,\
#      rasterio.open(cropped_REMA_stripes_path + selected_stripes_fnames[1],"r") as src1,\
#      rasterio.open(cropped_REMA_stripes_path + selected_stripes_fnames[2],"r") as src2,\
#      rasterio.open(cropped_REMA_stripes_path + selected_stripes_fnames[3],"r") as src3:
    
    print("1/4")
    
    out_image, out_transform =  rasterio.merge.merge(datasets=files,
                                                     bounds=None,
                                                     res=None,
                                                     nodata=None,
                                                     precision=7,
                                                     indexes=None)
    print("2/4")
    
    out_meta = files[0].meta.copy()
    data_type = out_image.dtype
        
    print("3/4")
    
    out_meta.update({"driver": "GTiff",
                     "height": out_image.shape[1],
                     "width": out_image.shape[2],
                     "transform": out_transform,
                     "dtype" : data_type
                     }
                    )

    out_name = "2m_REMA_mosaic_DEM"
    
    print("4/4")

with rasterio.open(output_path + out_name + "_dem.tif", "w", **out_meta) as dest:
    dest.write(out_image)

print("done")
    


# In[ ]:







# In[ ]:





# In[ ]:




