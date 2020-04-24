#!/usr/bin/env python
# coding: utf-8

# ## Exploring data visualization to look for subsetting
# 
# #### Credits
# * notebook by Jessica Scheick, derived from DEM example

# #### Setup
# ##### The Notebook was run on ICESat2 Hackweek 2019 pangeo image
# ##### For full functionality,
# - Please install [icepyx](https://github.com/icesat2py/icepyx), [topolib](https://github.com/ICESAT-2HackWeek/topohack), [contextily](https://github.com/darribas/contextily) using `git clone xxxxx`, `pip install -e .` workflow (see below; **you must restart your kernel after installing the packages**)
# - Download [NASA ASP](https://github.com/NeoGeographyToolkit/StereoPipeline) tar ball and unzip, we execute the commands from the notebook, using the path to the untared bin folder for the given commands.

# In[1]:


# get_ipython().run_cell_magic('bash', '', 'cd ~\n# git clone https://github.com/icesat2py/icepyx.git\n# git clone https://github.com/ICESAT-2HackWeek/topohack.git\n# git clone https://github.com/darribas/contextily.git\n\ncd contextily\npip install -e .\ncd ../topohack\npip install -e .\ncd ../icepyx\npip install -e .')


# In[1]:


#needs to be wherever icepyx, contextily, and topolib are installed in the previous step (ideally $HOME)
# %pwd

cd contextily
pip install -e .
cd ../topohack
pip install -e .
cd ../icepyx
pip install -e .


# #### ICESat-2 product being explored : [ATL08](https://nsidc.org/data/atl08)
# - Along track heights for canopy (land and vegitation) and  terrain
# - Terrain heights provided are aggregated over every 100 m along track interval, output contains "h_te_best_fit: height from best fit algorithm for all photons in the range", median height and others. Here we use h_te_best_fit.
# - See this preliminary introduction and quality assessment [paper](https://www.mdpi.com/2072-4292/11/14/1721) for more detail

# ### Import packages, including icepyx

# In[2]:

#from icepyx import is2class as ipd
import os
import shutil
import h5py
import xarray as xr 
# depedencies
import getpass
#from topolib.subsetDat import subsetBBox;
from topolib import icesat2_data 
import glob
import rasterio
from topolib import gda_lib
#from topolib import dwnldArctic
import numpy as np
import geopandas as gpd
from multiprocessing import Pool
import contextily as ctx
import pandas as pd
import matplotlib.pyplot as plt


# In[3]:


get_ipython().run_line_magic('load_ext', 'autoreload')
from icepyx import is2class as ipd
get_ipython().run_line_magic('autoreload', '2')
#in order to use "as ipd", you have to use autoreload 2, which will automatically reload any module not excluded by being imported with %aimport -[module]


# In[4]:


get_ipython().run_line_magic('cd', '~/software/icepyx/dev-notebooks')


# ## subset and non data objects

# In[6]:


# =============================================================================
#arr
# reg_a_poly = [(-64, 66), (-64, 72), (-55, 72), (-55, 66), (-64, 66)]

poly = gdf.iloc[0].geometry

# Simplify polygon. The larger the tolerance value, the more simplified the polygon.
poly = poly.simplify(0.05, preserve_topology=False)

# Orient counter-clockwise
poly = orient(poly, sign=1.0)

print(poly)

#Format dictionary to polygon coordinate pairs for CMR polygon filtering
polygon = ','.join([str(c) for xy in zip(*poly.exterior.coords.xy) for c in xy])
# =============================================================================


region_areg = ipd.Icesat2Data('ATL08', [-73.9, 10.7, -73.4, 11.1], ['2018-12-01','2019-09-01'],                           start_time='00:00:00', end_time='23:59:59')
#2019-01-04; 2019-01-06 works for subsetting


# In[7]:


region_asub = ipd.Icesat2Data('ATL08', [-73.9, 10.7, -73.4, 11.1], ['2018-12-01','2019-09-01'],                           start_time='00:00:00', end_time='23:59:59')
#2019-02-01; 2019-02-04 doesn't work for subsetting


# above: bounding box in Colombia
# below: shapefile in Antarctica

# In[7]:


region_areg = ipd.Icesat2Data('ATL06', '/home/jovyan/icepyx/doc/examples/supporting_files/data-access_PineIsland/glims_polygons.kml',                           ['2019-02-22','2019-02-28'],                           start_time='00:00:00', end_time='23:59:59')


# In[8]:


region_asub = ipd.Icesat2Data('ATL06', '/home/jovyan/icepyx/doc/examples/supporting_files/data-access_PineIsland/glims_polygons.kml',                           ['2019-02-22','2019-02-28'],                           start_time='00:00:00', end_time='23:59:59')


# In[18]:


region_areg=None
region_asub=None


# #### Log in to Earthdata

# In[9]:


earthdata_uid = 'Jessica.scheick'
email = 'jessica.scheick@maine.edu'
sessionr=region_areg.earthdata_login(earthdata_uid, email)
sessions=region_asub.earthdata_login(earthdata_uid, email)


# In[22]:


#search for available granules
region_areg.avail_granules()


# In[23]:


region_asub.avail_granules()


# In[21]:


print(region_areg.granule_info)
print(region_asub.granule_info)


# #### Place the order

# In[10]:


region_areg.order_granules(sessionr, subset=False)
#region_a.order_granules(session, verbose=True)


# In[16]:


region_asub.order_granules(sessions, subset=True, verbose=True)
#region_a.order_granules(session, verbose=True)


# #### Download the order

# In[17]:


wd = get_ipython().run_line_magic('pwd', '')
pathreg = wd + '/downloadreg'
pathsub = wd + '/downloadsub'


# In[18]:


region_areg.download_granules(sessionr, pathreg)


# In[19]:


region_asub.download_granules(sessions, pathsub)


# #### Clean up the download folder by removing individual order folders:

# In[22]:


#Clean up Outputs folder by removing individual granule folders 
path=pathsub
for root, dirs, files in os.walk(path, topdown=False):
    for file in files:
        try:
            shutil.move(os.path.join(root, file), path)
        except OSError:
            pass
        
for root, dirs, files in os.walk(path):
    for name in dirs:
        os.rmdir(os.path.join(root, name))


# ## Preprocess #2
# - Convert data into geopandas dataframe, which allows for doing basing geospatial opertaions

# In[23]:


get_ipython().run_line_magic('cd', '/home/jovyan/icepyx/dev-notebooks')


# In[24]:


# glob to list of files (run block of code creating wd and path variables if starting processing here)
ATL08_list = sorted(glob.glob(pathreg+'/*.h5'))
print(ATL08_list)


# In[25]:


# glob to list of files (run block of code creating wd and path variables if starting processing here)
ATL08_listsub = sorted(glob.glob(pathsub+'/*.h5'))
print(ATL08_listsub)


# ### Examine content of 1 ATLO8 hdf file

# In[27]:


# dict containing data entries to retrive (ATL08)
dataset_dict = {'land_segments':['delta_time','longitude','latitude','atl06_quality_summary','quality','terrain_flg'], 'land_segments/terrain':['h_te_best_fit']}


# In[38]:


#gda_lib.ATL08_to_dict(ATL08_list[0],dataset_dict)


# In[31]:


## the data can be converted to geopandas dataframe, see ATL08_2_gdf function in topolib gda_lib
temp_gdf = gda_lib.ATL08_2_gdf(ATL08_list[0],dataset_dict)


# In[29]:


## the data can be converted to geopandas dataframe, see ATL08_2_gdf function in topolib gda_lib
temp_gdfsub = gda_lib.ATL08_2_gdf(ATL08_listsub[0],dataset_dict)


# In[14]:


temp_gdf.head()


# In[32]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[41]:


temp_gdf.plot()
#plt.ylim(10.5,11.2)


# In[46]:


temp_gdfsub.plot()
#plt.ylim(10.5,11.2)


# In[42]:


temp_gdf.total_bounds


# In[43]:


region_areg.spatial_extent


# In[47]:


temp_gdfsub.total_bounds


# In[ ]:





# In[30]:


colombia_crs = {'init':'epsg:32618'}
plot_web = {'init':'epsg:3857'}


# In[15]:


temp_gdf.keys()


# In[28]:


gdf_list = [(gda_lib.ATL08_2_gdf(x,dataset_dict)) for x in ATL08_list]
gdf_colombia = gda_lib.concat_gdf(gdf_list)


# In[29]:


gdf_listsub = [(gda_lib.ATL08_2_gdf(x,dataset_dict)) for x in ATL08_listsub]
gdf_colombiasub = gda_lib.concat_gdf(gdf_listsub)


# # Plot Bounding box data (Colombia)
# - Visualise data footprints

# In[31]:


fig,ax = plt.subplots(figsize=(10,10))
temp_web = gdf_colombia.to_crs(plot_web)
clim = np.percentile(temp_web['h_te_best_fit'].values,(2,98))
temp_web.plot('h_te_best_fit',ax=ax,s=3,legend=True,cmap='inferno',vmin=clim[0],vmax=clim[1])
ctx.add_basemap(ax=ax)
ax.set_xticks([])
ax.set_yticks([])


# In[33]:


fig,ax = plt.subplots(figsize=(10,10))
temp_websub = gdf_colombiasub.to_crs(plot_web)
climsub = np.percentile(temp_websub['h_te_best_fit'].values,(2,98))
temp_websub.plot('h_te_best_fit',ax=ax,s=3,legend=True,cmap='inferno',vmin=climsub[0],vmax=climsub[1])
ctx.add_basemap(ax=ax)
ax.set_xticks([])
ax.set_yticks([])


# # Plot Polygon data (Antarctica)
# - Visualise data footprints

# ### Convert the list of hdf5 files into more familiar Pandas Dataframe

# In[34]:


# dict containing data entries to retrive (ATL06)
dataset_dict={'land_ice_segments':['atl06_quality_summary','delta_time','h_li','hli_sigma',           'latitude','longitude','segment_id','sigma_geo_h'], 'land_ice_segments/ground_track':['x_atc']}


# In[38]:


ant_crs = {'init':'epsg:3031'}
plot_web = {'init':'epsg:3857'}


# In[35]:


gdf_list = [(gda_lib.ATL06_2_gdf(x,dataset_dict)) for x in ATL08_list]
gdf_colombia = gda_lib.concat_gdf(gdf_list)


# In[36]:


gdf_listsub = [(gda_lib.ATL06_2_gdf(x,dataset_dict)) for x in ATL08_listsub]
gdf_colombiasub = gda_lib.concat_gdf(gdf_listsub)


# In[39]:


fig,ax = plt.subplots(figsize=(10,10))
temp_web = gdf_colombia.to_crs(plot_web)
clim = np.percentile(temp_web['x_atc'].values,(2,98))
temp_web.plot('x_atc',ax=ax,s=3,legend=True,cmap='inferno',vmin=clim[0],vmax=clim[1])
ctx.add_basemap(ax=ax)
ax.set_xticks([])
ax.set_yticks([])


# In[40]:


fig,ax = plt.subplots(figsize=(10,10))
temp_websub = gdf_colombiasub.to_crs(plot_web)
climsub = np.percentile(temp_websub['x_atc'].values,(2,98))
temp_websub.plot('x_atc',ax=ax,s=3,legend=True,cmap='inferno',vmin=climsub[0],vmax=climsub[1])
ctx.add_basemap(ax=ax)
ax.set_xticks([])
ax.set_yticks([])

