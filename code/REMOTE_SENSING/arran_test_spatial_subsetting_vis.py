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
from shapely.geometry.polygon import orient

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

shp_filepath = '~/study_area_buffer_geo.shp'



# #Return a GeoDataFrame object
# gdf = gpd.read_file(shp_filepath)

# poly = gdf.iloc[0].geometry

# # Simplify polygon. The larger the tolerance value, the more simplified the polygon.
# #poly = poly.simplify(0.05, preserve_topology=False)

# # Orient counter-clockwise
# poly = orient(poly, sign=1.0)

#Format dictionary to polygon coordinate pairs for CMR polygon filtering
# reg_a_poly = [xy for xy in zip(*poly.exterior.coords.xy)]

# reg_a_dates = ['2019-02-22','2019-02-28']
# =============================================================================

# region_a = ipd.Icesat2Data('ATL06', reg_a_poly, reg_a_dates)
# region_areg = ipd.Icesat2Data('ATL08', [-73.9, 10.7, -73.4, 11.1], ['2018-12-01','2019-09-01'],                           start_time='00:00:00', end_time='23:59:59')
#2019-01-04; 2019-01-06 works for subsetting


# In[7]:


# region_asub = ipd.Icesat2Data('ATL08', [-73.9, 10.7, -73.4, 11.1], ['2018-12-01','2019-09-01'],                           start_time='00:00:00', end_time='23:59:59')
#2019-02-01; 2019-02-04 doesn't work for subsetting


# above: bounding box in Colombia
# below: shapefile in Antarctica

# In[7]:


region_areg = ipd.Icesat2Data('ATL06', shp_filepath,['2019-02-22','2019-02-28'],start_time='00:00:00', end_time='23:59:59')


# In[8]:


# region_asub = ipd.Icesat2Data('ATL06', '/home/jovyan/icepyx/doc/examples/supporting_files/data-access_PineIsland/glims_polygons.kml',                           ['2019-02-22','2019-02-28'],                           start_time='00:00:00', end_time='23:59:59')


# In[18]:


# region_areg=None


# #### Log in to Earthdata

# In[9]:


earthdata_uid = 'arran'
email = 'arran.whiteford@vuw.ac.nz'
sessionr=region_areg.earthdata_login(earthdata_uid, email)


# In[22]:


#search for available granules
region_areg.avail_granules()



print(region_areg.granule_info)



# #### Place the order

# In[10]:


region_areg.order_granules(sessionr, subset=False)
#region_a.order_granules(session, verbose=True)


# In[16]:


# #### Download the order

# In[17]:


wd = get_ipython().run_line_magic('pwd', '')
pathreg = '~/Test1/'


# In[18]:


region_areg.download_granules(sessionr, pathreg)


# In[19]:


# #### Clean up the download folder by removing individual order folders:

# In[22]:


#Clean up Outputs folder by removing individual granule folders 
path=pathreg
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

# # In[23]:
# get_ipython().run_line_magic('pwd')

# get_ipython().run_line_magic('cd', '/Users/home/whitefar/PROJECTS/RADAR')


# In[24]:


# glob to list of files (run block of code creating wd and path variables if starting processing here)
ATL08_list = sorted(glob.glob(pathreg+'/*.h5'))
print(ATL08_list)


# In[25]:
# dict containing data entries to retrive (ATL08)


dataset_dict = {'land_segments':['delta_time','longitude','latitude','atl06_quality_summary','quality','terrain_flg'], 'land_segments/terrain':['h_te_best_fit']}

# ### Examine content of 1 ATLO8 hdf file

# In[27]:
# Earthdata Login credentials

# Setup a search session
search = IceSat2Data(user_id, password, variables=variables)

In [8]:

# Show available variables
search.show_variables()

# Show available formats
search.show_formats()


    
# =============================================================================
# 

def ATL08_to_dict(filename, dataset_dict):
    """
        Read selected datasets from an ATL06 file
        Input arguments:
            filename: ATl06 file to read
            dataset_dict: A dictinary describing the fields to be read
                    keys give the group names to be read,
                    entries are lists of datasets within the groups
        Output argument:
            D6: dictionary containing ATL06 data.  Each dataset in
                dataset_dict has its own entry in D6.  Each dataset
                in D6 contains a list of numpy arrays containing the
                data
    """

    D6=[]
    pairs=[1, 2, 3]
    beams=['l','r']
    # open the HDF5 file
    with h5py.File(filename,'r') as h5f:
        # loop over beam pairs
        for pair in pairs:
            # loop over beams
            for beam_ind, beam in enumerate(beams):
                # check if a beam exists, if not, skip it
                if '/gt%d%s/land_segments' % (pair, beam) not in h5f:
                    continue
                # loop over the groups in the dataset dictionary
                temp={}
                for group in dataset_dict.keys():
                    for dataset in dataset_dict[group]:
                        DS='/gt%d%s/%s/%s' % (pair, beam, group, dataset)
                        # since a dataset may not exist in a file, we're going to try to read it, and if it doesn't work, we'll move on to the next:
                        try:
                            temp[dataset]=np.array(h5f[DS])
                            # some parameters have a _FillValue attribute.  If it exists, use it to identify bad values, and set them to np.NaN
                            if '_FillValue' in h5f[DS].attrs:
                                fill_value=h5f[DS].attrs['_FillValue']
                                bad = temp[dataset]==fill_value
                                temp[dataset]=np.float64(temp[dataset])
                                temp[dataset][bad]=np.NaN
                        except KeyError as e:
                            pass
                if len(temp) > 0:
                    # it's sometimes convenient to have the beam and the pair as part of the output data structure: This is how we put them there.
                    #a = np.zeros_like(temp['h_te_best_fit'])
                    #print(a)
                    temp['pair']=np.zeros_like(temp['h_te_best_fit'])+pair
                    temp['beam']=np.zeros_like(temp['h_te_best_fit'])+beam_ind
                    #temp['filename']=filename
                    D6.append(temp)
    return D6


def ATL08_2_gdf(ATL06_fn,dataset_dict):
    """
    function to convert ATL06 hdf5 to geopandas dataframe, containing columns as passed in dataset dict
    Used Ben's ATL06_to_dict function
    """
    if ('latitude' in dataset_dict['land_segments']) != True:
        dataset_dict['land_segments'].append('latitude')
    if ('longitude' in dataset_dict['land_segments']) != True:
        dataset_dict['land_segments'].append('longitude')
    #use Ben's Scripts to convert to dict
    
    data_dict = ATL08_to_dict(ATL06_fn,dataset_dict)
    #this will give us 6 tracks
    i = 0
    
    df_final = []  
    
    for track in data_dict:
        #1 track
        #convert to datafrmae
        df = pd.DataFrame(track)
        df['p_b'] = str(track['pair'][0])+'_'+str(track['beam'][0])
        df['geometry'] = df.apply(point_covert,axis=1)
        if i==0:
            df_final = df.copy()
        else:
            df_final = df_final.append(df)
        i = i+1
    gdf_final = gpd.GeoDataFrame(df_final,geometry='geometry',crs={'init':'epsg:4326'})
    return gdf_final
# =============================================================================



# In[38]:


#gda_lib.ATL08_to_dict(ATL08_list[0],dataset_dict)


# In[31]:


## the data can be converted to geopandas dataframe, see ATL08_2_gdf function in topolib gda_lib
temp_gdf = ATL08_2_gdf(ATL08_list[0],dataset_dict)


# In[29]:




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

