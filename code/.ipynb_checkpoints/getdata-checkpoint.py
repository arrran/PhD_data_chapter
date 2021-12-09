#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 09:21:58 2019

@author: whitefar
"""


class data:
    
     def load_gis(self,perimeter_file):
            """
            Input filetype must be readable by geopandas .read_file, and must be a linestring, or polygon. E.g. .shp or .gpx
            
            see import fiona; help(fiona.open) for info.
            """
            self.perimeter_file = perimeter_file
                
            #Return a GeoDataFrame object using geopandas
            self.perimeter_gdf = gpd.read_file(perimeter_file)
            
            if self.perimeter_gdf.iloc[0].geometry.type =='Polygon':
                self.geometry = self.perimeter_gdf.iloc[0].geometry
                self.perimeter = LineString(list(self.geometry.exterior.coords))
                self.perimeter_coords = list(self.geometry.exterior.coords)
                self.perimeter_array = np.array(list(self.geometry.exterior.coords))
                self.num_points =  len(list(self.geometry.exterior.coords))
            elif self.perimeter_gdf.iloc[0].geometry.type =='LineString':
                self.geometry = self.perimeter_gdf.iloc[0].geometry
                self.perimeter = self.geometry
                self.perimeter_coords = list(self.perimeter.coords)
                self.perimeter_array = np.array(list(self.perimeter.coords))
                self.num_points = len(list(self.perimeter.coords))
            elif self.perimeter_gdf.iloc[0].geometry.type == 'Point':
                self.perimeter = LineString([self.perimeter_gdf.iloc[i].geometry.coords[:][0] for i in range(len(self.perimeter_gdf))])
                self.perimeter_coords = list(self.perimeter.coords)
                self.perimeter_array = np.array(list(self.perimeter.coords))
                self.num_points = len(list(self.perimeter.coords))
            else:
                raise ValueError("Input must have Point, LineString or Polygon geometry")
            
            print(self.perimeter_gdf.head())
            self.perimeter_gdf.plot()
            
        
        #def load_txt(self,text_file):
        #    """
        #    """
            

