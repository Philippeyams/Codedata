
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 15:04:44 2021

@author: Antoine, modifs: Carole le 28/04/21
"""
from osgeo import gdal,osr,ogr

#import pandas as pd
import numpy as np
import time
import sys
#import matplotlib.pyplot as plt
import os
import statistics

#J=['1j','2j']
name_DEM="C:/Users/Hp/Desktop/Dossier Stage HSM/Resultats SW2D et SW2D-DDP/LidarDrain1m.tif" 
InPath="C:/Users/Hp/Desktop/Dossier Stage HSM/Post traitement SW2D-DDP - Scenario2/Input/"                    # Les resultats SW2D
OutPath = "C:/Users/Hp/Desktop/Dossier Stage HSM/Post traitement SW2D-DDP - Scenario2/Output/"                 # DOSSIER DE SORTIE
MiscPath = "C:/Users/Hp/Desktop/Dossier Stage HSM/Post traitement SW2D-DDP - Scenario2/Misc/"                   
name_shp = "C:/Users/Hp/Desktop/Dossier Stage HSM/Resultats SW2D et SW2D-DDP/Resultats SW2D-DDP K=20/Modhydro7Cells.shp"

EPSG = 32630

print("Bonjour")


def cartes(FilePoros_resz):
    for fich in FilePoros_resz:
        print('Creating maps for ',fich,' ..............................')
        
        ZSimBrut = open(InPath + fich + ".txt", "r")
        MPcontents = ZSimBrut.readlines()
        MPr = []
        for x in MPcontents:
        	splitx = x.strip().split()
        	MPr.append(float(splitx[2]))
        ZSimBrut.close()
                
        #print("Opening SHP1: cells shapefile without Z info .........................")
        shpDriver = ogr.GetDriverByName("ESRI Shapefile")
        CellsSHP = shpDriver.Open(name_shp, 1)
        if CellsSHP == None:
            raise Exception("Unable to open SHP ....................................", name_shp)
        
        LayerSHP1 = CellsSHP.GetLayer()
        CellsCount = LayerSHP1.GetFeatureCount()
        
        #print("Getting spatial reference from SHP1 ..................................")
        
        SpatialReference = osr.SpatialReference()
        if (EPSG != 0) :
            SpatialReference.ImportFromEPSG(EPSG)
        else:
            spatialRef = LayerSHP1.GetSpatialRef() 
            SpatialReference.ImportFromWkt(spatialRef.ExportToWkt()) 		


        # ----------------------------------------------------------------------------
        # STEP 4: CREATE A NEW SHAPEFILE WITH Z INFO IN CELLS ------------------------	
        
        #print("Creating SHP2: new shapefile with Z info extracted from the results text file  .............................................................................")
        ZCellsSHP = shpDriver.CreateDataSource(MiscPath + fich + ".shp")
        LayerSHP2 = ZCellsSHP.CreateLayer('Cells',SpatialReference,geom_type=ogr.wkbMultiPolygon)		
        
        #print("Creating fields in SHP2 ..............................................")		
        FieldZSim = ogr.FieldDefn('ZSim', ogr.OFTReal)
        LayerSHP2.CreateField(FieldZSim)
        FieldID = ogr.FieldDefn('ID', ogr.OFTInteger)
        LayerSHP2.CreateField(FieldID)
        featureDefn = LayerSHP2.GetLayerDefn()
        #print('MP - creating the ZCells shapefile ...............................')
        
        # ----------------------------------------------------------------------------
        # STEP 5 : Fill THE NEW FIELD DB ----------------------------------------------
        
        i=0
        for elt in LayerSHP1:
            
            feature = ogr.Feature(featureDefn)
            feature.SetGeometry(elt.GetGeometryRef())
            i = i+1
            feature.SetField('ZSim', MPr[i-1])
            feature.SetField('ID', i)
            LayerSHP2.CreateFeature(feature)
        
    	
        # ----------------------------------------------------------------------------
        # STEP 6: READ THE DEM TO GET SPATIAL REFERENCE ------------------------------
        
        #print("Loading the DEM to get spatial reference .............................")
        try:
        	DEM = gdal.Open(name_DEM)
        except RuntimeError:
        	print("Unable to open ", name_DEM)
        	sys.exit(1)
        gt=DEM.GetGeoTransform()
        
        try:
        	DemBand=DEM.GetRasterBand(1)
        except RuntimeError:
        	print("Unable to read", name_DEM)
        	sys.exit(1)
    
        # STEP 7: CONVERTING THE SHAPEFILE INTO RASTER --------------------------------
        
        rasterDriver = gdal.GetDriverByName('GTiff')
        #print("Creating a new raster ................................................")
        # Create the destination data source
        
        RasZSimBrut = rasterDriver.Create(MiscPath + fich + "_Brut" + ".tif", 
                                         DEM.RasterXSize, 
                                         DEM.RasterYSize,
                                         1, #missed parameter (band)
                                         gdal.GDT_Float32)
        RasZSimBrut.SetGeoTransform(gt)
        RasZSimBrut.SetProjection(DEM.GetProjection())
        tempTile = RasZSimBrut.GetRasterBand(1)
        tempTile.SetNoDataValue(-9999)
        gdal.RasterizeLayer(RasZSimBrut, [1], LayerSHP2, None, options=['ATTRIBUTE=ZSim'])
        
        Z = tempTile.ReadAsArray()
    	
        # ----------------------------------------------------------------------------
        # STEP 8: FILTER THE Z VALUES -------------------------------------------------
    
        MNT = DemBand.ReadAsArray()
        #print(MNT)
        
        #GeotiffWrite(RasZSimBrut, "ZsimBrut.tif",gt,gdal.GDT_Float32,DEM.GetProjection())
        #print("Creating a new raster with Z values filtered .........................")
        
        
        # Create the raster with flood extent only
        RasFloodExtent = rasterDriver.Create(OutPath + "fem_" + fich + ".tif", 
                                         DEM.RasterXSize, 
                                         DEM.RasterYSize,
                                         1, #missed parameter (band)
                                         gdal.GDT_Byte)
        
        RasFloodExtent.SetGeoTransform(gt)
        RasFloodExtent.SetProjection(DEM.GetProjection())
        tempTile2 = RasFloodExtent.GetRasterBand(1)
        tempTile2.SetNoDataValue(-9999)	
        #print("Filtering the Z values ...............................................")
        
        #Set hmin ------------------------------------------------------------------
        hmin = 0.01								
        FloodExtent = np.where(Z>(MNT+hmin), 10.0, 0.0)
        FloodExtent[Z==-9999] = np.nan
        tempTile2.WriteArray(FloodExtent)
        
        ZRasFloodExtent = rasterDriver.Create(OutPath + "zfem_" + fich + ".tif", 
                                         DEM.RasterXSize, 
                                         DEM.RasterYSize,
                                         1, #missed parameter (band)
                                         gdal.GDT_Float32)
        
        ZRasFloodExtent.SetGeoTransform(gt)
        ZRasFloodExtent.SetProjection(DEM.GetProjection())
        ZtempTile2 = ZRasFloodExtent.GetRasterBand(1)
        ZtempTile2.SetNoDataValue(-9999)	
        #print("Filtering the Z values ...............................................")
        							
        ZFloodExtent = np.where(Z>(MNT+hmin), Z-MNT, 0.0)
        ZFloodExtent[Z==-9999] = np.nan
        ZFloodExtent[ZFloodExtent==0] = np.nan
        ZtempTile2.WriteArray(ZFloodExtent)	
        
    
        #------------------------------------------------------------------------------
        # -----------------------------------------------------------------------------
    
   
	


    def main():
        sys.stdout.write("Program done ..........................................")
		   
        if __name__ == '__main__':
            main()

    return

# ----------------------------------------------------------------------------
#FolderList = GetFolderList(Dossier)
ResFiles = []
for file in os.listdir(InPath):
    if file.endswith(".txt"):
        ResFiles.append(file[0:len(file)-4])

#FilePoros_resz=['alzette_z_j1h0.txt','alzette_z_j2h6.txt']
cartes(ResFiles)


