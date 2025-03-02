# -*- coding: utf-8 -*-
from osgeo import ogr,osr, gdal
import sys, struct
import os
import numpy as np
import Param
from shapely.geometry import shape, Point, LineString, Polygon
import time
import matplotlib.pyplot as plt

gdal.UseExceptions()

# Version: 2024-02 : loi -2 dans tous les cas, rivière ou pas

# TODO: gerer les interfaces de conditions aux limites

# Version: 2021-07 : law -2 (added in sw2d so that it is not modified)

# Version: 2021-04.2 : remove use of FILE Bathy 
#					+ remove river reach interfaces (thus BUFPlus BUFminus)
# 					test for the presence of RBpoly; get spatial ref from DEM

gdal.SetConfigOption('GTIFF_SRS_SOURCE', 'GEOKEYS')       #Pour fixer le probleme de la projection (YAMEOGO)

#creation des classes pour : les cellules, les noeuds, les interfaces
#on cree aussi une classe qui recupere les infos du maillage
class cell():
	def __init__(self):
		self.ID = 0
		self.IDSMS = 0
		self.Nnodes = 0
		self.nodesID = []
		self.poly = None
		self.area = 0
		self.interfaces = []
		self.center = [] 
		self.bounds = [] # xMin, xMax, yMin, yMax
		self.Type = "FP" 
		self.PorLaw = -2 # pour l'instant ne peut etre changé
		self.NbTab = 0
		self.zPor = []     
		self.zPhi = []
		self.PorVal = []
		self.SurfRatioInf = 0
		self.SurfRatioSup = 0


#
class node:
	def __init__(self, x, y, z):
		self.point = Point(x, y, z)
		self.ID = 0
		self.nbInt = 0
		self.interfaces = []

class interface:
	def __init__(self, Node1, Node2):
		self.line = LineString([Node1.point, Node2.point])
		self.ID = 0
		self.Node1 = Node1
		self.Node2 = Node2
		self.line = None
		self.celL = None
		self.celR = None
		self.Type = "FP"
		self.zPor = []
		self.zPhi = []
		self.PorVal = []
		

class Mesh:
	def __init__(self,name_2dm,name_RBshp,name_DEM):

		# read 2dM and define cells, nodes, interfaces

		
		self.listeCells = []
		self.listeNodes = []
		self.listeInterfaces = []

		self.nbCells = 0
		self.nbNodes = 0
		self.nbInterfaces = 0

		# todo: to be completed		
		self.checkparameters()

#Creation de 3 nouveaux fichiers shp en partant du .2dm
		name_shp = name_2dm[0:len(name_2dm)-4]+"Cells.shp"
		name_intershp = name_2dm[0:len(name_2dm)-4]+"Inter.shp"
		name_nodeshp = name_2dm[0:len(name_2dm)-4]+"Nodes.shp"

		# read cells and nodes from 2DM

		if Param.FIND_INTER == 1: 
			print("Read 2dm file and assign nodes elevation from DEM.....................")
			self.MeshFrom2dM(name_2dm,name_DEM)
		else:
			print("Read shp files .......................................................")
			self.MeshFromSHP(name_shp,name_intershp, name_nodeshp)
			
		
		# Read DEM
		try:
			imageDEM = gdal.Open(name_DEM)
		except RuntimeError:
			print("Unable to open ",name_DEM)
			sys.exit(1)
		gtDEM=imageDEM.GetGeoTransform()
		spatialRef = imageDEM.GetSpatialRef() 


		# initialise Cells with
		print("Assigning Cells and Interfaces attributes ..............................................")
		self.AssignAttributes(name_RBshp)
		print(".....................................................................done")
		
		#for c in self.listeCells:
		#	print(c.PorLaw)
		
		# create a shp with cells as polygons and law type
		if Param.FIND_INTER >=0 :
			print("Create the shapefiles for cells, nodes and interfaces ...................")
			self.MeshSHP(name_shp,spatialRef)
			self.InterSHP(name_intershp,spatialRef)
			self.NodeSHP(name_nodeshp,spatialRef)
		else:
			print("SHP files already exist")
		print(".....................................................................done")
		
		# get values from DEM
		print("Get values from DEM......................................................")
		self.TabFromDEM(name_2dm,name_DEM,name_RBshp)
		print(".....................................................................done")
		

	def checkparameters(self):  
		# check format
		if (Param.CPP == 1) :
			print("Creation of porosity files for SW2D CPP ======")
		if (Param.FOR == 1) :
			print("Creation of porosity files for SW2D Fortran ======")
		
        # check shape
    #verification du type de canal: trapézoidal, rect ou autre
		if (Param.SHAPE == 'trap') :
			print("Trapezoidal shape for riverbed")
		elif (Param.SHAPE == 'rect') :
			print("Rectangular shape for riverbed")
		else:
			print("Unknown shape for riverbed:", Param.SHAPE)




	# read 2dm file
	def MeshFrom2dM(self,name_2dm,name_DEM):
		
		# Read DEM to assign Z values -> TODO: a mettre ailleurs???
		try:
			imageDEM = gdal.Open(name_DEM)
		except RuntimeError:
			print("Unable to open ",name_DEM)
			sys.exit(1)
		gt=imageDEM.GetGeoTransform()

		try:
			rDEM=imageDEM.GetRasterBand(1)
		except RuntimeError:
			print("Unable to read DEM ",name_DEM)
			sys.exit(1)

		# read 2dM and define cells, nodes, interfaces
		SMS2dm = open(name_2dm,"r")
		line = SMS2dm.readline()
		Cells = []
		Nodes = []
		
		intcell = dict()
		listeInterfaces=[]

		Nc = 0
		Nn = 0

		for line in SMS2dm:
			words = line.split()


#parcours le 2DM ligne et par ligne et determine la nature des mailles: triangulaure, quadragulaire....
			#line = SMS2dm.readline()
			if (words[0]=="E3T") :
				Nc = Nc+1				
				C = cell()
				C.ID = Nc
				C.IDSMS = int(words[1])
				C.Nnodes = 3
				w2=int(words[2]); w3=int(words[3]); w4=int(words[4])
				C.nodesID = [w2,w3,w4,w2]
				intcell[Nc] = { 1:(w2,w3),2:(w3,w4),3:(w4,w2) }      #Un dictionnaire qui stocke les numeros des noeurs pour chaque maille
				Cells.append(C)


			if (words[0]=="E4Q") : 
				Nc = Nc+1				
				C = cell()
				C.ID = Nc
				C.IDSMS = int(words[1])
				C.Nnodes = 4
				w2=int(words[2]); w3=int(words[3]); w4=int(words[4]); w5=int(words[5])
				C.nodesID = [w2,w3,w4,w5,w2]
				intcell[Nc] = {1:(w2,w3),2:(w3,w4),3:(w4,w5),4:(w5,w2)}

				Cells.append(C)
				

			if (words[0]=="E8Q"):
				Nc = Nc+1				
				C = cell()
				C.ID = Nc
				C.IDSMS = int(words[1])
				C.Nnodes = 4
				w2=int(words[2]); w4=int(words[4]); w6=int(words[6]); w8=int(words[8])
				C.nodesID = [w2,w4,w6,w8,w2]
				intcell[Nc] = {1:(w2,w4),2:(w4,w6),3:(w6,w8),4:(w8,w2)}
				Cells.append(C)		
				

			if (words[0] =="ND"):
				Nn=Nn+1
				X = int((float(words[2]) - gt[0])/gt[1])
				Y = int((float(words[3]) - gt[3])/gt[5])

				if (Y >= imageDEM.RasterYSize or X >= imageDEM.RasterXSize or Y < 0 or X < 0):
					print(" exception x y not in DEM", X, Y)
					Z=0
				else: 
					z=rDEM.ReadAsArray(X,Y,1,1)
					if (z is None): 	
						print(" exception x y not in DEM", X, Y)
						Z=0
					elif (float(z[0][0]) < -100):
						print(" z = ", float(z[0][0]), "replaced by 9999")
						Z = 9999
					else :
						Z = float(z[0][0])
						
				N = node(float(words[2]),float(words[3]),Z)
				N.ID = int(words[1])
				Nodes.append(N)
				

		
		self.nbCells = Nc
		self.nbNodes = Nn
		self.listeCells = Cells
		self.listeNodes = Nodes
		print ("Number of cells in mesh: ", self.nbCells)
		print ("Number of nodes in mesh: ", self.nbNodes)

		
		print("........................................................... 2dm File read")
				
		

		print("Looking for interfaces...................................................")
		nbInt = 0
		Interfaces = dict()

		for k,v in intcell.items(): #every cells are stored here as a set of couples of nodes. k = cell.ID
			if (k % 100) == 0 : 
				print(".... cell ", k)
			
			for i in v.values() : # for each interface (couple of nodes) of c
	#			inter2 = interface(self.listeNodes[i[1]-1],self.listeNodes[i[0]-1])
				i2 = (i[1],i[0])

				# Create interface if it does not exist. If it does it is necessarily inter2
				if (i2 not in Interfaces):
					nbInt = nbInt+1
					inter = interface(self.listeNodes[i[0]-1],self.listeNodes[i[1]-1])
					inter.ID = nbInt
					inter.celL = k
					self.listeCells[k-1].interfaces.append(inter)
					self.listeInterfaces.append(inter)
					Interfaces[i] = inter
				else : # i2 already exists => k is inter2.celR
					Interfaces[i2].celR=k
							
		self.nbInterfaces = nbInt
		
		# TODO: faire une version pour les interfaces cpp: numérotées à partir des noeuds et pas des cellules
        # voir swMesh.cpp ligne 122 et PorosGenerator2023.py

		print(".....................................................................done")
		print ("Number of interfaces in mesh: ", self.nbInterfaces)




	def MeshFromSHP(self,name_shp,name_intershp, name_nodeshp):

		Cells = []
		Nodes = []
		Interfaces = []

		print("Open cells shapefile ...................", name_shp)
		driver = ogr.GetDriverByName("ESRI Shapefile")
		dataSource = driver.Open(name_shp, 1)
		layerCells = dataSource.GetLayer()
		Nc = 0
		for cc in layerCells:
			Nc = Nc + 1
			C = cell()
			C.ID = Nc
			C.Nnodes = cc.GetField("Nnodes")
			if (C.Nnodes == 3):
				C.nodesID = [cc.GetField("Node1"),cc.GetField("Node2"),cc.GetField("Node3"),cc.GetField("Node1")]
			elif (C.Nnodes == 4):
				C.nodesID = [cc.GetField("Node1"),cc.GetField("Node2"),cc.GetField("Node3"),cc.GetField("Node4"),cc.GetField("Node1")]
			else : 
				print("Problem with number of nodes in cell", C.ID)
				return

			Cells.append(C)

		self.nbCells = Nc
		print ("number of cells in mesh:", Nc)


		print("Open nodes shapefile ...................", name_nodeshp)
		driver = ogr.GetDriverByName("ESRI Shapefile")
		dataSource = driver.Open(name_nodeshp, 1)
		layerNodes = dataSource.GetLayer()
		Nn = 0
		for n in layerNodes:
			point = n.GetGeometryRef()
			Nn = Nn + 1
			Nod = node(point.GetPoint()[0],point.GetPoint()[1],point.GetPoint()[2])
			Nod.ID = Nn
			Nodes.append(Nod)

		self.Nnodes = Nn
		print ("Number of nodes in mesh: ", Nn)


		print("Open interfaces shapefile...............", name_intershp)
		driver = ogr.GetDriverByName("ESRI Shapefile")
		dataSource = driver.Open(name_intershp, 1)
		layerInter = dataSource.GetLayer()
		nbInt = 0
		for i in layerInter:
			nbInt = nbInt + 1
			line = i.GetGeometryRef()

			P1 = node(line.GetPoint(0)[0],line.GetPoint(0)[1],i.GetField("Z1"))
			P2 = node(line.GetPoint(1)[0],line.GetPoint(1)[1],i.GetField("Z2"))
			inter = interface(P1,P2)
			inter.ID = nbInt

			inter.celL = i.GetField("CelL")
			if inter.celL == 0 : inter.celL = None
			inter.celR = i.GetField("CelR")
			if inter.celR == 0 : inter.celR = None
			
			Interfaces.append(inter)

		self.nbInterfaces = nbInt

		self.listeCells = Cells
		self.listeNodes = Nodes
		self.listeInterfaces = Interfaces
		
		print ("Number of interfaces in mesh: ", nbInt)

	

#-----------------------------------------------				

	# Cells initialisation
	def AssignAttributes(self,name_RBshp):

		# get RiverBed polygon 
		RB = []
		if name_RBshp != "" :
			driver = ogr.GetDriverByName("ESRI Shapefile")
			try:
				dataSource = driver.Open(name_RBshp, 1)
				layerRiverBed = dataSource.GetLayer()
			except RuntimeError:
				print("Unable to open ",name_RBshp)
				sys.exit(1)
			
			for elt in layerRiverBed:
				RB.append(elt.GetGeometryRef())
		
		
		for c in self.listeCells:

			ring = ogr.Geometry(ogr.wkbLinearRing)
			for p in c.nodesID: # listeNodes includes the first point twice to close the ring
				pp = self.listeNodes[p-1]
				ring.AddPoint(pp.point.x,pp.point.y,pp.point.z)
			
			# bounds of the cell polygon	
			c.bounds = ring.GetEnvelope()
			

			c.poly = ogr.Geometry(ogr.wkbPolygon)
			c.poly.AddGeometry(ring)
			c.area = c.poly.GetArea()
			
			c.center = c.poly.Centroid()
			pt = ogr.Geometry(ogr.wkbPoint)
			pt.SetPoint(0,c.center.GetX(),c.center.GetY())
			
			for polyRB in RB:

				# TODO: ne fonctionne pas s'il y a plusieurs elements dans la couche RB
				polybuf = []
				if polyRB.Intersect(c.poly):
					inter = polyRB.Intersection(c.poly)
					area = inter.GetArea()
					# TODO: pose pb s'il y a un drain étroit dans une grande cellule => parametrer ou trouver autre chose
					if (area/c.area) > 0.05:
						c.SurfRatioSup = int(area/c.area*100)
						c.Type = "RB"
						
						polybuf=polyRB.Buffer(Param.BUFBanks)
						interbuf = polybuf.Intersection(c.poly)
						area = interbuf.GetArea()
						c.SurfRatioInf = int(area/c.area*100)

		# interface
		for i in self.listeInterfaces :	
			
			wkt = "LINESTRING ("+str(i.Node1.point.x)+" "+str(i.Node1.point.y)+","+str(i.Node2.point.x)+" "+str(i.Node2.point.y)+")"
			ring = ogr.CreateGeometryFromWkt(wkt)
			i.line = ring

			if (i.celL is None):
				# @bug: quand on ecrit cL=cR,il s'agit du meme objet en memoire
				# => si on change l'identifiant de cR ca change celui de cL
				# a priori aucune importance ici, mais a verifier.
				cR = self.listeCells[i.celR-1] 
				cL = cR # on prend les valeurs de l'autre cellule
				#cL.ID = -cR.ID
			if (i.celR is None):
				cL = self.listeCells[i.celL-1]
				cR = cL
				#cR.ID = -cL.ID

			if (i.celL is not None) and (i.celR is not None) :
				cL = self.listeCells[i.celL-1]
				cR = self.listeCells[i.celR-1]

			if cR.Type == "RB" or cL.Type == "RB":
				for polyRB in RB: # section en travers du lit mineur
					inter = ring.Intersection(polyRB)
					if str(inter) == "LINESTRING EMPTY" :
						i.Type = "R"
					else :
						i.Type = "RB"
					#print("inter", i.ID, i.Type, inter)
					

#-----------------------------------------------				


	# create a shp file from 2dm including porosity law type
	def MeshSHP(self,name_shp,spatialRef):
		
		# Create cells shp	
		SR = osr.SpatialReference()
		if (str(Param.EPSG) != "") :
			SR.ImportFromWkt(Param.EPSG)
		else:
			SR.ImportFromWkt(spatialRef.ExportToWkt()) 		
		
		driver = ogr.GetDriverByName("ESRI Shapefile")
		dataSource = driver.CreateDataSource(name_shp)
		layerCells = dataSource.CreateLayer('Cells',SR,geom_type=ogr.wkbMultiPolygon)
		
		# Attributes definition
		field_ID = ogr.FieldDefn('ID', ogr.OFTInteger)
		layerCells.CreateField(field_ID)
		field_IDSMS = ogr.FieldDefn('ID SMS', ogr.OFTInteger) 
		layerCells.CreateField(field_IDSMS)
		field_CenterX = ogr.FieldDefn('CenterX', ogr.OFTReal)
		layerCells.CreateField(field_CenterX)
		field_CenterY = ogr.FieldDefn('CenterY', ogr.OFTReal)
		layerCells.CreateField(field_CenterY)
		field_Area = ogr.FieldDefn('AREA', ogr.OFTInteger)
		layerCells.CreateField(field_Area)
		field_PorLaw = ogr.FieldDefn('PORLAW', ogr.OFTInteger)
		layerCells.CreateField(field_PorLaw)
		field_Type = ogr.FieldDefn('Type', ogr.OFTString)
		layerCells.CreateField(field_Type)
		field_RatioInf = ogr.FieldDefn('RATIO INF', ogr.OFTInteger)
		layerCells.CreateField(field_RatioInf)
		field_RatioSup = ogr.FieldDefn('RATIO SUP', ogr.OFTInteger)
		layerCells.CreateField(field_RatioSup)
		field_Nnodes = ogr.FieldDefn('Nnodes', ogr.OFTInteger)
		layerCells.CreateField(field_Nnodes)
		field_Node1 = ogr.FieldDefn('Node1', ogr.OFTInteger)
		layerCells.CreateField(field_Node1)
		field_Node2 = ogr.FieldDefn('Node2', ogr.OFTInteger)
		layerCells.CreateField(field_Node2)
		field_Node3 = ogr.FieldDefn('Node3', ogr.OFTInteger)
		layerCells.CreateField(field_Node3)
		field_Node4 = ogr.FieldDefn('Node4', ogr.OFTInteger)
		layerCells.CreateField(field_Node4)
		
		featureDefn = layerCells.GetLayerDefn()

		# FILL THE DATABASE
		for c in self.listeCells:
			feature = ogr.Feature(featureDefn)
			feature.SetGeometry(c.poly)
			feature.SetField('ID',c.ID)
			feature.SetField('ID SMS',c.IDSMS)
			feature.SetField('CenterX',c.center.GetX())			
			feature.SetField('CenterY',c.center.GetY())
			feature.SetField('AREA',c.area)
			feature.SetField('PORLAW',c.PorLaw)
			feature.SetField('Type',c.Type)
			feature.SetField('RATIO INF',c.SurfRatioInf)
			feature.SetField('RATIO SUP',c.SurfRatioSup)
			feature.SetField('Nnodes',c.Nnodes)

			feature.SetField('Node1',c.nodesID[0])
			feature.SetField('Node2',c.nodesID[1])
			feature.SetField('Node3',c.nodesID[2])
			feature.SetField('Node4',-999)
			if c.Nnodes == 4:
				feature.SetField('Node4',c.nodesID[3])
			layerCells.CreateFeature(feature)
#-----------------------------------------------


	def InterSHP(self,name_shp,spatialRef):

		# Create inter shp	
		SR = osr.SpatialReference()
		if (str(Param.EPSG) != "") :
			SR.ImportFromWkt(Param.EPSG)
		else:
			SR.ImportFromWkt(spatialRef.ExportToWkt()) 		

		driver = ogr.GetDriverByName("ESRI Shapefile")
		dataSource = driver.CreateDataSource(name_shp)
		layerInter = dataSource.CreateLayer('Interfaces',SR,geom_type=ogr.wkbLineString)
		
		# Attributes definition
		field_ID = ogr.FieldDefn('ID', ogr.OFTInteger)
		layerInter.CreateField(field_ID)
		field_Z1 = ogr.FieldDefn('Z1', ogr.OFTReal)
		layerInter.CreateField(field_Z1)
		field_Z2 = ogr.FieldDefn('Z2', ogr.OFTReal)
		layerInter.CreateField(field_Z2)
		
		field_PorLaw = ogr.FieldDefn('Type', ogr.OFTString)
		layerInter.CreateField(field_PorLaw)
		field_CelL = ogr.FieldDefn('CelL', ogr.OFTInteger)
		layerInter.CreateField(field_CelL)
		field_CelR = ogr.FieldDefn('CelR', ogr.OFTInteger)
		layerInter.CreateField(field_CelR)
				
		featureDefn = layerInter.GetLayerDefn()
		
		# FILL THE DATABASE
		for i in self.listeInterfaces:
			feature = ogr.Feature(featureDefn)
			feature.SetGeometry(i.line)
			feature.SetField('ID',i.ID)
			feature.SetField('Z1',i.Node1.point.z)
			feature.SetField('Z2',i.Node2.point.z)
			feature.SetField('Type',i.Type)
			feature.SetField('CelL',i.celL)
			feature.SetField('CelR',i.celR)
			
			layerInter.CreateFeature(feature)
#-----------------------------------------------

	def NodeSHP(self,name_shp,spatialRef):

		# Create node shp
		SR = osr.SpatialReference()
		if (str(Param.EPSG) != "") :
			SR.ImportFromWkt(Param.EPSG)
		else:
			SR.ImportFromWkt(spatialRef.ExportToWkt()) 		

		driver = ogr.GetDriverByName("ESRI Shapefile")
		dataSource = driver.CreateDataSource(name_shp)
		layerNode = dataSource.CreateLayer('Nodes',SR,geom_type=ogr.wkbPoint)
		
		# Attributes definition
		field_ID = ogr.FieldDefn('ID', ogr.OFTInteger)
		field_Z = ogr.FieldDefn('Z', ogr.OFTReal)
		
		layerNode.CreateField(field_ID)
		layerNode.CreateField(field_Z)
				
		featureDefn = layerNode.GetLayerDefn()
		
		# FILL THE DATABASE
		for n in self.listeNodes:
			point = ogr.Geometry(ogr.wkbPoint)
			point.AddPoint(n.point.x, n.point.y)
			feature = ogr.Feature(featureDefn)
			feature.SetGeometry(point)
			feature.SetField('ID',n.ID)
			feature.SetField('Z',n.point.z)

			layerNode.CreateFeature(feature)




	def TabFromDEM(self,name_2dm,name_DEM,name_RBshp):

		# Read DEM
		try:
			imageDEM = gdal.Open(name_DEM)
		except RuntimeError:
			print("Unable to open ",name_DEM)
			sys.exit(1)
		gtDEM=imageDEM.GetGeoTransform()

		try:
			rDEM=imageDEM.GetRasterBand(1)
		except RuntimeError:
			print("Unable to read DEM ",name_DEM)
			sys.exit(1)


		RB = []
		if name_RBshp != "" :
			# get RiverBed polygon 
			driver = ogr.GetDriverByName("ESRI Shapefile")
			dataSource = driver.Open(name_RBshp, 1)
			layerRiverBed = dataSource.GetLayer()
			for elt in layerRiverBed:
				RB.append(elt.GetGeometryRef())

		#--------------------------------
		# create PhiW files for Fortran and CPP
		#--------------------------------

		if (Param.CPP == 1) :
			name_phiW_CPP = Param.OUTPUTREP + "DDP_cell_porosity_map.txt"
			fcpp = open(name_phiW_CPP,"w")
			print("Unif	0",file=fcpp)
			print("NtabMax	", Param.NbTABW, file=fcpp)
			print("==== Default", file=fcpp)
			print("Type	Parameters(type)", file=fcpp)
			print("0	1	0", file=fcpp)
			print("====		Distrib", file=fcpp)
			print("Type 	Parameters(type)", file=fcpp)


		if (Param.FOR == 1) :
			name_phiW_FOR = Param.OUTPUTREP + name_2dm[name_2dm.rfind("/")+1:len(name_2dm)-4]+"_phiW.in"
			ffor = open(name_phiW_FOR,"w")
			print(self.nbCells,file=ffor)

			
		
		for c in self.listeCells:
			if (c.ID % 100) == 0 : 
				print(".... cell ", c.ID)
			listeZ = []
			bound = c.bounds # bounds = (minx, maxx, miny, maxy)
			# bounds in pixels BegX, BegY, Nx, Ny
			BegX = int((c.bounds[0] - gtDEM[0])/gtDEM[1])
			BegY = int((c.bounds[2] - gtDEM[3])/gtDEM[5])
			Nx = int((c.bounds[1]-c.bounds[0])/gtDEM[1])
			Ny = int((c.bounds[3]-c.bounds[2])/(-gtDEM[5]))

			#print("Cellule", c.ID)
			z2 = 0
			for x in range(0,Nx-1,Param.STEPDTM):

				for y in range(0,Ny-1,Param.STEPDTM):
					pt = ogr.Geometry(ogr.wkbPoint)
					pt.SetPoint(0,bound[0]+x*gtDEM[1]+gtDEM[1]/2,bound[2]+y*(-gtDEM[5])-gtDEM[5]/2)
					
					if c.poly.Intersect(pt) :
						# Pt in pixels: (begx+x, begy+y)
						if (BegY-y >= imageDEM.RasterYSize or 
							BegX+x >= imageDEM.RasterXSize or 
							BegY-y < 0 or BegX+x < 0):
							print(" NONE type exception (x, y)=",BegX+x,", ",BegY-y," not in DEM. Cell", c.ID)

						else: 
							z=rDEM.ReadAsArray( BegX+x,BegY-y,1,1)
							if (z is None): 
								print(" NONE type exception (x, y)=",BegX+x,", ",BegY-y," not in DEM. Cell", c.ID)
							elif (np.isnan(float(z[0][0]))):
								print(" NAN type exception (x, y)=",BegX+x,", ",BegY-y," not in DEM. Cell", c.ID)
							elif (float(z[0][0]) < -100):
								print(" z = ", float(z[0][0]), "in cell", c.ID)
							else:
								z = round(float(z[0][0]),2)
								listeZ.append(z)
							if Param.z2 == "banks": 
								for polyRB in RB:
									polybuf=polyRB.Buffer(gtDEM[1]*2.) # augmentation de la rivière pour chercher le max sur les berges
									if polybuf.Intersect(pt):
										if z > z2:
											z2 = z
			
			#print("Nb points dans la cellule", c.ID, "=", len(listeZ))
			if (len(listeZ)==0):
				print(" No value in cell", c.ID, " please correct mesh. TODO: add value of neighbouring cell")
				listeZ.append(9999)
				#sys.exit()


			# --- Law 0
			# if c.PorLaw == 0:
			# 	a = sorted(listeZ)
			#  	# keep a value each NbTABW 
			# 	n = int(len(a)/Param.NbTABW)
			# 	b = [0,Param.NbTABW]
			# 	for i in range(Param.NbTABW):
			# 		b.append(a[i*n])
			# 	c.PorVal = b[2:]
			# 	print(*b,sep=' ', file=fi)

			# A la place de faire le sort, on peut utiliser directement les quantiles
			# PB: ce n'est pas la philosophie utilisée dans sw2d => les valeurs sont modifiées
			#	phi = np.arange(0.,1,1./Param.NbTABW)
			#	b = [0,Param.NbTABW]
			#	zquant = np.quantile(listeZ,phi)
			#	for i  in range(Param.NbTABW):
			#		b.append(zquant[i])
							

			
			# --- Law -2
			# PB avec les lois 0 et -1 dans le code sw2d: il modifie ce qu'on lui donne !
			# création d'une nouvelle loi -2
			if c.Type == "FP":

				sortZ = sorted(listeZ)
				size = len(sortZ)
				if (size < Param.NbTABW): 
					print("not enough z in cell", c.ID, "nb tab for this cell=",len(listeZ))
					c.NbTab = size
				else: 
					c.NbTab = Param.NbTABW 
				zmax =sortZ[len(sortZ)-1]; zmin = sortZ[0] 
				dz = (zmax-zmin)/c.NbTab
				if dz != 0:
					z_tab = np.arange(zmin,zmax,dz) # TODO: do we need to keep the maximum?
					ranks = np.searchsorted(sortZ,z_tab,side="left")
				else:
					z_tab = zmax+np.zeros(size)
					ranks = [n for n in range(size+1)]

				b = [-2,c.NbTab]
				
				for i in range(c.NbTab-1):
					b.append(round(ranks[i+1]/size,4))
					b.append(round(z_tab[i],4))
				b.append(1.0000)
				b.append(round(z_tab[c.NbTab-1],4))
				

				c.PorVal = b[2:]

				if (Param.FOR == 1) :
					print(*b,sep=' ', file=ffor)
				if (Param.CPP == 1) :
					print(*b,sep='\t', file=fcpp)

				


			# --- Law 3
			elif c.Type == "RB": # within riverbed
				#print ("Riverbed cell: ",c.ID)
				z1 = min(listeZ)


				# rectangle z1, z2, phi1=phi2, phi3=1
				# z1 = min bathy 
				# z2 = moyenne des z de la maille
				# phi1 = largeur au radier fixée
				
				if Param.SHAPE == 'rect':
					# z1 already defined
					if Param.Phi1 > 0:
						phi1 = Param.phi1
					else:
						phi1 = c.SurfRatioSup/100.
					phi2=phi1
				
				if Param.SHAPE == 'trap':
					if Param.Phi1 > 0: 
						phi1 = Param.Phi1
					else:
						phi1 = c.SurfRatioInf/100.

					if Param.Phi2 > 0: 
						phi2 = Param.Phi2
					else:
						phi2 = c.SurfRatioSup/100.

				phi3 = 1
				
				# if Param.z2 == "banks": z2=z2

				if Param.z2 == 'max':
					z2 = round(max(listeZ),4)
				elif Param.z2 == 'mean':
					z2 = round((max(listeZ)+min(listeZ))/2.,4)
				# else:
				# 	print("Unknown option for parameter z2. Should be max or mean")

				if (Param.FOR == 1) :
					print("-2 3 {0} {1} {2} {3} {4} {5}".format(phi1, z1, phi2, z2, phi3, z2+0.01), file=ffor)
				
				if (Param.CPP == 1) :
					#print("3\t{0}\t{1}\t{2}\t{3}\t{4}".format(z1,z2,phi1,phi2,phi3), file = fcpp)				
					print("-2\t3\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(phi1, z1, phi2, z2, phi3, z2+0.01), file=fcpp)
				
				c.PorVal = [phi1, z1, phi2, z2, phi3, z2+0.01]


				#c.PorVal = [z1,z2,phi1,phi2,phi3]
				

				# trapeze z1 = min bathy, z2 =max cellule, phi1= largeur puis linéaire entre phi1 et phi2, phi2=phi3=1
				# z1 sera calibré ensuite
				# 2 choix possibles pour phi1:  a) phi1 calculé à partir d'une pente de berger fixée (45 par exemple)
				#								b) largeur au radier fixée
				
				
				


		#--------------------------------
		# create PhiG file
		#--------------------------------


		if (Param.CPP == 1) :
			fcpp.close()
			name_phiG_CPP = Param.OUTPUTREP + "DDP_edge_porosity_map.txt"
			fcpp = open(name_phiG_CPP,"w")
			print("Unif	0",file=fcpp)
			print("NtabMax	", Param.NbTABG, file=fcpp)
			print("==== Default", file=fcpp)
			print("Type	Parameters(type)", file=fcpp)
			print("0	1	0", file=fcpp)
			print("====		Distrib", file=fcpp)
			print("Type 	Parameters(type)", file=fcpp)


		if (Param.FOR == 1) :
			ffor.close()
			name_phiG_FOR = Param.OUTPUTREP + name_2dm[name_2dm.rfind("/")+1:len(name_2dm)-4]+"_phiG.in"
			ffor = open(name_phiG_FOR,"w")
			print(self.nbInterfaces,file=ffor)

		for i in self.listeInterfaces:
			
			iBC = False
			if (i.celL is None):
				# @bug: quand on ecrit cL=cR,il s'agit du meme objet en memoire
				# => si on change l'identifiant de cR ca change celui de cL
				# a priori aucune importance ici, mais a verifier.
				cR = self.listeCells[i.celR-1] 
				cL = cR # on prend les valeurs de l'autre cellule
				#cL.ID = -cR.ID
				iBC = True
			if (i.celR is None):
				cL = self.listeCells[i.celL-1]
				cR = cL
				iBC = True
				#cR.ID = -cL.ID

			if (i.celL is not None) and (i.celR is not None) :
				cL = self.listeCells[i.celL-1]
				cR = self.listeCells[i.celR-1]

			# check coherence between flood plain and river bed cells
			if (iBC == False):
				if (cL.Type == "FP") and (cR.Type == "FP"):
					if (i.Type != "FP"):
						print("WARNING: non floodplain type interface between two floodplain cells:")
						print("Interface={0}  -  CelL={1}  -  CelR={2}".format(i.ID, cL.ID, cR.ID))
				if ((cL.Type == "FP") and (cR.Type == "RB")):
					if (i.Type != "R"):
						print("WARNING: non reach type interface between floodplain and riverberd cells:")
						print("Interface={0}  -  CelL={1}  -  CelR={2} (Id CPP)".format(i.ID, cL.ID, cR.ID))
				if ((cL.Type == "RB") and (cR.Type == "FP")):
					if (i.Type != "R"):
						print("WARNING: non reach type interface between floodplain and riverberd cells:")
						print("Interface={0}  -  CelL={1}  -  CelR={2}".format(i.ID, cL.ID, cR.ID))
				if ((cL.Type == "RB") and (cR.Type == "RB")):
					if (i.Type != "RB"):
						print("WARNING: non riverbed type interface between two riverbed cells:")
						print("Interface={0}  -  CelL={1}  -  CelR={2}".format(i.ID, cL.ID, cR.ID))


			# take the min phi values of cL and cR 
			if (cL.Type == "FP") and (cR.Type == "FP"):
				NTab = max(cL.NbTab,cR.NbTab)
				b = [-2, NTab]
				if cL.NbTab < cR.NbTab: 
					c1=cL; c2=cR
				else: c1=cR;c2=cL
				
				j=0
				for n in range(NTab):
					if (c1.PorVal[j*2]<c2.PorVal[n*2]): 
						b.append(c1.PorVal[j*2]); j+1
						b.append(max(c1.PorVal[j*2+1],c2.PorVal[n*2+1]))
					else: 
						b.append(c2.PorVal[n*2])
						b.append(max(c1.PorVal[j*2+1],c2.PorVal[n*2+1]))
					#b.append(min(cL.PorVal[n*2],cR.PorVal[n*2]))
					#b.append(max(cL.PorVal[n*2+1],cR.PorVal[n*2+1]))
				
			# TODO: remettre ici ce qu'il y avait dans le code 2022??? ou garder ce qu'il y a pour le CPP???
			if ((cL.Type=="FP") and (cR.Type=="RB")):
				if (i.Type!="R"):  
					print("WARNING: non reach type interface between floodplain and riverberd cells:")
					print("Interface={0}  -  CelL={1}  -  CelR={2}".format(i.ID,cL.ID,cR.ID))
				
				NTab = max(cL.NbTab,cR.NbTab)
				for n in range(NTab):
					j = 0
					# pour chaque z de cellule de plaine, on cherche le z rivière correspondant
					while ((j<3) and  (cL.PorVal[n * 2 +1] > cR.PorVal[j*2+1])) : # on passe au dessus de la cellule riviere
						# keep the maximum z and minimum phi to avoid pb in cpp -> find a better way?
						zinter = max(cL.PorVal[n*2+1],cR.PorVal[j*2+1])
						phiinter = min(cL.PorVal[n * 2],cR.PorVal[j*2])
						j=j+1
					b.append(phiinter)  
					b.append(zinter)

			if ((cL.PorLaw=="RB") and (cR.PorLaw=="FP")):
				NTab = max(cL.NbTab,cR.NbTab)
				for n in range(NTab):
					j = 0
					# pour chaque z de cellule de plaine, on cherche le z rivière correspondant						
					while ((j<3) and (cR.PorVal[n * 2 +1] > cL.PorVal[j*2+1])) : # on passe au dessus de la cellule riviere
						# keep the maximum z and minimum phi to avoid pb in cpp -> find a better way?
						zinter = max(cR.PorVal[n*2+1],cL.PorVal[j*2+1])
						phiinter = min(cR.PorVal[n*2],cL.PorVal[j*2])
						j=j+1
					b.append(phiinter)  
					b.append(zinter)

			if ((cL.PorLaw == "RB") and (cR.PorLaw=="RB")):
				if (i.Type != "RB"):
					print("WARNING: non riverbed type interface between two riverbed cells:")
					print("Interface={0}  -  CelL={1}  -  CelR={2}".format(i.ID,cL.ID,cR.ID))

				b = [-2, 3]
				for n in range(3):
					b.append(min(cL.PorVal[n * 2], cR.PorVal[n * 2]))  
					b.append(max(cL.PorVal[n * 2 + 1], cR.PorVal[n * 2 + 1]))
				
			
			i.PorVal = b[2:]
			if (Param.FOR == 1) :
				print(*b,sep=' ', file=ffor)
			if (Param.CPP == 1) :
				print(*b,sep='\t', file=fcpp)
			
		if (Param.FOR == 1) :
			ffor.close()
		if (Param.CPP == 1) :
			fcpp.close()

				
		
def main(argv=sys.argv):

	debut = time.time()

	Poros=Mesh(Param.FILE_2DM,Param.FILE_shpRB,Param.FILE_DEM)

	duree = time.time() - debut

	print("Execution Time : {0} minute(s)".format(round(duree/60)))



if __name__ == "__main__":
	sys.exit(main())







