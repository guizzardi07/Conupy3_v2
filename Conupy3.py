# -*- coding: utf-8 -*-

# Setup
'''
python -m venv c:\pythonEnviroments\conupy
C:\pythonEnviroments\conupy\Scripts\Activate.bat

pip install --upgrade pip
pip install -r env_connpy3.txt

#  Librerias instaladas:
    pip install wheel
    pip install pipwin

    pip install shapely
    pip install gdal
    pip install fiona
    pip install pyproj
    pip install six
    pip install rtree
    pip install geopandas
    pip install rasterio

    pip install geovoronoi

    pip freeze -l > env_connpy3.txt 
    pip install -r env_connpy3.txt
'''

import pandas as pd
import numpy as np
import pickle
import  os

import geopandas as gpd
import rasterio
from shapely.geometry import LineString,Polygon
from shapely.ops import unary_union

### Pra ignorar un warnings de shapely. 
import warnings
from shapely.errors import ShapelyDeprecationWarning
warnings.filterwarnings('ignore', category=ShapelyDeprecationWarning) 

from geovoronoi import voronoi_regions_from_coords, points_to_coords

from itertools import tee

from munch import Munch as Bunch
from collections import OrderedDict

# Funciones en el archivo Funciones/script_funciones.py
from  Funciones.script_funciones import (pairwise, dist, distToSegment, 
                                        insertPoints, removePoints, addNode, atraviesaArroyo,
                                        intersection)

# Clases en el archivo Funciones/clases_conupy.py
from Funciones.clases_conupy import NodesList
from Funciones.clases_conupy import LinksDict

def createGutters(nodos, links):
    # Crear sumideros
    print ('Creando sumideros')
    for (i, nodo) in enumerate(nodos):
        if nodo.type != 'corner':
            continue

        # Buscar el nodo conducto mas cercano
        nearINodes = nodos.getINodesNear(nodo.p, params['maxDistGutter'])
        mindist, minj = params['maxDistGutter'], -1
        for j in nearINodes:
            nodo2 = nodos[j]
            if nodo2.type == 'conduit':
                d = dist(nodo.p, nodo2.p)
                if d < mindist:
                    mindist = d
                    minj = j
        if minj == -1:
            continue
        # Existe un nodo conducto cerca (< 80m)
        n0, n1 = i, minj
        # Si ya existe una conexión entre los nodos
        if (n0, n1) in links or (n1, n0) in links:
            continue
        # Crear un sumidero
        links[(n0, n1)] = {'type':'gutter', 'w':params['xsSumideroW']}

def createWeirs(nodos, links):
    # Crear vertederos
    print ('Creando vertederos')
    for (i, nodo) in enumerate(nodos):
        if nodo.type != 'corner':
            continue

        # Buscar el tramo de arroyo mas cercano
        nearLinks = links.getLinksNear(nodo.p, params['maxDistWeir'])

        # Skip node if it already has a weir
        def alreadyHasWeir():
            for (n10, n11), link  in nearLinks.items():
                if n10 == i and link['type'] == 'weir':
                    return True
            return False

        if alreadyHasWeir():
            continue

        mindist, minj = params['maxDistWeir'], -1
        for (n10, n11), link  in nearLinks.items():
            if link['type'] != 'channel':
                continue
            p10, p11 = nodos[n10].p, nodos[n11].p
            d = distToSegment(nodo.p, p10, p11)
            if d < mindist:
                mindist = d
                minj = n10 if dist(nodo.p, p10) <= dist(nodo.p, p11) else n11

        if minj == -1:
            continue

        # Existe un nodo arroyo cerca, conectar
        n0, n1 = i, minj
        # Si ya existe una conexión entre los nodos
        if (n0, n1) in links or (n1, n0) in links:
            continue
        # Crear un vertedero
        links[(n0, n1)] = {'type':'weir', 'w':params['xsVertederoW']}

def calculateElevations(nodos, links, shpFileNodos, rasterFileDEM, rasterFileSlope, rasterFileImpermeabilidad, spatial_ref,rasterFileCurvaNumero, rasterTipoCurva, rasterFileSaSPerv,rasterFileNPerv):
    print ('Elevacion, pendiente e impermeabilidad')
    id_n = []
    coord_n = []
    for i, nodo in enumerate(nodos):
        id_n.append(i)
        coord_n.append(nodo.p)

    Data= {'Id':id_n,'CoorXY':coord_n,}
    Nodos_streets = pd.DataFrame(Data)
    Nodos_streets[['x','y']] = pd.DataFrame(Nodos_streets.CoorXY.tolist(), index= Nodos_streets.index)
    del Nodos_streets['CoorXY']

    # Arma el shp
    Nodos_streets_shp = gpd.GeoDataFrame(Nodos_streets, geometry=gpd.points_from_xy(Nodos_streets.x,Nodos_streets.y),crs=spatial_ref)

    # Agrega del Elev Terreno
    raster_dem = rasterio.open(rasterFileDEM)
    coords = [(x,y) for x, y in zip(Nodos_streets_shp.x, Nodos_streets_shp.y)]
    # Sample the raster at every point location and store values in DataFrame
    Nodos_streets_shp['elev'] = [x[0] for x in raster_dem.sample(coords)]

    # Agrega Slope
    raster_slope = rasterio.open(rasterFileSlope)
    Nodos_streets_shp['slope'] = [x[0] for x in raster_slope.sample(coords)]

    # Agrega Impermeabilidad
    raster_imperm = rasterio.open(rasterFileImpermeabilidad)
    Nodos_streets_shp['imper'] = [x[0] for x in raster_imperm.sample(coords)]

    # Agrega CN
    raster_CN = rasterio.open(rasterFileCurvaNumero)
    Nodos_streets_shp['CN'] = [x[0] for x in raster_CN.sample(coords)]

    # Agrega TipoCurva
    raster_TipoCurva = rasterio.open(rasterTipoCurva)
    Nodos_streets_shp['TipoCurva'] = [x[0] for x in raster_TipoCurva.sample(coords)]

    # Agrega SaSPerv
    raster_SaSPerv = rasterio.open(rasterFileSaSPerv)
    Nodos_streets_shp['SaSPerv'] = [x[0] for x in raster_SaSPerv.sample(coords)]

    # Agrega saNPerv
    raster_saNPerv = rasterio.open(rasterFileNPerv)
    Nodos_streets_shp['saNPerv'] = [x[0] for x in raster_saNPerv.sample(coords)]

    Nodos_streets_shp.to_file(shpFileNodos, driver='ESRI Shapefile')

    for i, nodo in enumerate(nodos):
        nodo.elev = float(Nodos_streets_shp.loc[Nodos_streets_shp.Id == i,'elev'].iloc[0])
        nodo.slope = float(Nodos_streets_shp.loc[Nodos_streets_shp.Id == i,'slope'].iloc[0])
        nodo.imper = float(Nodos_streets_shp.loc[Nodos_streets_shp.Id == i,'imper'].iloc[0])
        nodo.curvaNumero = float(Nodos_streets_shp.loc[Nodos_streets_shp.Id == i,'CN'].iloc[0])
        nodo.TipoCurva = float(Nodos_streets_shp.loc[Nodos_streets_shp.Id == i,'TipoCurva'].iloc[0])
        nodo.SaSPerv = float(Nodos_streets_shp.loc[Nodos_streets_shp.Id == i,'SaSPerv'].iloc[0])
        nodo.saNPerv = float(Nodos_streets_shp.loc[Nodos_streets_shp.Id == i,'saNPerv'].iloc[0])
        nodo.offset = 0
        nodo.length = 0
    
    for (n0, n1), link in links.items():
        # Acumulate the length of segments incoming to the node
        if link['type'] in ['street', 'channel', '2d', 'culvert']:
            length = dist(nodos[n0].p, nodos[n1].p)
            nodos[n0].length += length
            nodos[n1].length += length

        # Ensure conduits are always below terrain level
        if link['type'] == 'conduit':
            link['levelIni'] = min(link['levelIni'],
                nodos[n0].elev - link['h'] - params['minCoverage'])
            link['levelFin'] = min(link['levelFin'],
                nodos[n1].elev - link['h'] - params['minCoverage'])

        # Update invert offset if appropiate
        if link['type'] in ['conduit', 'channel', 'culvert']:
            offset0 = link['levelIni'] - nodos[n0].elev
            offset1 = link['levelFin'] - nodos[n1].elev
            nodos[n0].offset = min(nodos[n0].offset, offset0)
            nodos[n1].offset = min(nodos[n1].offset, offset1)

    # Set elevations for streets, gutters and weirs
    for (n0, n1), link in links.items():
        if link['type'] == 'street':
            link['levelIni'] = nodos[n0].elev
            link['levelFin'] = nodos[n1].elev
        elif link['type'] in ['weir', 'gutter']:
            link['levelIni'] = nodos[n0].elev
            link['levelFin'] = nodos[n1].elev + nodos[n1].offset

            if link['type'] in 'weir':
                if link['levelFin'] > link['levelIni']:
                    print ('WARNING: Weir %s - channel node invert is higher than the weir crest. This creates instabilties in SWMM.')

def createOutfallNodes(nodos, shpFileNodosBorde):
    outfalls = gpd.read_file(shpFileNodosBorde)
    outfalls['points'] = outfalls.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
    
    nodosOutfall = []
    lineasOutfall = []
    for (i, nodo) in enumerate(nodos):
        # Buscar la posicion outfall mas cercana
        mindist, minj = params['maxDistConnectOutfallNodes'], -1
        for j, outfall in outfalls.iterrows():
            if dist(nodo.p, outfall.points) < mindist:
                mindist = dist(nodo.p, outfall.points)
                minj = j
        if minj == -1:
            continue

        coord_outfall = np.array(outfalls.loc[minj,'points'][0])
        nodoOutfall = Bunch(p = coord_outfall,
                            elev = nodo.elev,
                            offset = nodo.offset)
        nodosOutfall.append(nodoOutfall)

        lineasOutfall.append( [i, len(nodosOutfall)-1, params['outfallXsWidth']] )

    return nodosOutfall, lineasOutfall

def calculateDeadDepths(nodos, links, lineasOutfall):
    for nodo in nodos:
        nodo.tirante = 0#100
    for linea in lineasOutfall:
        nodo0 = linea[0]
        nodos[nodo0].tirante = 0

    # Collect links data
    linksData = []
    for n0, n1 in links:
        link = links[(n0, n1)]
        linklevel = max(link['levelIni'], link['levelFin'])
        linksData.append((n0, n1, linklevel))
    # Optimization: Sort links from lower to higher levels
    linksData.sort(key=lambda tup: tup[2])

    i = 0
    maxbajada = 0.001#100
    while maxbajada > 0.01:
        maxbajada = 0

        def igualar(nodo0, nodo1, linklevel):
            elev0 = nodos[nodo0].elev + nodos[nodo0].offset + nodos[nodo0].tirante
            elev1 = nodos[nodo1].elev + nodos[nodo1].offset + nodos[nodo1].tirante

            newEle0 = min(elev0, max(linklevel, elev1))
            newEle1 = min(elev1, max(linklevel, elev0))

            bajada = max(elev0 - newEle0, elev1 - newEle1)

            nodos[nodo0].tirante = newEle0 - (nodos[nodo0].elev + nodos[nodo0].offset)
            nodos[nodo1].tirante = newEle1 - (nodos[nodo1].elev + nodos[nodo1].offset)

            return bajada

        for n0, n1, linklevel in linksData:
            bajada = igualar(n0, n1, linklevel)
            maxbajada = max(maxbajada, bajada)

        i += 1
        print ('Iteracion %i - Max Bajada %f' % (i, maxbajada))


    # Elevate nodes to eliminate dead depths
    print ('Raising nodes to eliminate dead depths')
    for n0, n1 in links:
        link = links[(n0, n1)]
        link['levelIni'] += nodos[n0].tirante
        link['levelFin'] += nodos[n1].tirante
    for nodo in nodos:
        nodo.elev += nodo.tirante

def writeNetworkShapes(nodos, links, shpFileNodos, shpFileLineas, spatial_ref):
    print (u'Número de nodos: %d' % len(nodos))
    print (u'Número de links: %d' % len(links))

    # Escribir shape con la posicion de los nodos
    campos = OrderedDict()
    campos['name'] = ['NODO%d' % i for i, _ in enumerate(nodos)]
    campos['type'] = [nodo.type for nodo in nodos]
    campos['elev'] = [float(nodo.elev) for nodo in nodos]
    campos['offs'] = [float(nodo.offset) for nodo in nodos]
    campos['inve'] = [float(nodo.elev + nodo.offset) for nodo in nodos]
    campos['dead'] = [float(nodo.tirante) for nodo in nodos]
    campos['CoorXY'] = [nodo.p for nodo in nodos]

    df_nodos = pd.DataFrame(campos)
    df_nodos[['x','y']] = pd.DataFrame(df_nodos.CoorXY.tolist(), index= df_nodos.index)
    del df_nodos['CoorXY']
    shp_nodos = gpd.GeoDataFrame(df_nodos, geometry=gpd.points_from_xy(df_nodos.x,df_nodos.y),crs=spatial_ref)
    del shp_nodos['x']
    del shp_nodos['y']
    shp_nodos.to_file(shpFileNodos, driver='ESRI Shapefile')

    # Escribir shape con los links
    #polilineas = []
    campos = OrderedDict()
    campos['name'] = []
    campos['n0'] = []
    campos['n1'] = []
    campos['type']  = []
    campos['w'] = []
    campos['geometry'] = []
    campos['elev0'] = []
    campos['elev1'] = []
    for i, (n0, n1) in enumerate(links):
        link = links[(n0, n1)]
        #polilineas.append([nodos[n0].p, nodos[n1].p])
        campos['name'].append(str(link['type']) + str(i))
        campos['n0'].append(int(n0))
        campos['n1'].append(int(n1))
        campos['type'].append(str(link['type']))
        campos['geometry'].append(LineString([nodos[n0].p, nodos[n1].p]))
        
        try:
            campos['w'].append(float(link.get('w', -1.0)))
            campos['elev0'].append(float(link['levelIni']))
            campos['elev1'].append(float(link['levelFin']))
        except:
            print ('ERROR: Link (%d,%d) with type %s is missing w, levelIni or levelFin data.' % (n0, n1, link['type']))
            print (link)
    
    df_links = pd.DataFrame(campos)
    shp_links = gpd.GeoDataFrame(df_links, geometry='geometry',crs=spatial_ref)
    shp_links.to_file(shpFileLineas, driver='ESRI Shapefile')

def createSubcatchments(nodos, shpFileCuenca, spatial_ref):
        # Crear centros de cuencas
        centros = []
        for (i, nodo) in enumerate(nodos):
            if nodo.type != 'conduit':
                centros.append([nodo.p[0], nodo.p[1], i])

        xmax = max(nodo.p[0] for nodo in nodos)
        xmin = min(nodo.p[0] for nodo in nodos)
        ymax = max(nodo.p[1] for nodo in nodos)
        ymin = min(nodo.p[1] for nodo in nodos)

        centros2 = centros[:]
        centros2.append([0.5 * xmin + 0.5 * xmax, 2.0 * ymax - 1.0 * ymin, -1])
        centros2.append([0.5 * xmin + 0.5 * xmax, 2.0 * ymin - 1.0 * ymax, -1])
        centros2.append([2.0 * xmax - 1.0 * xmin, 0.5 * ymin + 0.5 * ymax, -1])
        centros2.append([2.0 * xmin - 1.0 * xmax, 0.5 * ymin + 0.5 * ymax, -1])
        print ('Numero de cuencas:\t%d' % len(centros))
        
        # Escribir shape con la posicion de los baricentros de subcuencas
        Data= {'CoorXY':centros2,}
        df_centros = pd.DataFrame(Data)
        df_centros[['x','y','id']] = pd.DataFrame(df_centros.CoorXY.tolist(), index= df_centros.index)

        shp_centros = gpd.GeoDataFrame(df_centros, geometry=gpd.points_from_xy(df_centros.x,df_centros.y),crs=spatial_ref)
        del shp_centros['CoorXY']
        #del shp_centros['x']
        #del shp_centros['y']
        shp_centros.to_file(shpFileCentros, driver='ESRI Shapefile')
        
        # Escribir shape con la posicion de los baricentros de subcuencas
        create_thiessen_polygons(shpFileCentros, subcuencasShpFile)
        clip_feature(subcuencasShpFile, shpFileCuenca, subcuencasClipShpFile)

        subcuencas = gpd.read_file(subcuencasClipShpFile)
        subcuencas['area'] = subcuencas['geometry'].area
        
        subcuencasDict = {}
        for idx, subcuenca in subcuencas.iterrows():
            fid = subcuenca['FID']
            area = subcuenca['area']

            if subcuenca['geometry'].type == 'MultiPolygon': 
                lst_coord_polig = []
                for x in subcuenca['geometry'].geoms:
                    lst_coord_polig.extend(list(x.exterior.coords))
            else:
                lst_coord_polig = list(subcuenca['geometry'].exterior.coords)
            poligono = []
            
            for pt_xy in lst_coord_polig:
                poligono.append(np.array(pt_xy))
            subcuencasDict[fid] = [poligono, area]

        # Completar si falta alguna y eliminar duplicadas si existieran
        subcuencasCompletas = []
        for i in range(0,len(centros)):
            subcuencasCompletas.append(subcuencasDict.get(i, [[], 0]))
        
        return centros, subcuencasCompletas

def create_thiessen_polygons(nodesFile, fileName):
    print ('Creando polígonos de Thiessen. Destino: ' + fileName)

    shp_Centros = gpd.read_file(nodesFile)
    spatial_ref = shp_Centros.crs

    area_max_lon = shp_Centros.x.max()
    area_min_lon = shp_Centros.x.min()
    area_max_lat = shp_Centros.y.max()
    area_min_lat = shp_Centros.y.min()

    lat_point_list = [area_min_lat, area_max_lat,area_max_lat,area_min_lat]
    lon_point_list = [area_min_lon, area_min_lon, area_max_lon, area_max_lon]
    polygon_geom = Polygon(zip(lon_point_list, lat_point_list))

    Data= {'Id':[1,],'geometry':polygon_geom}
    df_boundary = pd.DataFrame(Data)
    boundary = gpd.GeoDataFrame(df_boundary, geometry='geometry',crs=spatial_ref)

    boundary_shape = unary_union(boundary.geometry)
    coords = points_to_coords(shp_Centros.geometry)

    # Calculate Voronoi Regions
    region_polys, region_pts = voronoi_regions_from_coords(coords, boundary_shape)
    ''' region_polys -> id_region : shapely.polygon
        region_pts -> id_region : id_coords (ej: 'P0')  '''

    thiessen_polygons = gpd.GeoDataFrame(gpd.GeoSeries(region_polys), columns=['geometry'])
    thiessen_polygons['id'] = thiessen_polygons.index.map(region_pts)
    thiessen_polygons['id'] = thiessen_polygons.id.apply(lambda x: x[0])
    thiessen_polygons = thiessen_polygons.rename(columns={'id':'FID',})
    thiessen_polygons.to_file(fileName, crs=spatial_ref,driver='ESRI Shapefile')
    print ('Finalizada creación de polígonos.')

def clip_feature(originalFile, clipPolygonFile, fileName):
    print ('Recortando feature. Destino: ' + fileName)

    originalLayer = gpd.read_file(originalFile)
    clipLayer = gpd.read_file(clipPolygonFile)
    subcuencas_clip = gpd.clip(originalLayer, clipLayer)
    subcuencas_clip.to_file(fileName, driver='ESRI Shapefile')

    print ('Finalizado el recortado.')

def createRainGagesMethod0(centros, gageFileName, rasterFileCoeficiente, gagesFileName):
    print ('Leyendo pluviometro maestro...')
    print(gageFileName)

    with open(gageFileName, 'r') as iFile:
        lineas = iFile.readlines()


    print ('Creando estaciones...')
    gages = []
    with open(gagesFileName, 'w') as tF:
        for i in range(0, params['numDiscreteGages']):
            coef = i * (100/(params['numDiscreteGages']-1))

            gage = {}
            gage['name'] = 'GAGE'+str(coef)
            _, gage['file'] = os.path.split(gagesFileName)
            gage['interval'] = '0:10'
            gages.append(gage)

            for linea in lineas:
                datos = linea.split()
                if len(datos) == 0:
                    continue
                datos[0]  = gage['name']
                datos[-1] = float(datos[-1])*(float(coef)/100.0)
                tF.write(('').join([ str(x).ljust(15, ' ') for x in datos]))
                tF.write('\n')

    print ('Leyendo mapa de decaimiento')
    nodosDecaimiento = sample_raster_on_nodes(shpFileNodos, rasterFileCoeficiente)

    print ('Seleccionando pluviómetro para cada subcuenca')
    subcatchmentGages = []
    for (i, centro) in enumerate(centros):
        coef = int(nodosDecaimiento[centro[2]])
        gageName = 'GAGE' + str(coef - (coef%(100/(params['numDiscreteGages']-1))))
        subcatchmentGages.append(gageName)

    return gages, subcatchmentGages

def createRainGagesMethod1(centros, stationsFileName):
    print ('Leyendo lista de pluviometros...')
    gages = []
    with open(stationsFileName, 'r') as f:
        for i, line in enumerate(f):
            data = line.split()
            gage = Bunch()
            gage.coord = np.array([float(data[0]), float(data[1])])
            gage.name = data[2]
            gage.file = data[3]
            gage.interval = '0:10'
            gages.append(gage)

    print ('Seleccionando pluviómetro para cada subcuenca...')
    subcatchmentGages = []
    for (i, centro) in enumerate(centros):
        minDist, minGage = 1e10, None

        for gage in gages:
            gDist = dist(gage.coord, np.array(centro[0:2]))

            if gDist < minDist:
                minGage, minDist = gage, gDist

        subcatchmentGages.append(minGage.name)

    return gages, subcatchmentGages

def saveOnFile( workspace,data, name ):
    oF = open(workspace+name + '.pkl', 'wb')
    pickle.dump(data, oF)
    oF.close()

def readDrainageNetwork(nodos, links, shpFileDrainagePrepared):
        # Read the prepared drainage network
        streams = gpd.read_file(shpFileDrainagePrepared)
        streams = streams.rename(columns={'Alto':'h','Ancho':'w','Tipo':'typ','Transect':'transect'})
        streams['levelIni'] = pd.to_numeric(streams['levelIni'])
        streams['levelFin'] = pd.to_numeric(streams['levelFin'])

        streams['points'] = streams.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
        
        # Subdivide the streams in spans shorter than maxLengthForStreamSpanDivide
        for idx, stream in streams.iterrows():
            lst_points = stream['points']
            length  = insertPoints(lst_points, params['maxLengthForStreamSpanDivide'])
            streams.loc[idx,'points'] = [lst_points,] 
        
        snappedPoints = []
        # Snap the end nodes of each stream to nearby nodes of other streams
        for idx, stream in streams.iterrows():

            def snap(p, idx, streams):
                for idx2, targetStream in streams.iterrows():
                    if (idx == idx2):
                        continue
                    mindist, minj = params['maxDistSnapStreamNodes'], -1
                    lst_pts_targetStream = targetStream.points[0]
                    for (j, p2) in enumerate(lst_pts_targetStream):
                        d = dist(p, p2)
                        if d < mindist:
                            mindist, minj = d, j    ###### Porque cambia mindist? ######################################################################################################################################
                    if minj == -1:
                        continue
                    # Convert coordinates to tuple because numpy arrays don't
                    # work properly with the 'in' operand
                    snappedPoints.append(tuple(lst_pts_targetStream[minj]))
                    return lst_pts_targetStream[minj]
                return p

            lst_pts_stream = stream.points[0]
            lst_pts_stream[0] = snap(lst_pts_stream[0],  idx, streams)
            lst_pts_stream[-1] = snap(lst_pts_stream[-1], idx, streams)
            streams.loc[idx,'points'] = [lst_pts_stream,]
        
        # Join streams spans shorter than minLengthForStreamSpanJoin unless they've been 'snapped'
        for idx, stream in streams.iterrows():
            lst_pts_Stream = stream.points[0]
            lst_pts_Stream_update = removePoints(lst_pts_Stream, params['minLengthForStreamSpanJoin'], snappedPoints)
            streams.loc[idx,'points'] = [lst_pts_Stream_update,] 
        
        for idx, stream in streams.iterrows():
            tipoTramo = 'conduit' if stream.typ in ['entubado', 'conducto', 'conduit'] else 'channel'

            lst_pts_Stream = stream.points[0]

            nodesn = [addNode(nodos, p, tipoTramo) for p in lst_pts_Stream]

            length = sum(dist(nodos[n0].p, nodos[n1].p) for n0, n1 in pairwise(nodesn))

            progFin = 0
            for n0, n1 in pairwise(nodesn):
                progIni = progFin
                progFin = progIni + dist(nodos[n0].p, nodos[n1].p)
                if n0 == n1:
                    continue
                
                # Create a new link
                links[(n0, n1)] = {'type': tipoTramo,
                                   'w': stream.w,
                                   'h': stream.h,
                                   'levelIni': stream.levelIni + (stream.levelFin - stream.levelIni) * progIni / length,
                                   'levelFin': stream.levelIni + (stream.levelFin - stream.levelIni) * progFin / length,
                                   'transect': stream.transect}

        print ('\tNumber of nodes for the drainage network: %i' % len(nodos))
        print ('\tNumber of links for the drainage network: %i' % len(links))
        print ('FINISHED: Drainage network construction')

def readStreets(nodos, links, shpFileCalles):
    streets = gpd.read_file(shpFileCalles)
    streets = streets.rename(columns={'ANCHO':'w',})
    streets = streets.rename(columns={'Ancho':'w',})
    streets = streets.rename(columns={'ancho':'w',})
    streets = streets.rename(columns={'manning':'n_manning',})
    
    streets['w'] = pd.to_numeric(streets['w'])
    streets['points'] = streets.apply(lambda x: [y for y in x['geometry'].coords], axis=1)    

    for idx, street in streets.iterrows():
        lst_street_points = street['points']
        if len(lst_street_points) < 2:
            continue
        
        # Subdivide the street in spans shorter than maxLengthForStreetSpanDivide
        length  = insertPoints(lst_street_points, params['maxLengthForStreetSpanDivide'])
        #street.loc[idx,'points'] = [lst_street_points,]

        # Join streams spans shorter than minLengthForStreetSpanJoin
        lst_street_points_update = removePoints(lst_street_points, params['minLengthForStreetSpanJoin'], [])
        streets.loc[idx,'points'] = [lst_street_points_update,]
    
    for idx, street in streets.iterrows():
        if idx % 100 == 0:
            print ('Procesando calle ' + str(idx) + ' de ' + str(len(streets)))
        
        lst_street_points = street['points'][0]

        # Si el ancho es nulo
        if (street.w == 0):
            continue

        if len(lst_street_points) < 2:
            continue

        for p0, p1 in pairwise(lst_street_points):
            # Calculate
            # Crear u obtener los nodos extremos
            n0 = addNode(nodos, p0, 'corner')
            n1 = addNode(nodos, p1, 'corner')

            # Si la polilinea es demasiado corta
            if n0 == n1:
                continue

            # Verificar si la calle atraviesa un arroyo
            channel_data = atraviesaArroyo(nodos[n0].p, nodos[n1].p, nodos, links)

            if channel_data is not None:
                nch0, nch1, channel_link = channel_data
                p4 = intersection(nodos[n0].p, nodos[n1].p, nodos[nch0].p, nodos[nch1].p)

                nchannel = nch0 if dist(nodos[nch0].p, p4) <= dist(nodos[nch1].p, p4) else nch1

                # Dividir primero el arroyo si vale la pena (distancia al nodo más
                # cercano > a 15% maxLengthForStreamSpanDivide)
                if (min(dist(nodos[nch0].p, p4), dist(nodos[nch1].p, p4)) >
                    params['maxLengthForStreamSpanDivide'] * 0.15):

                    nchannel = addNode(nodos, p4, 'channel', 0)

                    alpha = dist(nodos[nch0].p, nodos[nchannel].p) / dist(nodos[nch0].p, nodos[nch1].p)
                    levelMid = (1 - alpha) * channel_link['levelIni'] + alpha * channel_link['levelFin']

                    links[(nch0, nchannel)] = {'type': channel_link['type'],
                                                'w': channel_link['w'],
                                                'h': channel_link['h'],
                                                'levelIni': channel_link['levelIni'],
                                                'levelFin': levelMid,
                                                'transect': channel_link['transect']}
                    links[(nchannel, nch1)] = {'type': channel_link['type'],
                                                'w': channel_link['w'],
                                                'h': channel_link['h'],

                                                'levelIni': levelMid,
                                                'levelFin': channel_link['levelFin'],
                                                'transect': channel_link['transect']}
                    del links[(nch0, nch1)]

                # Atraviesa un arroyo --> crear conexión entre el las dos esquinas y el arroyo
                if n0 != nchannel:
                    # Crear un nuevo vertedero
                    links[(n0, nchannel)] = {'type':'weir', 'w':street.w,'Nmann':street.n_manning}
                if n1 != nchannel:
                    # Crear un nuevo vertedero
                    links[(n1, nchannel)] = {'type':'weir', 'w':street.w,'Nmann':street.n_manning}
            else:
                # No atraviesa --> crear una calle comun
                links[(n0, n1)] = {'type':'street', 'w':street.w,'Nmann':street.n_manning,'Id_calle':street.Id_calle,'transect':street.Id_calle}

def writeSWMMFile(nodos, links, centros, subcuencas, nodosOutfall, lineasOutfall, gages, subcatchmentGages, swmmInputFileName, modelFolder):
    for nodo in nodos:
        nodo.area = 1.167
    for (i, centro) in enumerate(centros):
        nodos[centro[2]].area = params['juApondPer'] * subcuencas[i][1]

    print('write SWMMFile')

    with open(swmmInputFileName, 'w') as tF:
        tF.write('[TITLE]\n')
        tF.write(';;Project Title/Notes\n')
        tF.write('Conurbano\n')
        tF.write('\n')
        tF.write('[OPTIONS]\n')
        tF.write(';;Option             Value\n')
        tF.write('FLOW_UNITS           CMS\n')
        tF.write('INFILTRATION         CURVE_NUMBER\n') #HORTON
        tF.write('FLOW_ROUTING         DYNWAVE\n')
        tF.write('LINK_OFFSETS         ELEVATION\n')
        tF.write('MIN_SLOPE            0.0001\n')
        tF.write('ALLOW_PONDING        YES\n')
        tF.write('SKIP_STEADY_STATE    NO\n')
        tF.write('\n')
        tF.write('START_DATE           03/10/2023\n')
        tF.write('START_TIME           18:00\n')
        tF.write('REPORT_START_DATE    03/10/2023\n')
        tF.write('REPORT_START_TIME    18:00:00\n')
        tF.write('END_DATE             03/11/2023\n')
        tF.write('END_TIME             06:00\n')
        tF.write('SWEEP_START          01/01\n')
        tF.write('SWEEP_END            12/31\n')
        tF.write('DRY_DAYS             0\n')
        tF.write('REPORT_STEP          00:15:00\n')
        tF.write('WET_STEP             00:00:30\n')
        tF.write('DRY_STEP             00:00:30\n')
        tF.write('ROUTING_STEP         00:00:01\n')
        tF.write('\n')
        tF.write('INERTIAL_DAMPING     FULL\n')
        tF.write('NORMAL_FLOW_LIMITED  BOTH\n')
        tF.write('FORCE_MAIN_EQUATION  H-W\n')
        tF.write('VARIABLE_STEP        0.75\n')
        tF.write('LENGTHENING_STEP     5\n')
        tF.write('MIN_SURFAREA         1.167\n')
        tF.write('MAX_TRIALS           20\n')
        tF.write('HEAD_TOLERANCE       0.002\n')
        tF.write('SYS_FLOW_TOL         5\n')
        tF.write('LAT_FLOW_TOL         5\n')
        tF.write('MINIMUM_STEP         0.5\n')
        tF.write('THREADS              4\n')
        tF.write('\n')

        tF.write('[FILES]\n')
        tF.write(';;Interfacing Files\n')
        tF.write('SAVE RAINFALL        rainfall.rff\n')
        tF.write(';;SAVE RUNOFF          runoff.rof\n')
        tF.write('SAVE OUTFLOWS        outflows.txt\n')
        tF.write('\n')

        tF.write('[EVAPORATION]\n')
        tF.write(';;Data Source    Parameters\n')
        tF.write(';;-------------- ----------------\n')
        tF.write('CONSTANT         %.3f\n' % params['evap'])
        tF.write('DRY_ONLY         NO\n')
        tF.write('\n')

        tF.write('[RAINGAGES]\n')
        tF.write(';;Name         Format         Interval       SCF            DataSrc        SourceName     Sta            Units\n')
        tF.write(';;============================================================================================================\n')
        for gage in gages:
            list = [gage['name'], 'INTENSITY', gage['interval'], 1.0, 'FILE', gage['file'], gage['name'], 'MM']
            tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            tF.write('\n')

        tF.write('\n')
        tF.write('[SUBCATCHMENTS]\n')
        tF.write(';;Name         Raingage       Outlet         Area           %ImperV        Width          Slope          Curve Length\n')
        tF.write(';;======================================================================================================================\n')
        for (i, centro) in enumerate(centros):
            n1 = centro[2]
            nodo = nodos[n1]
            gageName = subcatchmentGages[i]
            list = ['CUENCA'+str(i), gageName, 'NODO'+str(n1), '%.3f' % (float(subcuencas[i][1])/10000.0), '%.3f' % nodo.imper, '%.3f' % (nodo.length/2), '%.3f' % (nodo.slope), '%.3f' % (subcuencas[i][1]**0.5)]
            tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            tF.write('\n')


        """ 
        tF.write('\n')
        tF.write('[SUBAREAS]\n')
        tF.write(';;Name         N_Imp          N_Perv         S_Imp          S_Perv         %ZER           RouteTo\n')
        tF.write(';;=======================================================================================================\n')
        for (i, centro) in enumerate(centros):
            list = ['CUENCA'+str(i), params['saNImp'], params['saNPerv'], params['saSImp'], params['saSPerv'], params['saZero'], 'OUTLET']
            tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            tF.write('\n') 
        """

        tF.write("\n")
        tF.write("[SUBAREAS]\n")
        tF.write(";;Name         N_Imp          N_Perv         S_Imp          S_Perv         %ZER           RouteTo\n")
        tF.write(";;=======================================================================================================\n")
        for (i, centro) in enumerate(centros):
            n1 = centro[2]
            nodo = nodos[n1]
            #params["saNPerv"] =  (1-nodo.UrbanoRural)*0.025+0.025 #Para la zona rural donde UrbanoRural=1, n=0.025. En cambio para la zona rural donde UrbanoRural=0, n=0.05
            list = ['CUENCA'+str(i), params["saNImp"], "%.3f" % nodo.saNPerv, params["saSImp"],  "%.3f" % nodo.SaSPerv, params["saZero"], 'OUTLET']

            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")

        ## Inf por HORTON
        # tF.write('\n')
        # tF.write('[INFILTRATION]\n')
        # tF.write(';;Subcat       MaxRate        MinRate        Decay          DryTime        Max Inf\n')
        # tF.write(';;========================================================================================\n')
        # for (i, centro) in enumerate(centros):
            # list = ['CUENCA'+str(i), params['inF0'], params['inFf'], params['inCoefDecaim'], params['inDryTime'], params['inMaxInf']]
            # tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            # tF.write('\n')

        ## Inf por CN
        tF.write("\n")
        tF.write("[INFILTRATION]\n")
        tF.write(";;Subcatchment   CurveNum              DryTime\n")   #cuando se usa CN
		#tF.write(";;Subcat       MaxRate        MinRate        Decay          DryTime        Max Inf\n")  #cuando se usa HORTON
        tF.write(";;========================================================================================\n")
        for (i, centro) in enumerate(centros):
            n1 = centro[2]
            nodo = nodos[n1]
            list = ['CUENCA'+str(i), "%.3f" % nodo.curvaNumero, params["Conductividad"], params["inDryTime"]]
            tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
            tF.write("\n")

        # Depth from ground to invert elevation (ft or m) (default is 0).
        # juYmax = 1.0
        # maximum additional head above ground elevation that manhole junction can sustain under surcharge conditions (ft or m) (default is 0).
        # juYsur = 0
        # % del area de la subcuenca en que se almacena el agua en el nodo
        params['juApondPer'] = 0.75
        tF.write('\n')
        tF.write('[JUNCTIONS]\n')
        tF.write(';;Name         Elev           Ymax           Y0             Ysur           Apond) \n')
        tF.write(';;========================================================================================\n')
        for (i, nodo) in enumerate(nodos):
            if not nodo.type in params['nodeTypesAsJunctions']:
                continue
            list = ['NODO%d' % i,
                    '%.3f' % (nodo.elev + nodo.offset),
                    '%.3f' % (-nodo.offset + 20),
                    '%.3f' % params['juY0'],
                    '%.3f' % 0,
                    '%.3f' % nodo.area]
            tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            tF.write('\n')


        tF.write('\n')
        tF.write('[STORAGE]\n')
        tF.write(';;Name         Elev           Ymax           Y0             TABULAR        Apond          ) \n')
        tF.write(';;========================================================================================\n')
        for (i, nodo) in enumerate(nodos):
            if nodo.type in params['nodeTypesAsJunctions']:
                continue
            tF.write(('').join([ str(x).ljust(15, ' ') for x in [
                'NODO%d' % i,
                '%.3f' % (nodo.elev + nodo.offset),
                '%.3f' % (-nodo.offset+20),
                '%.3f' % params['juY0'],
                'TABULAR',
                'STORAGE%d' % i,
                '%.3f' % nodo.area]]))
            tF.write('\n')

        tF.write("\n")

        tF.write("[CURVES]\n")                                                                                      #En esta sección se cargan las curvas de almacenamiento según uso de suelo.
        tF.write(";;Name         Type           x-value        y-value\n")
        tF.write(";;========================================================================================\n")
        for (i, nodo) in enumerate(nodos):
            if (nodo.TipoCurva == 1):                                                                             #El raster que tiene la zonificación de uso de suelo es UrbanoRural
                if (nodo.offset < 0):
                    list = [
                        'STORAGE%d' % i,
                        'STORAGE',
                        "%.3f" % 0,                   "%.3f" % 0,
                        "%.3f" % (-nodo.offset + 1.5),  "%.3f" % (nodo.area*1/100),
                        "%.3f" % (-nodo.offset + 3),  "%.3f" % (nodo.area*10/100),
                        "%.3f" % (-nodo.offset + 9.5),  "%.3f" % (nodo.area*90/100),
                        "%.3f" % (-nodo.offset + 16), "%.3f" % nodo.area]
                else:
                    list = [
                        'STORAGE%d' % i,
                        'STORAGE',
                        "%.3f" % 0,                   "%.3f" % 0,
                        "%.3f" % 1.5,                   "%.3f" % (nodo.area*1/100),
                        "%.3f" % 3,                   "%.3f" % (nodo.area*10/100),
                        "%.3f" % 9.5,                   "%.3f" % (nodo.area*90/100),
                        "%.3f" % 16,                  "%.3f" % nodo.area]
                tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                tF.write("\n")
            elif(nodo.TipoCurva == 2):

                if (nodo.offset < 0):
                    list = [
                        'STORAGE%d' % i,
                        'STORAGE',
                        "%.3f" % 0,                   "%.3f" % 0,
                        "%.3f" % (-nodo.offset + 1),                   "%.3f" % (nodo.area*1/100),
                        "%.3f" % (-nodo.offset + 3.8),                   "%.3f" % (nodo.area*30/100),
                        "%.3f" % (-nodo.offset + 7.5),                   "%.3f" % (nodo.area*93/100),
                        "%.3f" % (-nodo.offset + 10.5),                  "%.3f" % nodo.area]
                else:
                    list = [
                        'STORAGE%d' % i,
                        'STORAGE',
                        "%.3f" % 0,                   "%.3f" % 0,
                        "%.3f" % 1,                   "%.3f" % (nodo.area*1/100),
                        "%.3f" % 3.8,                   "%.3f" % (nodo.area*30/100),
                        "%.3f" % 7.5,                   "%.3f" % (nodo.area*93/100),
                        "%.3f" % 10.5,                  "%.3f" % nodo.area]
                tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                tF.write("\n")
                
            else:#elif(nodo.TipoCurva == 3):
                if (nodo.offset < 0):
                    list = [
                        'STORAGE%d' % i,
                        'STORAGE',
                        "%.3f" % 0,                   "%.3f" % 0,
                        "%.3f" % (-nodo.offset+0.5),      "%.3f" % (nodo.area),
                        "%.3f" % (-nodo.offset + 20), "%.3f" % (nodo.area)]
                else:
                    list = [
                        'STORAGE%d' % i,
                        'STORAGE',
                        "%.3f" % 0,                   "%.3f" % 0,
                        "%.3f" % 0.5,                 "%.3f" % (nodo.area),
                        "%.3f" % 20,                  "%.3f" % (nodo.area)]
                    
                tF.write(("").join([ str(x).ljust(15, ' ') for x in list]))
                tF.write("\n")

        tF.write('\n')

        tF.write('[OUTFALLS]\n')
        tF.write(';;Name         Elev           Type           Gate\n')
        tF.write(';;==========================================================\n')
        for (i, nodo) in enumerate(nodosOutfall):
            list = [
                'NODOOUT'+str(i),
                round(nodo.elev + nodo.offset,2),
                'FREE',
                'NO']
            tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            tF.write('\n')

        tF.write('\n')

        tF.write('[CONDUITS]\n')
        tF.write(';;Name         Node1          Node2          Length         N              Z1             Z2             Q0\n')
        tF.write(';;======================================================================================================================\n')
        for i, ((in0, in1), link) in enumerate(links.items()):
            name = link['type'] + str(i)
            length = dist(nodos[in0].p, nodos[in1].p)
            if link['type'] == 'street':
                list = [
                    name,
                    'NODO'+str(in0),
                    'NODO'+str(in1),
                    '%.3f' % length,
                    '%.3f' % params['coN'],
                    '%.3f' % nodos[in0].elev,
                    '%.3f' % nodos[in1].elev,
                    0]
            elif link['type'] == '2d':
                list = [
                    name,
                    'NODO'+str(in0),
                    'NODO'+str(in1),
                    '%.3f' % length,
                    '%.3f' % params['coN'],
                    '%.3f' % nodos[in0].elev,
                    '%.3f' % nodos[in1].elev,
                    0]
            elif link['type'] in ['channel', 'conduit', 'culvert']:
                list = [
                    name,
                    'NODO'+str(in0),
                    'NODO'+str(in1),
                    '%.3f' % length,
                    '%.3f' % params['coN'],
                    '%.3f' % link['levelIni'],  #nodos[in0].elev,#
                    '%.3f' % link['levelFin'],  #nodos[in1].elev,#
                    0]
            else:
                continue
            tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            tF.write('\n')
        for (i, linea) in enumerate(lineasOutfall):
            in0, in1, ancho = linea
            name = 'SALIDA' + str(i)
            length = dist(nodos[in0].p, nodosOutfall[in1].p)
            list = [
                name,
                'NODO'+str(in0),
                'NODOOUT'+str(in1),
                '%.3f' % length,
                '%.3f' % params['coN'],
                '%.3f' % (nodos[in0].elev + nodos[in0].offset),
                '%.3f' % (nodos[in0].elev + nodos[in0].offset),
                0]
            tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            tF.write('\n')


        tF.write('\n')
        tF.write('[WEIRS]\n')
        tF.write(';;Name           From Node        To Node          Type         CrestHt    Qcoeff     Gated    EndCon   EndCoeff   Surcharge \n')
        tF.write(';;-------------- ---------------- ---------------- ------------ ---------- ---------- -------- -------- ---------- ----------\n')
        for i, ((in0, in1), link) in enumerate(links.items()):
            if link['type'] == 'weir':
                name = 'weir' + str(i)
                list = [
                    name,
                    'NODO'+str(in0),
                    'NODO'+str(in1),
                    'SIDEFLOW',
                    '%.3f' % (nodos[in0].elev + params['weAlturaCordon']),
                    params['weCd'],
                    'NO',
                    0,
                    0,
                    'NO']
            elif link['type'] == 'culvert':
                name = 'overpass' + str(i)
                list = [
                    name,
                    'NODO'+str(in0),
                    'NODO'+str(in1),
                    'ROADWAY',
                    '%.3f' % (link['levelIni'] + link['htop']),
                    params['culvertRoadwayCd'],
                    'NO',
                    0,
                    0,
                    'NO',
                    0,
                    0,
                    '%.3f' % link['l'],
                    'PAVED'
                    ]
            else:
                continue
            tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            tF.write('\n')


        tF.write('\n')
        tF.write('[ORIFICES]\n')
        tF.write(';;Name         Node1          Node2          Type           Offset         Cd             Flap           Orate\n')
        tF.write(';;======================================================================================================================\n')
        for i, ((in0, in1), link) in enumerate(links.items()):
            if link['type'] != 'gutter':
                continue
            name = 'gutter' + str(i)
            list = [
                name,
                'NODO'+str(in0),
                'NODO'+str(in1),
                'SIDE',
                '%.3f' % nodos[in0].elev,
                params['weCd'],
                'NO',
                0]
            tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            tF.write('\n')


        tF.write('\n')
        tF.write('[XSECTIONS]\n')
        tF.write(';;Link         Type           G1             G2             G3             G4\n')
        tF.write(';;========================================================================================\n')
        transectas = {}
        customTransects = set()
        
        for i, ((in0, in1), link) in enumerate(links.items()):
            name = link['type'] + str(i)

            if link['type'] == 'street':
                tname = link['type'] + str(int(link['w']))+str(int(link['Nmann']*1000))
                transectas[tname] = [link['type'], tname, link['w'], link.get('h',0),link['Nmann']]
                #tname = link['Id_calle']
                list = [name, 'IRREGULAR', tname, 0, 0, 0]
                customTransects.add(tname)
                #list = [name, 'RECT_OPEN', link.get('h',20), link['w'], 0, 0]
            elif link['type'] == '2d':
                list = [name, 'RECT_OPEN', 100, 20.0, 2, 0, 1]
            elif link['type'] == 'channel':
                #print(link)
                if link['transect'] != '':
                    tname = link['transect']
                    customTransects.add(tname)
                else:
                    tname = link['type'] + str(int(link['w'])) + 'x' + str(int(link['h']))
                    transectas[tname] = [link['type'], tname, link['w'], link.get('h',0)]
                list = [name, 'IRREGULAR', tname, 0, 0, 0]
            elif link['type'] == 'conduit':
                if link['transect'] == '':
                    link['transect'] = 'RECT_CLOSED'
                list = [name, link['transect'], link['h'], link['w'], 0, 0]
            elif link['type'] == 'culvert':
                # Write cross section for overpass weir
                list = ['overpass' + str(i), 'RECT_OPEN', params['xsVertederoH'], link['wtop'], 0, 0]
                tF.write(('').join([ str(x).ljust(14, ' ') + ' ' for x in list]))
                # Define cross section for underpass conduit
                list = [name, 'RECT_CLOSED', link['h'], link['w'], 0, 0, 1, params['culvertCode']]
                tF.write('\n')
            elif link['type'] == 'weir':
                list = [link['type']+str(i), 'RECT_OPEN', params['xsVertederoH'], params['xsVertederoW'], 0, 0]
            elif link['type'] == 'gutter':
                list = [link['type']+str(i), 'RECT_CLOSED', params['xsSumideroH'], params['xsSumideroW'], 0, 0]
            else:
                continue
            tF.write(('').join([ str(x).ljust(14, ' ') + ' ' for x in list]))
            tF.write('\n')
        for (i, linea) in enumerate(lineasOutfall):
            in0, in1, ancho = linea
            list = ['SALIDA'+str(i), 'RECT_OPEN', params['xsG1'], ancho, params['xsG3'], params['xsG4']]
            tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            tF.write('\n')


        tF.write('\n')
        tF.write('[TRANSECTS]\n')
        tF.write(';;NC           Nleft          Nright         Nchanl         x              x \n')
        tF.write(';;X1           Name           Nsta           Xleft          Xright         0 0 0 0 0 0 \n')
        tF.write(';;GR           Elev           GR             Elev           GR             Elev....\n')
        tF.write(';;========================================================================================\n')
        for key, value in transectas.items():
            tipo, tname, ancho, alto, manning = value
            if (tipo == 'channel'):
                traTiranteArroyo = alto
                list = ['NC', params['traNArroyoPlanicie'], params['traNArroyoPlanicie'], params['traNArroyoCauce']]
                tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
                tF.write('\n')
                list = ['X1', tname, 8, -ancho * 0.5, ancho * 0.5, 0, 0, 0, 0, 0, 0]
                tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
                tF.write('\n')
                list = ['GR',
                        traTiranteArroyo + 3.0,        -ancho*0.5 - params['traAnchoMargenArroyo'],
                        traTiranteArroyo + 0.5,        -ancho*0.5 - params['traAnchoMargenArroyo'],
                        traTiranteArroyo,              -ancho*0.5,
                        0,                             -ancho*0.5 + 0.25,
                        0,                              ancho*0.5 - 0.25,
                        traTiranteArroyo,               ancho*0.5,
                        traTiranteArroyo + 0.5,         ancho*0.5 + params['traAnchoMargenArroyo'],
                        traTiranteArroyo + 3.0,         ancho*0.5 + params['traAnchoMargenArroyo']]
                tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
                tF.write('\n')
                tF.write(';;-------------------------------------------\n')
            elif (tipo == 'conduit'):
                list = ['NC', params['traNConducto'], params['traNConducto'], params['traNConducto']]
                tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
                tF.write('\n')
                list = ['X1', tname, 6, -ancho * 0.5, ancho * 0.5, 0, 0, 0, 0, 0, 0]
                tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
                tF.write('\n')
                list = ['GR', -ancho*0.5, alto, -ancho*0.5, 0, ancho*0.5, 0, ancho*0.5, alto]
                tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
                tF.write('\n')
                tF.write(';;-------------------------------------------\n')
            elif (tipo == 'street'):
                #list = ['NC', params['traNCalle'], params['traNCalle'], params['traNCalle']]                                                        #Indica los mannning de vereda-calle-vereda
                list = ['NC', manning, manning, manning] 
                tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
                tF.write('\n')
                list = ['X1', tname, 9, -ancho * 0.5 + params['traAnchoVereda'], ancho * 0.5 - params['traAnchoVereda'], 0, 0, 0, 0, 0, 0]          #Define los bankstations. 9 son la cantidad de puntos que componen la secc.transversal
                tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
                tF.write('\n')
                list = ['GR',
                        params['traAltoCordon'] + 3.0, -ancho*0.5,                                                                                  #Define la secc. transversal como y,x
                        params['traAltoCordon'] + 0.2, -ancho*0.5,
                        params['traAltoCordon'],       -ancho*0.5 + params['traAnchoVereda'],
                        0,                             -ancho*0.5 + params['traAnchoVereda'],
                        params['traAltoCentroCalle'],   0,
                        0,                              ancho*0.5 - params['traAnchoVereda'],
                        params['traAltoCordon'],        ancho*0.5 - params['traAnchoVereda'],
                        params['traAltoCordon'] + 0.2,  ancho*0.5,
                        params['traAltoCordon'] + 3.0,  ancho*0.5]
                tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
                tF.write('\n')
                tF.write(';;-------------------------------------------\n')
        print(customTransects)
        for transect in customTransects:
            try:
                with open(modelFolder + '/archivos_dat/' + transect + '.dat', 'r') as iF:
                    for line in iF:
                        tF.write(line)
                #print ('Transect file "%s" loaded' % (modelFolder + '/' + transect + '.dat'))
            except:
                print ('WARNING: Transect file "%s" not found' % (modelFolder + '/' + transect + '.dat'))

            tF.write(';;-------------------------------------------\n')

        tF.write('\n')
        tF.write('[LOSSES]\n')
        tF.write(';;Link           Kentry     Kexit      Kavg       Flap Gate  Seepage   \n')
        tF.write(';;=========================================================================\n')
        for i, ((in0, in1), link) in enumerate(links.items()):
            name = link['type'] + str(i)
            if link['type'] == 'culvert':
                list = [name, params['culvertKentry'], params['culvertKexit'], 0, 'NO', 0]
            else:
                continue
            tF.write(('').join([ str(x).ljust(14, ' ') + ' ' for x in list]))
            tF.write('\n')

        tF.write('\n')
        tF.write('[REPORT]\n')
        tF.write('CONTINUITY               YES\n')
        tF.write('FLOWSTATS                YES\n')
        tF.write('CONTROLS                 NO\n')
        tF.write('SUBCATCHMENTS            ALL\n')
        tF.write('NODES                    ALL\n')
        tF.write('LINKS                    ALL\n')


        tF.write('\n\n\n\n')
        tF.write('[MAP]\n')
        minx = round(min(nodo.p[0] for nodo in nodos),4)
        maxx = round(max(nodo.p[0] for nodo in nodos),4)
        miny = round(min(nodo.p[1] for nodo in nodos),4)
        maxy = round(max(nodo.p[1] for nodo in nodos),4)
        distx, disty = maxx - minx, maxy - miny
        tF.write(';;             minx           miny           maxx           maxy\n')
        tF.write(';;=========================================================================\n')
        list = ['DIMENSIONS', round(minx - 0.1 * distx,4), round(miny - 0.1 * disty,4), round(maxx + 0.1 * distx,4), round(maxy + 0.1 * disty,4)]
        tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
        tF.write('\nUNITS METERS\n')

        tF.write('\n')
        tF.write('[COORDINATES]\n')
        tF.write(';;Node         xcoord         ycoord         \n')
        tF.write(';;===========================================\n')
        for (i, nodo) in enumerate(nodos):
            list = ['NODO%d' % i, round(nodo.p[0], 4), round(nodo.p[1], 4)]
            tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            tF.write('\n')
        for (i, nodo) in enumerate(nodosOutfall):
            list = ['NODOOUT%d' % i, round(nodo.p[0], 4), round(nodo.p[1], 4)]
            tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
            tF.write('\n')

        tF.write('\n')
        tF.write('[VERTICES]\n')
        tF.write(';;LInk         xcoord         ycoord         \n')
        tF.write(';;===========================================\n')

        tF.write('\n')
        tF.write('[POLYGONS]\n')
        tF.write(';;Subcat       xcoord         ycoord         \n')
        tF.write(';;===========================================\n')
        for (i, subcuenca) in enumerate(subcuencas):
            for punto in subcuenca[0]:
                list = ['CUENCA'+str(i), np.round(punto[0], 4), np.round(punto[1], 4)]
                tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
                tF.write('\n')

        tF.write('\n')
        tF.write('[SYMBOLS]\n')
        tF.write(';;Gage         xcoord         ycoord         \n')
        tF.write(';;===========================================\n')
        list = ['GAGE1', 0.5*(minx+maxx), 0.5*(miny+maxy)]
        tF.write(('').join([ str(x).ljust(15, ' ') for x in list]))
        tF.write('\n')
    
    print('Fin write SWMMFile')

def mainPrepareDrainageNetwork(shpFileDrainageOriginal,shpFileDrainagePrepared, rasterFileDEM):
    streams = gpd.read_file(shpFileDrainageOriginal)
    crs_proy = streams.crs
    raster_dem = rasterio.open(rasterFileDEM)

    # Crea una columna con los pts inicial y final del trama
    streams['points'] = streams.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
    #streams.loc[:,'levelIni'] = None
    #streams.loc[:,'levelFin'] = None

    for idx, stream in streams.iterrows():
        # puntos del tramo
        
        pt_i = stream['points'][0]  # Pt Inicio
        pt_f = stream['points'][-1]  # Pt Fin
        lst_pts = [pt_i,pt_f]
        #streams.loc[idx,'Pt_i'] = pt_i
        #streams.loc[idx,'Pt_f'] = pt_f
        
        # Sample the terrain on the nodes
        lst_pts_elev = []
        for x in raster_dem.sample(lst_pts):
            lst_pts_elev.append(x)
        
        streams.loc[idx,'Elev_i'] = lst_pts_elev[0]
        streams.loc[idx,'Elev_f'] = lst_pts_elev[-1]

        transect = streams.loc[idx,'Transect']
        if transect is None:
            transect = ''
        streams.loc[idx,'Transect'] = str(transect)

        w = streams.loc[idx,'Ancho']
        if w is None:
            print ('ERROR: Missing w value on stream')
            #return
        h = streams.loc[idx,'Alto']
        if h is None:
            print ('ERROR: Missing h value on stream')
            #return
        typ = streams.loc[idx,'Tipo']
        if typ is None:
            print ('ERROR: Missing type value on stream')
            #return
        if str(typ).lower() in ['entubado', 'conducto', 'conduit']:
            typ = 'conduit'
        else:
            typ = 'channel'
        streams.loc[idx,'Tipo'] = typ
                
        levelIni = streams.loc[idx,'levelIni']
        levelFin = streams.loc[idx,'levelFin']
        depthIni = streams.loc[idx,'depthIni']
        depthFin = streams.loc[idx,'depthFin']

        if levelIni is None:
            if depthIni is None:
                if typ == 'conduit':
                    depthIni = h + params['minCoverage']
                else:
                    depthIni = h
            levelIni = streams.loc[idx,'Elev_i'] - depthIni
        else:
            depthIni = streams.loc[idx,'Elev_i'] - levelIni

        if levelFin is None:
            if depthFin is None:
                if typ == 'conduit':
                    depthFin = h + params['minCoverage']
                else:
                    depthFin = h
            levelFin = streams.loc[idx,'Elev_f'] - depthFin
        else:
            depthFin = streams.loc[idx,'Elev_f'] - levelFin

        streams.loc[idx,'depthIni'] = float(depthIni)
        streams.loc[idx,'depthFin'] = float(depthFin)
        streams.loc[idx,'levelIni'] = float(levelIni)
        streams.loc[idx,'levelFin'] = float(levelFin)

    streams.loc[:,'longitud'] = streams.length

    streams.loc[:,'slope'] = (streams['levelIni'] - streams['levelFin'])/streams['longitud']

    del streams['points']
    streams.to_file(shpFileDrainagePrepared, driver='ESRI Shapefile')

    # Create points along streams
    streams['points'] = streams.apply(lambda x: [y for y in x['geometry'].coords], axis=1)
    pointsAlongNetwork = []
    levelsAlongNetwork = []

    for idx, stream in streams.iterrows():
        lst_points = stream['points']
        length  = insertPoints(lst_points, 50)

        levelIni = streams.loc[idx,'levelIni']
        levelFin = streams.loc[idx,'levelFin']

        prog = 0
        levels = [float(levelIni)]
        for p0, p1 in pairwise(lst_points):
            prog += dist(p0,p1)
            levels.append(float(levelIni + prog / length * (levelFin - levelIni)))

        pointsAlongNetwork.extend(lst_points)
        levelsAlongNetwork.extend(levels)

    Data= {'CoorXY':pointsAlongNetwork,'levelBot':levelsAlongNetwork}
    NodesAlongNetwork = pd.DataFrame(Data)

    NodesAlongNetwork[['x','y']] = pd.DataFrame(NodesAlongNetwork.CoorXY.tolist(), index= NodesAlongNetwork.index)
    del NodesAlongNetwork['CoorXY']

    NodesAlongNetwork_shp = gpd.GeoDataFrame(NodesAlongNetwork, geometry=gpd.points_from_xy(NodesAlongNetwork.x,NodesAlongNetwork.y),crs=crs_proy)
    
    # Sample the terrain on the nodes
    coords = [(x,y) for x, y in zip(NodesAlongNetwork_shp.x, NodesAlongNetwork_shp.y)]
    NodesAlongNetwork_shp['levelTer'] = [x[0] for x in raster_dem.sample(coords)]
    NodesAlongNetwork_shp['depth'] = NodesAlongNetwork_shp['levelTer'] - NodesAlongNetwork_shp['levelBot']

    NodesAlongNetwork_shp.to_file(shpFileNodesNetworkDepth, driver='ESRI Shapefile')

def mainCreateSWMMModel(shpFileDrainagePrepared, shpFileCulverts,
            shpFileCalles, shpFileCuenca, rasterFileDEM, rasterFileSlope,
            rasterFileImpermeabilidad, shpFileNodosBorde, gageMethod, gageFileName,
            rasterFileCoeficiente, gagesFileName, stationsFileName, swmmInputFileName,
            rasterFileCurvaNumero, rasterTipoCurva, rasterFileSaSPerv,rasterFileNPerv,
            modelFolder,workspace):
        # Create nodes and links for the model
        nodos = NodesList()
        links = LinksDict(nodos)

        # Read spatial reference for the project
        shp_calles = gpd.read_file(shpFileCalles)
        spatial_ref = shp_calles.crs

         # Read the drainage network and create nodes and links
        readDrainageNetwork(nodos, links, shpFileDrainagePrepared)       

        ''' NO ESTABA HABILITADA 
        # Read culvert shape and and create culvert nodes and links
        ##readCulverts(nodos, links, shpFileCulverts)        
        '''

        # Read the street network and create nodes and links
        readStreets(nodos, links, shpFileCalles)

        # Create gutters and weirs
        createGutters(nodos, links)
        createWeirs(nodos, links)

        ''' NO ESTABA HABILITADA 
        # Read 2d zones and create nodes and links
        # generate2dZones(nodos, links, 'F:/Trabajo/Dropbox/Federico/ConuPy_version_27-4-2016/Ejemplo/Zona2D/zona2D.shp')
        '''

        # Calculate elevations
        calculateElevations(nodos, links, shpFileNodosSample, rasterFileDEM, rasterFileSlope, rasterFileImpermeabilidad, spatial_ref,rasterFileCurvaNumero, rasterTipoCurva, rasterFileSaSPerv,rasterFileNPerv)

        # Create outfalls
        nodosOutfall, lineasOutfall = createOutfallNodes(nodos, shpFileNodosBorde)

        # Calculate dead depths
        calculateDeadDepths(nodos, links, lineasOutfall)

        # Flip links
        linksList = [(n0, n1) for n0, n1 in links]
        for n0, n1 in linksList:
            link = links[(n0, n1)]
            if link['levelFin'] > link['levelIni']:
                del links[(n0, n1)]
                link['levelIni'], link['levelFin'] = link['levelFin'], link['levelIni']
                links[(n1, n0)] = link

        
        # Write the network
        writeNetworkShapes(nodos, links, shpFileNodos, shpFileLineas, spatial_ref)

        # Create subcatchments
        centros, subcuencas = createSubcatchments(nodos, shpFileCuenca, spatial_ref)
        # print(np.array(list(subcuencas[0][0].exterior.coords)))

        # Create raingages and map them to each subcatchment
        if gageMethod == 'createRainGagesMethod0':
            gages, subcatchmentGages = createRainGagesMethod0(centros, gageFileName, rasterFileCoeficiente, gagesFileName)
        else:
            gages, subcatchmentGages = createRainGagesMethod1(centros, stationsFileName)

        saveOnFile(workspace,nodos, 'nodos')
        saveOnFile(workspace,links, 'links')
        saveOnFile(workspace,centros, 'centros')
        saveOnFile(workspace,lineasOutfall, 'lineasOutfall')

        # Write the model file
        writeSWMMFile(nodos, links, centros, subcuencas, nodosOutfall, lineasOutfall, gages, subcatchmentGages, swmmInputFileName, modelFolder)

################################################################################################################################################################

# Captura la ruta
dataFolder = os.path.abspath(os.getcwd()) + '/Inputs_ConuPy/'
dataFolder = dataFolder.replace('\\', '/')

# Abre el csv con los nombres de los archivos y los parametros del modelo
archivo_parametros =  dataFolder + '/Parametros_conupy.csv'
df_params = pd.read_csv(archivo_parametros,index_col='variable',encoding = 'ISO-8859-1')

# Archivos de entrada
nf_Red = df_params.loc['nf_Red','valor']
shpFileDrainageOriginal     = dataFolder + '02_RedDeDrenaje/'+nf_Red+'.shp'
shpFileDrainagePrepared     = dataFolder + '03_RedDeDrenajePrep/'+nf_Red+'.shp'
rasterFileDEM               = dataFolder + '04_MDT/' + df_params.loc['rasterFileDEM','valor']

nf_calles = df_params.loc['nf_calles','valor']
shpFileCalles               = dataFolder + '01_Calles/'+nf_calles+'.shp'

shpFileCulverts             = dataFolder + 'Alcantarillas/' + df_params.loc['shpFileCulverts','valor']
shpFileCuenca               = dataFolder + '00_Cuenca/' + df_params.loc['shpFileCuenca','valor']
rasterFileSlope             = dataFolder + '05_Pendientes/' + df_params.loc['rasterFileSlope','valor']
rasterFileImpermeabilidad   = dataFolder + '06_Impermeabilidad/' + df_params.loc['rasterFileImpermeabilidad','valor']
rasterFileCoeficiente       = dataFolder + '08_Precipitacion/' + df_params.loc['rasterFileCoeficiente','valor']
shpFileNodosBorde           = dataFolder + '07_NodosBorde/' + df_params.loc['shpFileNodosBorde','valor']

# CN
rasterFileCurvaNumero   = dataFolder + "06_Impermeabilidad/" + df_params.loc['rasterFileCurvaNumero','valor']
rasterFileSaSPerv   = dataFolder + "06_Impermeabilidad/" + df_params.loc['rasterFileSaSPerv','valor']
rasterFileNPerv   = dataFolder + "06_Impermeabilidad/" + df_params.loc['rasterFileNPerv','valor']

# Urbano - Rural
rasterTipoCurva   = dataFolder + "06_Impermeabilidad/" + df_params.loc['rasterTipoCurva','valor']

# Directorio del modelo
modelFolder                 = dataFolder + '10_Modelo/'
swmmInputFileName           = modelFolder + df_params.loc['swmmInputFileName','valor']
defaultGageFileName         = modelFolder + df_params.loc['defaultGageFileName','valor']
defaultGagesFileName        = modelFolder + df_params.loc['defaultGagesFileName','valor']
swmmOuputFileName           = modelFolder + df_params.loc['swmmOuputFileName','valor']
gageMethod                  = df_params.loc['gageMethod','valor']
gageFileName                = modelFolder + df_params.loc['gageFileName','valor']
gagesFileName               = modelFolder + df_params.loc['gagesFileName','valor']
stationsFileName            = dataFolder +'08_Precipitacion/' + df_params.loc['stationsFileName','valor']

# workspace temporal
workspace                   = dataFolder + '09_WS/'

# Output files
shpFileNodesNetworkDepth    = workspace +'nodesNetworkDepth.shp'
shpFileNodosSample          = workspace +'nodosSample.shp'
shpFileNodos                = workspace +'nodos.shp'
shpFileLineas               = workspace +'lineas.shp'
shpFileCentros              = workspace +'centros.shp'
subcuencasShpFile           = workspace +'subcuencas.shp'
subcuencasClipShpFile       = workspace +'subcuencasClip.shp'

# Parametros del Modelo
target_Dx = float(df_params.loc['targetDx','valor'])
params = {
    'nodeTypesAsJunctions': ['conduit', '2d', 'channel'], #'corner'
    'maxDistSnapStreamNodes': float(df_params.loc['maxDistSnapStreamNodes','valor']),                          # Max dist for which channel and stream nodes are snapped together
    'targetDx': target_Dx,
    'maxLengthForStreamSpanDivide': target_Dx * float(df_params.loc['maxLengthForStreamSpanDivide','valor']),  # Max length stream spans can have without been subdivided
    'minLengthForStreamSpanJoin': target_Dx * float(df_params.loc['minLengthForStreamSpanJoin','valor']),      # Min length stream spans can have without been joined
    'maxLengthForStreetSpanDivide': target_Dx * float(df_params.loc['maxLengthForStreetSpanDivide','valor']),  # Max length street spans can have without been subdivided
    'minLengthForStreetSpanJoin': target_Dx * float(df_params.loc['minLengthForStreetSpanJoin','valor']),      # Min length street spans can have without been joined
    'minCoverage': float(df_params.loc['minCoverage','valor']),                                                # Default coverage for conduits when no depths or levels were defined
    'maxDistCulvert': float(df_params.loc['maxDistCulvert','valor']),                                          # Max distance from a stream required for a culvert to be created
    'maxDistGutter': target_Dx * float(df_params.loc['maxDistGutter','valor']),                                # Max distance required for gutters to be created from a corner to a conduit node
    'maxDistWeir': target_Dx * float(df_params.loc['maxDistWeir','valor']),                                    # Max distance required for weirs to be created from a corner to a channel segment
    'maxDistConnectOutfallNodes': float(df_params.loc['maxDistConnectOutfallNodes','valor']),                  # Max dist for which outfall nodes are connected to regular nodes
    'juApondPer': float(df_params.loc['juApondPer','valor']),                                                  # % of the basin area in which water is stored in a node
    'juY0': float(df_params.loc['juY0','valor']),                                                              # Water depth at start of simulation (ft or m) (default is 0).
    'evap': float(df_params.loc['evap','valor']),                                                              # mm/dia
    'saNImp': float(df_params.loc['saNImp','valor']),                                                          # Manning's n for overland flow over the impervious sub-area.
    'saNPerv': float(df_params.loc['saNPerv','valor']),                                                         # Manning's n for overland flow over the pervious sub-area.
    'saSImp': float(df_params.loc['saSImp','valor']),                                                           # depression storage for impervious sub-area (inches or mm). [mm]
    'saSPerv': float(df_params.loc['saSPerv','valor']),                                                         # depression storage for pervious sub-area (inches or mm). [mm]
    'saZero': float(df_params.loc['saZero','valor']),                                                           # percent of impervious area with no depression storage. [%]
    'inF0': float(df_params.loc['inF0','valor']),                                                               # Infiltration Horton  [mm/hr]
    'inFf': float(df_params.loc['inFf','valor']),                                                               # Infiltration Horton  [mm/hr]
    'inCoefDecaim': float(df_params.loc['inCoefDecaim','valor']),                                               # Infiltration Horton
    'Conductividad': float(df_params.loc['Conductividad','valor']),
    'inDryTime': float(df_params.loc['inDryTime','valor']),                                                     # Time it takes for fully saturated soil to dry  (days).
    'inMaxInf': float(df_params.loc['inMaxInf','valor']),                                                       # Maximum infiltration volume possible (0 if not applicable) (in or mm)
    'coN': float(df_params.loc['coN','valor']),                                                                 # Value of n (i.e., roughness parameter) in Manning's equation for conduits.
    'culvertRoadwayCd': float(df_params.loc['culvertRoadwayCd','valor']),                                       # Culvert parameters
    'culvertCode': float(df_params.loc['culvertCode','valor']),                                                 # Culvert parameters
    'culvertKentry': float(df_params.loc['culvertKentry','valor']),                                             # Culvert parameters
    'culvertKexit': float(df_params.loc['culvertKexit','valor']),                                               # Culvert parameters
    'weAlturaCordon': float(df_params.loc['weAlturaCordon','valor']),                                           # Weirs parameters
    'weCd': float(df_params.loc['weCd','valor']),                                                               # Weirs parameters
    'orCd': float(df_params.loc['orCd','valor']),                                                               # Orifice parameters
    'outfallXsWidth': float(df_params.loc['outfallXsWidth','valor']),                                           # XSection parameters
    'xsG1': float(df_params.loc['xsG1','valor']),                                                               # XSection parameters
    'xsG3': float(df_params.loc['xsG3','valor']),                                                               # XSection parameters
    'xsG4': float(df_params.loc['xsG4','valor']),                                                               # XSection parameters
    'xsSumideroH': float(df_params.loc['xsSumideroH','valor']),                                                 # XSection parameters
    'xsSumideroW':float(df_params.loc['xsSumideroW','valor']),                                                  # XSection parameters
    'xsVertederoH': float(df_params.loc['xsVertederoH','valor']),                                               # XSection parameters
    'xsVertederoW': float(df_params.loc['xsVertederoW','valor']),                                               # XSection parameters
    'traNConducto': float(df_params.loc['traNConducto','valor']),                                               # Transect parameters
    'traNArroyoPlanicie': float(df_params.loc['traNArroyoPlanicie','valor']),                                   # Transect parameters
    'traNArroyoCauce': float(df_params.loc['traNArroyoCauce','valor']),                                         # Transect parameters
    'traNCalle': float(df_params.loc['traNCalle','valor']),                                                     # Transect parameters
    'traAnchoMargenArroyo': float(df_params.loc['traAnchoMargenArroyo','valor']),                               # Transect parameters
    'traAnchoVereda': float(df_params.loc['traAnchoVereda','valor']),                                           # Transect parameters
    'traAltoCordon': float(df_params.loc['traAltoCordon','valor']),                                             # Transect parameters
    'traAltoCentroCalle': float(df_params.loc['traAltoCentroCalle','valor']),                                   # Transect parameters
    'cellSize2D': float(df_params.loc['cellSize2D','valor']),                                                   # 2d zones parameters
    'maxDist2DConnection': float(df_params.loc['maxDist2DConnection','valor']) * float(df_params.loc['cellSize2D','valor']),    # 2d zones parameters
    'numDiscreteGages': float(df_params.loc['numDiscreteGages','valor'])                                        # Number of discrete gages

    }

mainPrepareDrainageNetwork(shpFileDrainageOriginal, shpFileDrainagePrepared, rasterFileDEM)

mainCreateSWMMModel(shpFileDrainagePrepared, shpFileCulverts, shpFileCalles,
        shpFileCuenca, rasterFileDEM, rasterFileSlope, rasterFileImpermeabilidad,
        shpFileNodosBorde, gageMethod, gageFileName, rasterFileCoeficiente,
        gagesFileName, stationsFileName, swmmInputFileName,
        rasterFileCurvaNumero, rasterTipoCurva, rasterFileSaSPerv,rasterFileNPerv,
        modelFolder,workspace)
