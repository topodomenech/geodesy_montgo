# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 20:21:47 2016

@author: Vicent
"""

import math
import numpy as np
import scipy.optimize as optimization
from sympy import *
import copy
import itertools

def intial_coordinates_and_parameters():
    #Coordenadas actuales UTM
    vertex_utm_84 =[["Morella",409200.482,4572241.097,31,646.822],
                 ["Mont_Serrat",706923.489,4359905.813,30,330.106],
                 ["Montagut",476251.574,4610986.154,31,266.808],
                 ["Saint_Jean",421454.775,4580177.693,31,127.528],
                 ["Lleberia_Maybe",296088.314,4512810.419,31,146.959],
                 ["Mont_Sia",293431.056,4506444.148,31,341.979],
                 ["Prosch_de_Espina",277632.792,4528716.903,31,1231.119],
                 ["Tossal_de_Encanade",257344.048,4511992.791,31,1444.801],
                 ["Ares",743452.486,4483708.349,30,1373.263],
                 ["Desierto_de_las_Palmas",246865.005,4441459.997,31,786.629],
                 ["Espadan",723857.605,4420603.824,30,1094.790],
                 ["Cullera_Maybe",737464.022,4339919.082,30,285.072],
                 #["Torre_de_Cullera__Maybe",738505.363,4336252.883,30,64.843],
                 ["Montgo",250717.485,4298856.202,31,801.130],
                 ["La_Mola",372599.872,4280689.297,31,246.530],
                 ["Camp_Vell",357605.551,4324689.368,31,451.233],
                 ["Clop_de_Galazo_Maybe",452460.591,4386324.378,31,976.054]]
    
    #Ellipsoid GRS 80
    a = 6378137
    b = 6356752.314140347
    #Geometría del elipsoide
    '''
    Excentricidad_1 = math.sqrt(a**2+-b**2)/a
    aplanamiento = (a-b)/a
    '''    
    Excentricidad_2 = math.sqrt(a**2+-b**2)/b
    #Coordenades obtingudes dle llibre de 
    #Oeuvres completes de Francois Arago pag 98

    ### FALTA COMPROVAR QUE NO SIGUEN COORDENADES EN FORMATO DECIMAL... O ALGÚN ALTRE
    vertex_geograficas_1800 =  [["Morella",41.29696389,1.948,593.2],
                                ["Mont_Serrat",41.60597778,1.811,1237.2],
                                ["Montagut",41.40695,1.422,952.3],
                                ["Saint_Jean",41.13425278,1.352,85.7],
                                ["Lleberia_Maybe",41.09286667,0.864,918.1],
                                ["Mont_Sia",40.61442778,0.529,762.3],
                                ["Prosch_de_Espina",40.87995278,0.360,1178],
                                ["Tossal_de_Encanade",40.72314722,0.126,726.4],
                                ["Ares",40.46688611,-0.133,1317.6],
                                ["Desierto_de_las_Palmas",40.08625,0.030,726.4],
                                ["Espadan",39.9068,-0.383,1038.7],
                                ["Cullera_Maybe",39.17698056,-0.253,219.8],
                                ["Montgo",38.80750556,0.123,741.9],
                                ["La_Mola",38.66577778,1.534,966],
                                ["Camp_Vell",39.06001667,1.353,396.4],
                                #["Mont_Jou",41.36393611,2.052],
                                ["Clop_de_Galazo_Maybe",39.62166667,2.542]]
    '''
    achatamiento = 1/310
    Excentricidad_1_1880 = math.sqrt(a**2+-b**2)/a #0.9998185
    '''
    a_1800 = 6375653
    b_1800 = a_1800-(1/334.)*a_1800
    
    return a, b, Excentricidad_2, vertex_utm_84, a_1800, b_1800, vertex_geograficas_1800
    
def from_UTM_to_geograficas(vertex_utm_84, a, b, Excentricidad_2):
    radio_polar_curvatura = a**2/b
    vertex_geograficas_84=[]
    for i in vertex_utm_84:
        geograficas=[]
        X = i[1]
        Y = i[2]
        H = i[3]
        #Tratamiento previo de X e Y
        X=X-500000
        #Si el punto a tranformar estuviese en el hemisferio sur:
        #Y=Y-10000000
        # Cálculo del meridiano central
        landa0=H*6-183
        #Ecuaciones de Coticchia-Surace para el Problema Inverso
        latitud_prima = Y/(6366197.724*0.9996)
        v = radio_polar_curvatura*0.9996/(1+Excentricidad_2**2*math.cos(latitud_prima)**2)**0.5
        a = X/v
        A1 = math.sin(2*latitud_prima)
        A2 = A1*math.cos(latitud_prima)**2
        J2 = latitud_prima+A1/2
        J4 = (3*J2+A2)/4
        J6 = (5*J4+A2*math.cos(latitud_prima)**2)/3
        alfa = 3*Excentricidad_2**2/4
        beta = 5*alfa**2/3
        fi = 35*alfa**3/27
        B= 0.9996*radio_polar_curvatura*(latitud_prima-alfa*J2+beta*J4-fi*J6)
        b=(Y-B)/v
        bicho1 = Excentricidad_2**2*a**2*math.cos(latitud_prima)**2/2
        bicho2 = a*(1-bicho1)
        bicho3 = b*(1-bicho1)+latitud_prima
        bichos4 = (np.exp(bicho2)-np.exp(-bicho2))/2
        inc_landa = math.atan(bichos4/math.cos(bicho3))
        bicho5 = math.atan(math.cos(inc_landa)*math.tan(bicho3))
        longitud = inc_landa*180/math.pi +landa0
        latitud_rad = latitud_prima+(1+Excentricidad_2**2*math.cos(latitud_prima)**2-3*Excentricidad_2**2*math.sin(latitud_prima)*math.cos(latitud_prima)*(bicho5-latitud_prima)/2)*(bicho5-latitud_prima)
        latitud = latitud_rad*180/math.pi
        geograficas=[i[0],latitud,longitud,i[4]]
        vertex_geograficas_84.append(geograficas)
    return vertex_geograficas_84
    
    
def from_geod_to_geoc(vertex_geografica,a,b):
    # Paso de geodésicas a geocéntricas en datum origen.
    vertex_geocentricas=[]
    for i in vertex_geografica:
        lat = i[1]*math.pi/180
        lon = i[2]*math.pi/180
        h = i[len(i)-1]
        N = a*a/(math.sqrt((a*a*math.cos(lat)*math.cos(lat))+(b*b*math.sin(lat)*math.sin(lat))))
        X = (N+h)*math.cos(lat)*math.cos(lon)
        Y = (N+h)*math.cos(lat)*math.sin(lon)
        Z = ((b*b/(a*a))*N+h)*math.sin(lat)
        vertex_geocentricas.append([i[0],X,Y,Z])
    return vertex_geocentricas


def vcross(v):
    x, y, z = v
    mat = np.zeros((3,3))
    mat[0] = [ 0, -z,  y]
    mat[1] = [ z,  0, -x]
    mat[2] = [-y,  x,  0]
    return mat

def block(v):
    return np.hstack((np.eye(3), -vcross(v), v[:, np.newaxis]))
    
#vertex_geocentricas_1800=vertex_geocentricas_1800_subset
#vertex_geocentricas_84=vertex_geocentricas_84_subset
def func1(vertex_geocentricas_1800, vertex_geocentricas_84,lista_subset, lista_conjutos, lista_error, lista_parametros, cx, cy,cz,s,rx,ry,rz):
    lista1 = []
    lista2 = []
    for i in range(len(vertex_geocentricas_1800)):
        X= vertex_geocentricas_1800[i][1]
        Y= vertex_geocentricas_1800[i][2]
        Z= vertex_geocentricas_1800[i][3]
        X_b= vertex_geocentricas_84[i][1]
        Y_b= vertex_geocentricas_84[i][2]
        Z_b= vertex_geocentricas_84[i][3]
        lista1.append([X,Y,Z])
        lista2.append([X_b,Y_b,Z_b])
    
    pt2=np.array(lista1)
    pt1=np.array(lista2)


    A = []
    rhs = []
    for i in range(3):
        A.append(block(pt1[i]))
        rhs.append((pt2[i] - pt1[i])[:, np.newaxis])
    
    A = np.vstack(A)
    rhs = np.vstack(rhs)
    
    ttt = np.linalg.lstsq(A, rhs)[0] # parameters

#quoated pdf results

#pdf = np.array([-9.256,-23.701,16.792,-0.0001990982,0.0001778762,0.00015,0.00000046])
    
    total = 0
    for i in range(len(vertex_geocentricas_1800)):
        X= vertex_geocentricas_1800[i][1]
        Y= vertex_geocentricas_1800[i][2]
        Z= vertex_geocentricas_1800[i][3]
        X_b= vertex_geocentricas_84[i][1]
        Y_b= vertex_geocentricas_84[i][2]
        Z_b= vertex_geocentricas_84[i][3]
        
        
        X_ = (ttt[0] +(1+ttt[6]*1e-6)*(X_b-ttt[5]*Y_b+ttt[4]*Z_b))
        Y_ = (ttt[1] +(1+ttt[6]*1e-6)*(ttt[5]*X_b+Y_b-ttt[3]*Z_b))
        Z_ = (ttt[2] +(1+ttt[6]*1e-6)*(-ttt[4]*X_b+ttt[3]*Y_b+Z_b))
             
        total += (X-X_)**2+(Y-Y_)**2+(Z-Z_)**2
    error = math.sqrt(int(total))
    if error<20:
        print error
        print lista_subset
    lista_error.append(error)
    lista_conjutos.append(lista_subset)
    lista_parametros.append(ttt)
    return lista_error, lista_conjutos, lista_parametros

#conjunto_comprobacion =['Desierto_de_las_Palmas', 'Espadan', 'Cullera_Maybe', 'Montgo', 'Camp_Vell']
#ttt = [-3.12657067e+03,-9.94553511e+02,  1.07933898e+03, -4.80442441e-04,  7.62315536e-05, -2.19644782e-04,  1.65212866e-05]
def comprobacion(vertex_geocentricas_1800,vertex_geocentricas_84, conjunto_comprobacion, parametros_comprobacion):
    for i in range(len(vertex_geocentricas_1800)):
        if vertex_geocentricas_1800[i][0] in conjunto_comprobacion:
            X= vertex_geocentricas_1800[i][1]
            Y= vertex_geocentricas_1800[i][2]
            Z= vertex_geocentricas_1800[i][3]
            X_b= vertex_geocentricas_84[i][1]
            Y_b= vertex_geocentricas_84[i][2]
            Z_b= vertex_geocentricas_84[i][3]
            X_ = (ttt[0] +(1+ttt[6]*1e-6)*(X_b-ttt[5]*Y_b+ttt[4]*Z_b))
            Y_ = (ttt[1] +(1+ttt[6]*1e-6)*(ttt[5]*X_b+Y_b-ttt[3]*Z_b))
            Z_ = (ttt[2] +(1+ttt[6]*1e-6)*(-ttt[4]*X_b+ttt[3]*Y_b+Z_b))            
            print "X: {} -- {}//dif = {}".format(X,X_,X-X_)         
            print "Y: {} -- {}//dif = {}".format(Y,Y_,Y-Y_)         
            print "Z: {} -- {}//dif = {}".format(Z,Z_,Z-Z_)   
            print math.sqrt((X-X_)**2+(Y-Y_)**2+(Z-Z_)**2)
            print "--------------------------------"
    
def main():
    a, b, Excentricidad_2, vertex_utm_84, a_1800, b_1800, vertex_geograficas_1800 = intial_coordinates_and_parameters()
    vertex_geograficas_84 = from_UTM_to_geograficas(vertex_utm_84, a, b, Excentricidad_2)
    vertex_geocentricas_84 = from_geod_to_geoc(vertex_geograficas_84,a,b)
    vertex_geocentricas_1800 = from_geod_to_geoc(vertex_geograficas_1800,a_1800,b_1800)
    
    cx = Symbol('cx')
    cy = Symbol('cy')
    cz = Symbol('cz')
    s = Symbol('s')
    rx = Symbol('rx')
    ry= Symbol('ry')
    rz= Symbol('rz')
    
    lista_error = []
    lista_conjutos = []  
    lista_parametros = []    
    
    for L in range(1, len(vertex_geocentricas_84)+1):
        for subset in itertools.combinations(vertex_geocentricas_84, L):
            lista_subset=[]
            if len(subset)>4:
                vertex_geocentricas_84_subset = list(subset)
                vertex_geocentricas_1800_subset=[]
                for i in vertex_geocentricas_84_subset:
                    lista_subset.append(i[0])
                for i in vertex_geocentricas_1800:
                    if i[0] in lista_subset:
                        vertex_geocentricas_1800_subset.append(i)
                lista_error, lista_conjutos, lista_parametros =func1(vertex_geocentricas_1800_subset, vertex_geocentricas_84_subset, lista_subset, lista_conjutos, lista_error, lista_parametros, cx, cy,cz,s,rx,ry,rz)
    print lista_error[np.argmin(lista_error)]
    print lista_conjutos[np.argmin(lista_error)]
    print lista_parametros[np.argmin(lista_error)]
    #conjunto_comprobacion = lista_conjutos[np.argmin(lista_error)] 
    #parametros_comprobacion = lista_parametros[np.argmin(lista_error)][0]
    #comprobacion(vertex_geocentricas_1800,vertex_geocentricas_84, conjunto_comprobacion, parametros_comprobacion)
#https://mail.scipy.org/pipermail/scipy-user/attachments/20100302/113cb4b3/attachment-0001.py
if __name__ == '__main__':
    main()