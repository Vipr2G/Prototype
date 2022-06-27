#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 20 22:06:07 2020

@author: russelmiller
"""
import numpy as np
from math import sin,cos,atan2,pi,sqrt
from  scipy.linalg import eig
dlat2nmi = 60.0
I = np.diag([1,1])

def z_score(lat,lon,sma,smi,orient,slat,slon,ssma,ssmi,sorient,flat_earth=True):
    # All inputs in DB units
    P = cov(ssma,ssmi,sorient)
    R = cov(sma,smi,orient)
    if flat_earth == True:
         # have to convert X and Y to nmi ("flat-earth" based on site lat)
         X = np.array([[slat * dlat2nmi],[slon * dlon2nmi(slat)]])
         Y = np.array([[lat * dlat2nmi],[lon * dlon2nmi(slat)]])
         delta = np.subtract(Y,X)

         #Y = np.array([[lat],[lon]])
         #X = np.array([[slat],[slon]])
         #delta = XYdeg2delta_nmi(X,Y)  #Y - X = INT-SITE (X-mu) move mu (site) to (0,0)
    else:
        lat_nmi = lat*dlat2nmi
        lon_nmi = dlon2nmi(lat)*lon
        slat_nmi = slat*dlat2nmi
        slon_nmi = dlon2nmi(slat)*slon
        delta = np.array([[lat_nmi-slat_nmi],[lon_nmi-slon_nmi]])
    delta_transp = np.ndarray.transpose(delta)
    covmx = np.add(P,R)
    inv_cov = np.linalg.inv(covmx)
    z2 = np.dot(delta_transp,np.dot(inv_cov,delta))
    z = sqrt(z2)
    return z

def kf_update(rlat,rlon,rsma,rsmi,rorient,plat,plon,psma,psmi,porient,phc,hc_weighting,w):
    # KF update based on intercept and site data, all inputs in DB units
    #   rlat,..rorient - intercept parms
    #   plat,...porient - site parms
    #   phc = site heard count
    #   hc_weighting (logical) True--> added emphasis on site location
    #   w (scalar) - KF noise parameter
    ###################################################################
    
    # prepare inputs for KF update
    P = cov(psma,psmi,porient)
    R = cov(rsma,rsmi,rorient)

    # have to convert X and Y to nmi ("flat-earth" based on site lat)
    X = np.array([[plat * dlat2nmi],[plon * dlon2nmi(plat)]])
    Y = np.array([[rlat * dlat2nmi],[rlon * dlon2nmi(plat)]])
    
    # execute the KF update
    K = np.dot(P,np.linalg.inv(np.add(P,R)))
    if hc_weighting == True and phc > 0:
        X = X + (1/phc)*np.dot(K,np.subtract(Y,X))
    else:
        X = X + np.dot(K,np.subtract(Y,X))
    P = np.dot(np.subtract(I,K),P)
    P = P + noise_gen(w,P)
    
    # extract new site parameters
    newlat = X[0,0]/dlat2nmi
    newlon = X[1,0]/dlon2nmi(newlat)
    newsma,newsmi,neworient = get_ellipse(P)
    return newlat,newlon,newsma,newsmi,neworient
    
def cov(sma,smi,orient,double=False):
    axes = np.diag([sma,smi])/2
    if double == True:
        axes = axes * 2
    T = rot_mx(orient)
    s_dev = np.dot(T,axes)
    cov = np.dot(s_dev,np.ndarray.transpose(s_dev))
    return cov
    
def rot_mx(orient):  #returns transformation matrix for rotation orient (deg)
    orient = np.radians(orient)
    T = np.array([[cos(orient),-sin(orient)],[sin(orient),cos(orient)]])
    return T

def dlon2nmi(lat):
    # Conversion factor to convert degrees longitude to nmi
    #    Input is latitude in decimal degrees
    #    Example use:  lon_nmi = lon_deg * dlon2nmi(lat_deg)
    ########################################################
    
    lat = np.radians(lat)
    Re = 6378 - 21*np.sin(lat)
    B15 = min(1,np.abs(np.cos(lat)*np.sin(np.radians(1/120))))
    B16 = 2 * np.arcsin(B15)
    d = B16*Re/1.852 #dist in nmi for 1 minute of arc
    return d*60  #converts lon (deg) --> nmi

def noise_gen(w,P):
    pn_ratio = 1
    sma,smi,orient = get_ellipse(P)
    T = rot_mx(orient)
    Tinv = np.linalg.inv(T)
    p2x2 =np.diag([w,pn_ratio*w])
    p2x2 = np.dot(T,np.dot(p2x2,Tinv))
    return p2x2

def get_ellipse(P):
    # Given cov matrix in lat/lon ref, Return SMA, SMI, Orient
    
    diag = eig(P)
    eig_vals = np.real(diag[0])
    eig_vecs = np.real(diag[1])

    s1 = eig_vals[0]**0.5
    s2 = eig_vals[1]**0.5
    theta = atan2(eig_vecs[1,0],eig_vecs[0,0])*180/pi
    if s1 > s2 :
        newsig1 = s1
        newsig2 = s2
        if theta < 0:
            neworient = theta + 180
        else:
            neworient = theta
    else:
        newsig1 = s2
        newsig2 = s1
        if theta < 90:
            neworient = theta + 90
        else:
            neworient = theta - 90
    if neworient < 0:
        neworient = neworient + 180
    
    newsma = 2 * newsig1
    newsmi = 2 * newsig2
    return newsma,newsmi,neworient

def XYdeg2delta_nmi(x,y):
    # "Flat-Earth" delta (nmi) from 2 lat,lon column vectors in degrees
    avg_lat = (x[0,0] + y[0,0])/2
    delta_lat_nmi = (y[0,0]-x[0,0])*dlat2nmi
    delta_lon_nmi = dlon2nmi(avg_lat)*(y[1,0]-x[1,0])
    delta_nmi = np.array([[delta_lat_nmi],[delta_lon_nmi]])
    return delta_nmi

def to_nmi(x):
    #Converts 2-dim lat/lon array from decimal degrees to nmi
    #   convert arrays independently for non-"Flat-Earth" computations
    lat_nmi = x[0,0]*dlat2nmi
    k = dlon2nmi(x[0,0])
    lon_nmi = k*x[1,0]
    y = np.array([[lat_nmi],[lon_nmi]])
    return y

def to_deg(x):
    # converts 2-dim lat/lon vector from nmi to degrees
    lat_deg = x[0,0]/dlat2nmi
    lon_deg = float(x[1,0]/dlon2nmi(lat_deg))
    y = np.array([[lat_deg],[lon_deg]])
    return y

def convolver(latar,lonar,smaar,smiar,orientar):
    # Center point
    v = np.array([np.sum(np.sin(orientar*pi/180)), np.sum(np.cos(orientar*pi/180))])
    u = v/np.sqrt(np.dot(v,v))
    conv_lon = np.sum(lonar/smaar)/np.sum(1/smaar)
    conv_lat = np.sum(latar/smaar)/np.sum(1/smaar)
    conv_sma = np.dot(u,v)/np.sum(1/smaar)
    conv_smi = np.dot(u,v)/np.sum(1/smiar)
    conv_orient = np.arctan2(u[0],u[1])*180/pi
    return conv_lat,conv_lon,conv_sma,conv_smi,conv_orient

def geo2dist(lat1,lon1,lat2,lon2):
    lat1r = lat1*pi/180
    lon1r = lon1*pi/180
    lat2r = lat2*pi/180
    lon2r = lon2*pi/180
    Re = 6378- 21 * np.sin((lat1r + lat2r)/2)
    B15 =( ( np.sin( 0.5*(lat1r-lat2r) ) )**2 + np.cos(lat1r)*np.cos(lat2r)*( np.sin( .5*(lon1r-lon2r) ) )**2 )**0.5
    if B15 < 1:
        B16 = 2*np.arcsin(B15)
    else:
        B16 = pi
    d = B16*Re
    return d



    


   
