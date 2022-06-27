#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 31 09:00:18 2021

@author: russelmiller
"""
import numpy as np
from VL2db_utils import db_con, db_geo_insert, tune_geo
from math import pi, trunc, sin, cos, sqrt
from numpy import arcsin, sign, arctan2
from VL2shared_utils import feql,rip,add_unique
from VL2config import fq_idb, axis_min_nm, min_temp_pct, min_geo_corr_ints  #, max_SMA
from VL2geo_associator import associate_geo
from datetime import datetime
nmi2k = 1.852

def proc_geo(elnot):
    print("starting geo miner")
    start = datetime.now()
    db = db_con()
    cur = db.cursor()
    # pre-associate data
    cur.execute("delete from comparison_geo_snapshot")
    cur.execute("insert into comparison_geo_snapshot select intercept_id from " +fq_idb+ """intercepts where elnot = %s and intercept_id in
                (select intercept_id from residue_ints where r_type = 'G')""",(elnot,))
    cur.execute("commit")
    associate_geo(False)
    
    
    
    
    #cur.execute("""select a.intercept_id, latitude,longitude,area_dist_maj,area_dist_min,orientation from """ + fq_idb + """intercepts a, residue_ints b where
    #               a.intercept_id = b.intercept_id and a.elnot = %s and b.r_type = 'G'
    #               and a.latitude is not null and a.longitude is not null and a.area_dist_maj is not null and a.area_dist_min is not null
    #               and a.latitude >= 0 and a.longitude >= 0 and a.area_dist_maj > 0 and a.area_dist_min > 0 and a.orientation >= 0 
    #                       and a.orientation is not null and a.is_emitter = 'Y'""",(elnot,))
    cur.execute("""select a.intercept_id, latitude,longitude,area_dist_maj,area_dist_min,orientation from """ + fq_idb + """intercepts a, residue_ints b where
                   a.intercept_id = b.intercept_id and a.elnot = %s and b.r_type = 'G'
                   and a.latitude is not null and a.longitude is not null and a.area_dist_maj is not null and a.area_dist_min is not null
                   and a.area_dist_maj > 0 and a.area_dist_min > 0 and a.orientation >= 0 
                           and a.orientation is not null and a.is_emitter = 'Y'""",(elnot,))
    heat_input = cur.fetchall()
    cur.close()
    db.close()
    if len(heat_input)<=1:
        print("No data in residue....exiting proc_geo")
        return [],[]
    # repace any null orientations with 0.0
    #heat_input = repl_nulls(heat_input,5,0.0) -- added coalesce to sql query instead
    
    next_id = 1   # one-up index on new sites...db_geo_insert() will reconcile new site_ids with db
    
    grid_res,roa = tune_geo(elnot)
    heat_summary,heat_out = heat_map_gen(heat_input,grid_res)
    site_summary = heat_hunter(heat_summary)
    new_sites,new_site_ints = auto_convolve(site_summary,heat_out,next_id)
    db_geo_insert(elnot,new_sites,new_site_ints)
    end = datetime.now()
    #plot_all('sh')
    print("geo_miner start: ",start)
    print("geo_miner end: ",end)
    print("geo_miner total:",end-start)
    
    return new_sites,new_site_ints
    
def heat_map_gen(heat_input,grid_res):
    
    pad_nmi = 5.0 #nmi pad for the processing grid

    lats = rip(heat_input,1)
    lons = rip(heat_input,2)
    min_lat = float(min(lats))
    max_lat = float(max(lats))
    min_lon = float(min(lons))
    max_lon = float(max(lons))
    
    min_lat,max_lat,min_lon,max_lon = get_grid(pad_nmi,min_lat,max_lat,min_lon,max_lon)
    
    
    hit = 0
    km2lat_deg = (max_lat-min_lat)/geo2dist(min_lat,min_lon,min_lat,max_lon)
    delta_lat = grid_res*km2lat_deg
    avg_lat = (min_lat+max_lat)/2
    km2lon_deg = 1/geo2dist(avg_lat,min_lon,avg_lat,min_lon+1)
    delta_lon = km2lon_deg*grid_res
    Re = 6378-21*sin(((max_lat + min_lat)/2)*pi/180)  #converted to radians
    N_lat = trunc((max_lat-min_lat)/delta_lat) + 1
    N_lon = trunc((max_lon-min_lon)/delta_lon) + 1

  
    # heat_input is array of intercept_id,latitude,longitude,area_dist_maj,area_dist_min,orientation
    heat_out = []
        
    for row in heat_input:
        int_id = row[0]
        lat1 = float(row[1]) * pi/180
        lon1 = float(row[2]) * pi/180
        phi1 = float(row[5]) * pi/180
        SMA1 = max(axis_min_nm, float(row[3]) * nmi2k)
        SMI1 = max(axis_min_nm, float(row[4]) * nmi2k)
        d2Nlon = geo2dist(lat1,lon1,lat1,lon1 + 10*delta_lon)/10
        d2Nlat = geo2dist(lat1,lon1,lat1 + 10*delta_lat,lon1)/10
        dNN = SMA1*abs(cos(phi1))+SMI1*abs(sin(phi1))
        dEE = SMA1*abs(sin(phi1)) + SMI1*abs(cos(phi1))

        lon1a = lon1-(dEE/d2Nlon*delta_lon) * pi/180
        lon1b = lon1+(dEE/d2Nlon*delta_lon) * pi/180
        NLON_start = max(0,trunc((lon1a*180/pi-min_lon)/delta_lon))
        NLON_stop = min(N_lon,trunc((lon1b*180/pi-min_lon)/delta_lon +1))

        lat1a = lat1-(dNN/d2Nlat*delta_lat) * pi/180
        lat1b = lat1+(dNN/d2Nlat*delta_lat) * pi/180
        NLAT_start = max(0,trunc((lat1a*180/pi-min_lat)/delta_lat))
        NLAT_stop = min(N_lat,trunc((lat1b*180/pi-min_lat)/delta_lat +1))
        for j in range(NLON_start, NLON_stop + 1):
            lon2 = (min_lon + j*delta_lon)*pi/180
            hit =0
            for k in range(NLAT_start, NLAT_stop + 1):
                lat2 = (min_lat + k*delta_lat)*pi/180
                a1 = (sin((lat1-lat2)/2))**2 + cos(lat1)*cos(lat2)*(sin((lon1-lon2)/2))**2
                #sqrt_a = sqrt(a1)
                c1 = 2*arcsin(min(1,sqrt(a1)))
                d = c1*Re
                dn = sign(lat2-lat1)*Re*2*arcsin(min(1,sqrt((sin((lat2-lat1)/2))**2)))
                B25 = ((cos(lat1)*sin((lon1-lon2)/2)))**2
                B26 = 2*arcsin(min(1,sqrt(B25)))
                de = sign(lon2-lon1)*Re*B26
                if feql(de,0) and feql(dn,0):  # THEN --28 May mod
                    heat_out.append([int_id,j,k])
                    #insert into heat_out(intercept_id,N_lon,N_lat) values(ellipses.intercept_id,j,k); -- 28 May mod
                    hit = 1
                else:
                    gamma = arctan2(de,dn)  #note: convention is opposite of Excel

                    beta1 = gamma-phi1

                    r1 = SMA1*SMI1/sqrt( ( (SMI1*cos(beta1)))**2 + ( (SMA1*sin(beta1)))**2 )
                    if r1 > d:
                        heat_out.append([int_id,j,k])
                        #insert into heat_out(intercept_id,N_lon,N_lat) values(ellipses.intercept_id,j,k);
                        hit = 1  #--????? Addt'l risky mod...if probs on testing remove for debug
                    elif hit == 1:
                        continue
    #Populate heat_summary
    indx_temps = []
    for row in heat_out:
        N_lon = row[1]
        N_lat = row[2]
        indx_temps = add_unique([N_lon,N_lat],indx_temps)
    heat_summary = []
    for row in indx_temps:
        N_lon = row[0][0]
        N_lat = row[0][1]
        temp = row[1]
        lat = min_lat + N_lat * delta_lat
        lon = min_lon + N_lon *  delta_lon
        heat_summary.append([lat,lon,N_lat,N_lon,temp])
    
    return heat_summary, heat_out

def heat_hunter(heat_summary):
    #if heat_summary == []:
    #   return []
    site_summary = []
    #TYPE grid_index_type is table of INTEGER;
    #lon_index_var grid_index_type;
    s_id = 1

    temps = rip(heat_summary,4)
    if temps ==[]:
        max_temp = 0
    else:
        max_temp = max(temps)
    min_temp = max(min_temp_pct * max_temp, min_geo_corr_ints)
    if min_temp > max_temp:
        print("No temps above threshold....exiting")
        return[]

    for row in heat_summary.copy():
        if row[4] <= min_temp:
            heat_summary.remove(row)
    if len(heat_summary) == 0:
        print("HEAT SUMMARY EMPTY....RETURNING NULL")
        return []
    n_lons = rip(heat_summary,3)  #.sort()  here returns NoneType...
    n_lons.sort()
    #nlats = rip(heat_summary,2).sort()
    
    

    print('temp = ',min_temp, max_temp)

    grp_min = min(n_lons)
    grp_end = max(n_lons)
    prev_indx = grp_min - 1
    lon_groups = []
    k = 1
    for lon in n_lons:
        curr_indx = lon
        if curr_indx > prev_indx + 1:
            grp_max = prev_indx
            lon_groups.append([k,grp_min,grp_max])
            grp_min = curr_indx
            k +=1
        elif curr_indx == grp_end and curr_indx-prev_indx == 1:
            grp_max = grp_end
            lon_groups.append([k, grp_min,grp_max])
        if curr_indx == grp_end and curr_indx - prev_indx > 1:
            grp_max = grp_end
            lon_groups.append([k, grp_min,grp_max])
        prev_indx = curr_indx


    k = 1
    for lon_grp in lon_groups:
        #grp_id = lon_grp[0]
        min_nlon = lon_grp[1]
        max_nlon = lon_grp[2]
        nlats = []
        for row in heat_summary:
            hs_lon = row[3]
            if hs_lon >= min_nlon and hs_lon <= max_nlon:
                n_lat = row[2]
                nlats = add_unique(n_lat,nlats)
                #nlats.append(n_lat)
        nlats = rip(nlats,0)
        nlats.sort()
        grp_min = min(nlats)
        grp_end = max(nlats)
        prev_indx = grp_min-1
        for lat in nlats:
            curr_indx = lat
            #grp_min = 0
            #grp_max = 0
            if curr_indx > prev_indx + 1:
                grp_max = prev_indx
                #lat_groups.append([s_id,grp_id,k,grp_min,grp_max])
                site_summary.append([s_id,min_nlon,max_nlon,grp_min,grp_max])
                s_id +=1
                grp_min = curr_indx
                k +=1
            elif curr_indx == grp_end and (curr_indx - prev_indx == 1 or curr_indx - prev_indx == 0):
                grp_max = grp_end
                #lat_groups.append(s_id,grp_id,k,grp_min,grp_max)
                site_summary.append([s_id,min_nlon,max_nlon,grp_min,grp_max])
                s_id +=1
                k +=1
            if curr_indx == grp_end and curr_indx-prev_indx >1:
                grp_min = curr_indx
                grp_max = grp_end
                #lat_groups.append([s_id,grp_id,k,grp_min,grp_max])
                site_summary.append([s_id,min_nlon,max_nlon,grp_min,grp_max])
                s_id +=1
                k +=1
            prev_indx = curr_indx
            #if grp_min != 0 and grp_max !=0:
            #site_summary.append([s_id,min_nlon,max_nlon,grp_min,grp_max])
    return site_summary
    # End Heat_hunter

def auto_convolve(site_summary, heat_out, next_id):
    db = db_con()
    cur = db.cursor()
    new_sites = []
    new_site_ints = []
    for row in site_summary:
            #for row in site_summary:
            #s_id = row[0]
            min_lon_idx = row[1]
            max_lon_idx = row[2]
            min_lat_idx = row[3]
            max_lat_idx = row[4]

            
            site_ints = []
            for int in heat_out:
                int_id = int[0]
                int_lon_idx = int[1]
                int_lat_idx = int[2]
                if between(int_lon_idx,min_lon_idx,max_lon_idx) and between(int_lat_idx,min_lat_idx,max_lat_idx):
                    site_ints = add_unique(int_id,site_ints)
            site_ints = rip(site_ints,0)  # unique ints contributing
            
            int_str = str(tuple(site_ints))
            cur.execute("""select latitude,longitude,area_dist_maj,area_dist_min,orientation from intercepts where intercept_id in """ + int_str)
            dat = cur.fetchall()
            latar = np.array(rip(dat,0)).astype(float)
            lonar = np.array(rip(dat,1)).astype(float)
            smaar = np.array(rip(dat,2)).astype(float)
            smiar = np.array(rip(dat,3)).astype(float)
            orientar = np.array(rip(dat,4)).astype(float)
            slat,slon,ssma,ssmi,sorient = convolver(latar,lonar,smaar,smiar,orientar)
            new_sites.append([next_id,slat,slon,ssma,ssmi,sorient])
            for int_id in site_ints:
                new_site_ints.append([next_id, int_id])
            next_id +=1
    cur.close()
    db.close()
    return new_sites,new_site_ints
            
            
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


def get_grid(pad_nmi,min_lat,max_lat,min_lon,max_lon):
    # pads the lat/lon box by pad_nmi (nmi)
    # lat/lon inputs and outputs all in deg
    
    nmi2k = 1.852
    
    km2lat_deg = (max_lat-min_lat)/geo2dist(min_lat,min_lon,max_lat,min_lon)
    avg_lat = (min_lat+max_lat)/2
    km2lon_deg = 1/geo2dist(avg_lat,min_lon,avg_lat,min_lon+1)

    min_lon = min_lon - nmi2k*pad_nmi*km2lon_deg
    max_lon = max_lon + nmi2k*pad_nmi*km2lon_deg
    min_lat = min_lat - nmi2k*pad_nmi*km2lat_deg
    max_lat = max_lat + nmi2k*pad_nmi*km2lat_deg

    if (min_lon < -179.99):
       min_lon = -179.99

    if (max_lon > 179.99):
        max_lon = 179.99
        
    return min_lat,max_lat,min_lon,max_lon
    ## End get_grid
    

def between(x,x_min,x_max):
    if x_min <= x and x_max >= x:
        return True
    else:
        return False
    

def geo2dist(lat1,lon1,lat2,lon2):
    lat1r = float(lat1)*pi/180
    lon1r = float(lon1)*pi/180
    lat2r = float(lat2)*pi/180
    lon2r = float(lon2)*pi/180
    Re = 6378- 21 * sin((lat1r + lat2r)/2)
    B15 =( ( sin( 0.5*(lat1r-lat2r) ) )**2 + cos(lat1r)*cos(lat2r)*( sin( .5*(lon1r-lon2r) ) )**2 )**0.5
    if B15 < 1:
        B16 = 2*arcsin(B15)
    else:
        B16 = pi
    d = B16*Re
    return d
