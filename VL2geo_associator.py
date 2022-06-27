#import psycopg2
from VL2db_utils import db_con
from VL2KF_utils import z_score, kf_update
from datetime import datetime
from math import sqrt
from scipy.stats import norm
#from GeoViewer3 import plot_all
from VL2config import LOC, track_deltas, save_tracks, w, hc_weighting, hc_max, use_sqrt, fq_idb, max_SMA


def associate_geo(pop_residue = True):
    print("starting geo association")
    start = datetime.now()

    crit_val = norm.ppf((1+LOC)/2) # z-score
    print("LOC,crical value",LOC,crit_val)


    # Initialize:
    connection = db_con()

    intCur = connection.cursor()
    siteCur = connection.cursor()
    parmCur = connection.cursor()
    #out_x = []
    #out_y=[]
    if track_deltas == True:
        siteCur.execute("delete from site_history")
        siteCur.execute("insert into site_history(site_id,lat,lon,sma,smi,orient) select site_id,lat,lon,sma,smi,orient from site_geos")
        siteCur.execute("update site_history set iter = 0")
    if save_tracks:
        siteCur.execute("delete from residue_ints where r_type='G'")
    
    # safeguard against null heard counts
    siteCur.execute("update site_geos set heard_count = 0 where heard_count is null")

    # Loop through intercepts
    #intCur.execute("delete from site_ints")
    #intCur.execute("delete from site_ints2")
    
    intCur.execute("""select intercept_id,elnot,latitude,longitude,area_dist_maj,area_dist_min,orientation from """ + fq_idb + """intercepts
                           where intercept_id in (select intercept_id from comparison_geo_snapshot)
                           and substr(elnot,1,1) not in ('A','C','K','M','O')
                           and latitude is not null and longitude is not null and area_dist_maj is not null and area_dist_min is not null
                           and latitude >= 0 and longitude >= 0 and area_dist_maj > 0 and area_dist_min > 0 and orientation >= 0 
                           and orientation is not null and is_emitter = 'Y' and area_dist_maj <= %s
                           order by intercept_id""",(max_SMA,))
    for row in intCur:
        int_id = row[0]
        elnot = row[1]
        lat = float(row[2])
        lon = float(row[3])
        sma = float(row[4])
        smi = float(row[5])
        if row[6] == None:
            orient = 0.0
        else:
            orient = float(row[6])
    
   
        # Loop through sites and find best match, if any
        # for broader application constrain retrieved sites to realistic candidates (elnot, distance, etc.)
        siteCur.execute("select site_id,elnot,lat,lon,sma,smi,orient,heard_count,roa from site_geos where elnot = %s order by site_id",(elnot,))
        low_score = 9999
        best_site = 0
        for site in siteCur:
            site_id = site[0]
            selnot = site[1]
            slat = float(site[2])
            slon = float(site[3])
            ssma = float(site[4])
            ssmi = float(site[5])
            sorient = float(site[6])
            shc = int(site[7])
            sroa = site[8]
            
            if sroa == None:
                score = z_score(lat,lon,sma,smi,orient,slat,slon,ssma,ssmi,sorient,True)
            else:
                score = z_score(lat,lon,sma,smi,orient,slat,slon,float(sroa),float(sroa),sorient,True)
            
            #print("int_id,score",int_id,score)
            if score <= crit_val and score < low_score:
                best_site = site_id
                low_score = score
                rlat = lat
                rlon = lon
                rsma = sma
                rsmi = smi
                rorient = orient
                plat = slat
                plon = slon
                psma = ssma
                psmi = ssmi
                porient = sorient
                if use_sqrt == True:
                    phc = sqrt(min(shc,hc_max))  
                else:
                    phc = min(shc,hc_max)
        
        if best_site != 0:
            newlat,newlon,newsma,newsmi,neworient = kf_update(rlat,rlon,rsma,rsmi,rorient,plat,plon,psma,psmi,porient,phc,hc_weighting,w)
            newlat = round(newlat,3)
            newlon = round(newlon,3)
            newsma = round(newsma,3)
            newsmi = round(newsmi,3)
            neworient = round(neworient,1)
            siteCur.execute('''update site_geos set lat=%s,lon=%s,sma=%s,smi=%s,orient=%s where site_id=%s''',(newlat,newlon,newsma,newsmi,neworient,best_site,))
            siteCur.execute('''update site_geos set heard_count = heard_count+1 where site_id=%s''',(best_site,))
            siteCur.execute('''insert into site_ints(intercept_id,site_id) values(%s,%s)''',(int_id,best_site,))
            if track_deltas == True:
                siteCur.execute('''select max(iter) from site_history where site_id=%s''',(best_site,))
                prev_iter = siteCur.fetchone()[0]
                if prev_iter == None:
                    prev_iter = 0
                #print("prev_iter",prev_iter)
                siteCur.execute("insert into site_history(site_id,iter,lat,lon,sma,smi,orient) values(%s,%s,%s,%s,%s,%s,%s)",(best_site,prev_iter+1,newlat,newlon,newsma,newsmi,neworient,))
            siteCur.execute("Commit")
        else:
            if not save_tracks:   # populate residue
                if pop_residue:
                    siteCur.execute("insert into residue_ints(intercept_id,r_type) values(%s,'G')",(int_id,))
            else:   # store a new track
                siteCur.execute("select max(site_id) from site_geos")
                max_id = siteCur.fetchone()[0]
                if max_id == None:
                    max_id=0
                next_id=max_id+1
                # kosher to add sites inside a loop on sites? do subsequent ints have opportunity to associate to new sites or will gen redundant new sites?
                siteCur.execute("insert into site_geos(elnot,site_id,lat,lon,sma,smi,orient,heard_count) values(%s,%s,%s,%s,%s,%s,%s,%s)",(elnot,next_id,lat,lon,sma,smi,orient,1,))
                siteCur.execute("insert into site_ints(intercept_id,site_id) values(%s,%s)",(int_id,next_id,))
                if track_deltas == True:
                    siteCur.execute("insert into site_history(site_id,iter,lat,lon,sma,smi,orient) values(%s,%s,%s,%s,%s,%s,%s)",(next_id,0,lat,lon,sma,smi,orient,))
                
    siteCur.execute("commit")
    
    if not pop_residue:   # pre-processing residue so need to delete anything that associated
        siteCur.execute("delete from residue_ints where r_type = 'G' and intercept_id in (select intercept_id from site_ints)")
        siteCur.execute("commit")
    
    intCur.execute("delete from comparison_geo_snapshot")
    intCur.execute("commit")
        


    parmCur.close()
    intCur.close()
    siteCur.close()
    connection.close()
    end = datetime.now()
    #plot_all('sh')
    print("geo_assoc start: ",start)
    print("geo_assoc end: ",end)
    print("geo_assoc total:",end-start)


