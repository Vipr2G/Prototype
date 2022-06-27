import psycopg2
from VL2db_utils import db_con

try:
   connection = db_con()
   cur = connection.cursor()
   
   # #### Intercept related tables created via the "official" SQL scripts
   #cur.execute("drop table intercepts")
   #cur.execute("""Create table if not exists intercepts(
   #                                    intercept_id integer,
   #                                    mod_type text,
   #                                    elnot text,
   #                                    rf numeric,
   #                                    pd numeric,
   #                                    sp numeric,
   #                                    st text,
   #                                    sma numeric,
   #                                    smi numeric,
   #                                    orient numeric,
   #                                    lat numeric,
   #                                    lon numeric)""")
   #cur.execute("""Create table if not exists intercept_pris(
   #                                    intercept_id integer,
   #                                    pri_value numeric,
   #                                    pri_number integer)""")
   #
   
   ################################### PRIMARY IPE TABLES:  #####################################
   cur.execute("""Create table if not exists parameter_groups(
                                       group_id integer,
                                       elnot varchar(5),
                                       mod_type varchar(1),
                                       group_descr text,
                                       source text,
                                       status text)""")
   cur.execute("""Create table if not exists global_parm_clstrs(
                                       clstr_id integer,
                                       elnot varchar(5),
                                       mod_type varchar(1),
                                       parm varchar(3),
                                       parm_min numeric,
                                       parm_max numeric,
                                       bin_counts integer[])""")

   cur.execute("""Create table if not exists local_parm_clstrs(
                                       clstr_id integer,
                                       group_id integer,
                                       parm varchar(3),
                                       parm_min numeric,
                                       parm_max numeric)""")
   cur.execute("""Create table if not exists global_scans(
                                       elnot varchar(5),
                                       mod_type varchar(1),
                                       scan_type varchar(1))""")
   
   cur.execute("""Create table if not exists local_scans(
                                       group_id integer,
                                       scan_type varchar(1))""")

   
   cur.execute("""Create table if not exists sequences(
                                       seq_id integer,
                                       elnot varchar(5),
                                       mod_type varchar(1),
                                       group_id integer,
                                       heard_count integer,
                                       is_parent varchar(1),
                                       parent_id integer)""")
   cur.execute("""Create table if not exists sequence_elements(
                                       seq_id integer,
                                       clstr_id integer,
                                       position integer)""")
   cur.execute("""Create table if not exists sequence_ints(
                                       intercept_id integer,
                                       seq_id integer)""")
                                       
                                          
   ############################ Geo-related IPE tables:
   cur.execute("""Create table if not exists site_geos(
                                       site_id integer,
                                       elnot varchar(5),
                                       sma numeric,
                                       smi numeric,
                                       orient numeric,
                                       lat numeric,
                                       lon numeric,
                                       heard_count integer,
                                       ROA numeric)""")
   cur.execute("Create table if not exists site_ints(site_id integer,intercept_id integer)")



   ############################# Residue Management Tables:
    
   cur.execute("create table if not exists residue_ints(intercept_id integer, r_type text)")
   cur.execute("""Create table if not exists parm_residue(elnot varchar(5), 
                                                          mod_type varchar(1), 
                                                          parm varchar(3), 
                                                          idx integer,
                                                          idx_cnt integer)""")


   ########################### DATA FLOW MANAGEMENT:
   cur.execute("Create table if not exists comparison_queue(intercept_id integer)")
   cur.execute("""Create table if not exists comparison_snapshot(
                                       intercept_id integer)""")
   cur.execute("""Create table if not exists comparison_geo_snapshot(
                                       intercept_id integer)""")


   ########################### TUNING-RELATED TABLES:
   cur.execute("""Create table if not exists direct_tune_parms(
                                       elnot varchar(5),
                                       parm varchar(3),
                                       incr numeric,
                                       horizon numeric,
                                       p numeric,
                                       thresh numeric)""")

   
   
   ############################ ADMIN/LOGGING-RELATED TABLES:
   cur.execute("""Create table if not exists wob_event_log(
                                       event_date timestamptz,
                                       event_txt text)""")



   ########################## OPTIONAL UTILITY SUPPORT TABLES:
   cur.execute("Create table if not exists site_ints2(site_id integer,intercept_id integer)")

   cur.execute("""Create table if not exists site_history(
                                       site_id integer,
                                       iter integer,
                                       sma numeric,
                                       smi numeric,
                                       orient numeric,
                                       lat numeric,
                                       lon numeric)""")
   
   cur.execute("""Create table if not exists geo_tuning_parms(
                                       elnot varchar(5),
                                       roa numeric,
                                       grid_res numeric)""")


   cur.execute("Commit")

   cur.close()
   connection.close()

except (Exception, psycopg2.Error) as error :
    print ("Error while fetching data from PostgreSQL", error)

