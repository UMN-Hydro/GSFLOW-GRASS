class BuildINI(object):

    def __init__(self):
            
        # [paths]
        self.proj_name='gridTestCannon'
        self.gsflow_exe='/home/awickert/models/GSFLOW_1.2.0/bin/gsflow'
        self.gsflow_ver='1.2.0'
        self.gsflow_path_simdir='/home/awickert/GSFLOW2018/'
        
        # [GRASS_core]
        self.DEM_file_path_to_import=''
        self.gisdb='/PATH/TO/YOUR/HOME/DIRECTORY/PROBBALY/grassdata'
        self.version='74'

        # [run_mode]
        self.sw_1spinup_2restart='1'
        self.init_PRMSfil='/PATH/TO/RESTART/OUTFILE'
        self.init_MODfil='/PATH/TO/RESTART/OUTFILE.out'

        # [time]
        self.start_date='1938-05-12'
        self.end_date='1943-11-05'
        self.init_start_date=''

        # [GRASS_drainage]
        self.flow_weights=''
        self.threshold_drainage_area_meters2='36000000'
        self.MODFLOW_grid_resolution_meters='1000'
        self.outlet_point_x='532267.96306'
        self.outlet_point_y='4938984.83574'

        # [GRASS_hydraulics]
        self.icalc='1'
        self.channel_Mannings_n='0.035'
        self.channel_Mannings_n_grid=''
        self.channel_Mannings_n_vector=''
        self.channel_Mannings_n_vector_col=''
        self.overbank_Mannings_n='0.060'
        self.channel_width='5'
        self.channel_width_vector=''
        self.channel_width_vector_col=''
        self.floodplain_width='0'
        self.floodplain_width_vector=''
        self.floodplain_width_vector_col=''

        # [climate_inputs]
        self.fl_print_climate_hru='1'
        self.climate_data_file='/home/awickert/GSFLOW/Old2/CannonRiver/Precip_Zumbrota.txt'
        self.climate_hru_dir=''

        # [hydrogeologic_inputs]
        self.NLAY='2'
        self.DZ='20,200'
        self.fl_create_hydcond='0'
        self.hydcond='1E-5,1'
        self.finf='0.0002'

    def writeINI(self, filename):
        
        f = file(filename, 'w')
        
        def writeline(string):
            f.write(string + '\n')
        
        writeline("[paths]")
        writeline("proj_name="+self.proj_name)
        writeline("gsflow_exe="+self.gsflow_exe)
        writeline("gsflow_ver="+self.gsflow_ver)
        writeline("gsflow_path_simdir="+self.gsflow_path_simdir)
        writeline("")
        writeline("[GRASS_core]")
        writeline("DEM_file_path_to_import="+self.DEM_file_path_to_import)
        writeline("gisdb="+self.gisdb)
        writeline("version="+self.version)
        writeline("")
        writeline("[run_mode]")
        writeline("sw_1spinup_2restart="+self.sw_1spinup_2restart)
        writeline("init_PRMSfil="+self.init_PRMSfil)
        writeline("init_MODfil="+self.init_MODfil)
        writeline("")
        writeline("[time]")
        writeline("start_date="+self.start_date)
        writeline("end_date="+self.end_date)
        writeline("init_start_date="+self.init_start_date)
        writeline("")
        writeline("[GRASS_drainage]")
        writeline("flow_weights="+self.flow_weights)
        writeline("threshold_drainage_area_meters2="+self.threshold_drainage_area_meters2)
        writeline("MODFLOW_grid_resolution_meters="+self.MODFLOW_grid_resolution_meters)
        writeline("outlet_point_x="+self.outlet_point_x)
        writeline("outlet_point_y="+self.outlet_point_y)
        writeline("")
        writeline("[GRASS_hydraulics]")
        writeline("icalc="+self.icalc)
        writeline("channel_Mannings_n="+self.channel_Mannings_n)
        writeline("channel_Mannings_n_grid="+self.channel_Mannings_n_grid)
        writeline("channel_Mannings_n_vector="+self.channel_Mannings_n_vector)
        writeline("channel_Mannings_n_vector_col="+self.channel_Mannings_n_vector_col)
        writeline("overbank_Mannings_n="+self.overbank_Mannings_n)
        writeline("channel_width="+self.channel_width)
        writeline("channel_width_vector="+self.channel_width_vector)
        writeline("channel_width_vector_col="+self.channel_width_vector_col)
        writeline("floodplain_width="+self.floodplain_width)
        writeline("floodplain_width_vector="+self.floodplain_width_vector)
        writeline("floodplain_width_vector_col="+self.floodplain_width_vector_col)
        writeline("")
        writeline("[climate_inputs]")
        writeline("fl_print_climate_hru="+self.fl_print_climate_hru)
        writeline("climate_data_file="+self.climate_data_file)
        writeline("climate_hru_dir="+self.climate_hru_dir)
        writeline("")
        writeline("[hydrogeologic_inputs]")
        writeline("NLAY="+self.NLAY)
        writeline("DZ="+self.DZ)
        writeline("fl_create_hydcond="+self.fl_create_hydcond)
        writeline("hydcond="+self.hydcond)
        writeline("finf="+self.finf)

    
