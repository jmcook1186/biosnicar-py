 
% %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% The Snow, Ice, and Aerosol Radiative Model with Adding-Doubling
% solver, version 4.0 (SNICAR-ADv4)
%
% References cited and open-source license are at end of file.
%
% This code calculates snow and ice spectral albedo and radiative fluxes
% through multi-layer snowpacks using a Delta-Eddington
% Adding-Doubling two-stream approximation (Briegleb and Light, 2007;
% Dang et al., 2019).
%
% The code reads particle optical properties from an external library
% of NetCDF files, based on user-defined properties of each snow
% layer.  The root directory containing this collection of files must
% be specified in variable "dir_op_root".
%
% Input parameters are described below. The user can run this code as
% a script instead of a function by commenting out the function line
% and setting '1==1' below the function definition.
%
% Please contact cwhicker@umich.edu with questions or concerns relating to
% SNICAR-ADv4
%
%%%%%%%%%   FEATURES IN THIS RELEASE   %%%%%%%%%%%%
%  - Pure ice properties (Whicker et al., 2021)
%  - Fresnel Layer and snells law applied to ice surfaces (Brigleb and
%    Light, 2007) and extended to be spectrally varying (Whicker et al.,
%    2021)
%  - Glacier algae optical properties (Whicker et al., 2021; Cook et al.,
%    2020) 
%  - Adding-Doubling solver (Dang et al., 2019)
%  - Snow algae (Cook et al., 2017)
%  - Glacier Algae (Cook et al., 2020)
%  - Non-spherical ice particles (He et al., 2017)
%  - CO2 ice (Hansen, 2005; Singh and Flanner, 2016)
%  - Additional dust optical properties:
%     - Saharan dust (Balkanski et al., 2007)
%     - San Juan Mountains, Colorado (Skiles et al., 2017)
%     - Greenland (Polashenski et al., 2015)
%     - Martian dust (Wolff et al., 2009;2010; Singh and Flanner, 2016)
%  - Brown carbon (Kirchstetter et al., 2004)
%  - Volcanic ash (Flanner et al., 2014)
%  - A larger size bin (5-50um radius) for dust and ash particles
%  - Multiple options for H2O ice refractive indices
%    - Warren (1984)
%    - Warren and Brandt (2008)
%    - Picard et al. (2016)
%  - Several new options for surface spectral irradiance weighting
%  - Extension of the simulated spectrum down to 0.2um (Flanner et al.,
%    2021)


%%%%%%%%%%  Input parameters: %%%%%%%%%%%
%
% direct_beam:   Direct or diffuse incident radiation
%                  1 = direct-beam (clear-sky)
%                  0 = diffuse (cloudy)
%
% atm:           Atmospheric profile used to obtain surface-incident spectral flux distribution
%                and subsequent broadband albedo:
%                  1 = mid-latitude winter
%                  2 = mid-latitude summer
%                  3 = sub-Arctic winter
%                  4 = sub-Arctic summer
%                  5 = Summit,Greenland (sub-Arctic summer, surface pressure of 796hPa)
%                  6 = High Mountain (summer, surface pressure of 556 hPa)
%                  7 = Top-of-atmosphere
%
% flx_dwn_bb:    Broadband surface-incident solar flux [W/m2]. 
%                  Spectral flux fractions, determined by the
%                  atmospheric profile and direct_beam/cloudy flag, 
%                  are scaled by this constant

% coszen:        Cosine of solar zenith angle (only applies when direct_beam=1)
%                  Range: 0-1
%
% ice_ri:        Flag for ice refractive index data to use
%                  1 = Warren (1984)
%                  2 = Warren and Brandt (2008)
%                  3 = Picard et al (2016)
%                  4 = CO2 ice, Hansen (2005)
%
% R_sfc_all_wvl  Broadband albedo of underlying surface 
%                  Range: 0-1
%                  (user can also define a spectrally-dependent albedo below)
%
% dz:            Array of snow layer thicknesses [meters]. Top layer has index 1. 
%                  *** THE LENGTH OF THIS ARRAY DEFINES THE NUMBER OF SNOW LAYERS ***
%
% rho_snw:       Array of snow layer densities [kg m-3].
%                  Must have same length as dz
%
% lyr_typ:       Array of layer type indicators
%                   Must have the same length as dz
%                        1 = snow 
%                        2 = ice
% 
% rds_snw:       Array of snow layer effective grain radii [microns]
%                  Must have same length as dz
%
% options for asperical ice particles added by Cenlin He:
% sno_shp:      Snow grain shape
%               1=sphere (original SNICAR scheme)
%               2=spheroid; 3=hexagonal plate; 4=koch snowflake
%               air bubbles within ice are always spheres
% sno_fs:       Shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
%               0=use recommended default value (He et al., 2017);
%               others(0<fs<1)= use user-specified value
%               only activated when sno_shp > 1 (i.e. nonspherical)
% sno_ar:       Aspect ratio: ratio of grain width to length
%               0=use recommended default value (He et al., 2017);
%               others(0.1<fs<20)= use user-specified value
%               only activated when sno_shp > 1 (i.e. nonspherical)
%
% mss_cnc_sot1:  Layer-specific array of mass mixing ratio of black carbon species 1 (uncoated BC)
%                  Units of parts per billion, ng g-1
%                  NOTE: This is mass of impurity per mass of snow
%                  (i.e., mass of impurity / mass of ice+impurity)
%
% mss_cnc_sot2:  Layer-specific array of mass mixing ratio of black carbon species 2 (sulfate-coated)
%                  Units of parts per billion, ng g-1
%
% mss_cnc_brc1:  Layer-specific array of mass mixing ratio of brown carbon species 1
%                  Units of parts per billion, ng g-1
%
% mss_cnc_brc2:  Layer-specific array of mass mixing ratio of brown carbon species 2 (sulfate-coated)
%                  Units of parts per billion, ng g-1
%
% dust_type:     Type of dust to use:
%                  1 = Saharan dust (Balkanski et al., 2007, central hematite)
%                  2 = San Juan Mountains, CO (Skiles et al, 2017)
%                  3 = Greenland (Polashenski et al., 2015, central absorptivity)
%                  4 = Martian dust (Wolff et al., 2009;2010, Singh and Flanner, 2016)
%
% mss_cnc_dst1:  Layer-specific array of mass mixing ratio of dust species 1 (radii of 0.05-0.5um)
%                  Units of parts per million, ug g-1
%
% mss_cnc_dst2:  Layer-specific array of mass mixing ratio of dust species 2 (radii of 0.5-1.25um)
%                  Units of parts per million, ug g-1
%
% mss_cnc_dst3:  Layer-specific array of mass mixing ratio of dust species 3 (radii of 1.25-2.5um)
%                  Units of parts per million, ug g-1
%
% mss_cnc_dst4:  Layer-specific array of mass mixing ratio of dust species 4 (radii of 2.5-5.0um)
%                  Units of parts per million, ug g-1
%
% mss_cnc_dst5:  Layer-specific array of mass mixing ratio of dust species 5 (radii of 5.0-50um)
%                  Units of parts per million, ug g-1
%
% ash_type:      Type of volcanic ash to use
%                  1 = Eyjafjallajokull ash (Flanner et al.,2014)
%
% mss_cnc_ash1:  Layer-specific array of mass mixing ratio of volcanic ash species 1 (radii of 0.05-0.5um)
%                  Units of parts per million, ug g-1
%
% mss_cnc_ash2:  Layer-specific array of mass mixing ratio of volcanic ash species 2 (radii of 0.5-1.25um)
%                  Units of parts per million, ug g-1
%
% mss_cnc_ash3:  Layer-specific array of mass mixing ratio of volcanic ash species 3 (radii of 1.25-2.5um)
%                  Units of parts per million, ug g-1
%
% mss_cnc_ash4:  Layer-specific array of mass mixing ratio of volcanic ash species 4 (radii of 2.5-5.0um)
%                  Units of parts per million, ug g-1
%
% mss_cnc_ash5:  Layer-specific array of mass mixing ratio of volcanic ash species 5 (radii of 5.0-50um)
%                  Units of parts per million, ug g-1
%
% nw_alg_cell_nbr_conc: Layer-specific array of snow algal cell abundance
%                       Units of cells/mL or cells/gH20
%
% alg_rds:       Layer-specific mean cell radius of Gaussian distribution [um]
%                  Valid values are: [1 2 5 10 15 20 25 30 40 50] um
%
% dcmf_pig_chla: Dry cell mass fraction of chlorophyll-a
%                  Valid values are: [0.000 0.005 0.010 0.015 0.020 0.025 0.030]
%
% dcmf_pig_chlb: Dry cell mass fraction of chlorophyll-b
%                  Valid values are: [0.000 0.005 0.010 0.015 0.020 0.025 0.030]
%
% dcmf_pig_cara: Dry cell mass fraction of photoprotective carotenoids
%                  Valid values are: 0.00-0.15 by 0.01 increments
%
% dcmf_pig_carb: Dry cell mass fraction of photosynthetic carotenoids  
%                  Valid values are: 0.00-0.15 by 0.01 increments
%
% glc_alg_mss_cnc: Dry cell mass of glacier algae [ng of dry cell mass/g of ice]
%
% input_args.glc_alg_rds: Radius of the glacier algae [um]
%                            Valid values are: 4
%
% input_args.glc_alg_len: Length of the glacier algae [um]
%                            Valid values are: 40
%
%%%%%  Output data: %%%%%
%
% All output data is contained in the structure "data_out",
% described at end of this file


%%%%%%%%%%%%% BEGIN CODE: %%%%%%%%%%%%%%%

function data_out = snicarAD_v4(input_args)

% ROOT DIRECTORY FOR ALL OPTICAL AND SPECTRAL IRRADIANCE FILES
dir_op_root = '/data/jcook/SNICAR_ADv4/snicar_480band/';


% IF USER PREFERS TO RUN THIS AS A SCRIPT INSTEAD OF A FUNCTION:
%   (1) COMMENT OUT THE FUNCTION CALL ABOVE
%   (2) SET 1==1 BELOW AND DEFINE ALL INPUT IN THIS BLOCK

if (0==1)

    % RADIATIVE TRANSFER CONFIGURATION:
    direct_beam   = 1;   % 1= Direct-beam incident flux, 0= Diffuse incident flux
                         % NOTE that cloudy-sky spectral fluxes are loaded when direct_beam=0
    
    % ATMOSPHERIC PROFILE for surface-incident flux:
    %    1 = mid-latitude winter
    %    2 = mid-latitude summer
    %    3 = sub-Arctic winter
    %    4 = sub-Arctic summer
    %    5 = Summit,Greenland (sub-Arctic summer, surface pressure of 796hPa)
    %    6 = High Mountain (summer, surface pressure of 556 hPa)
    %    7 = Top-of-atmosphere
    % NOTE that clear-sky spectral fluxes are loaded when direct_beam=1,
    % and cloudy-sky spectral fluxes are loaded when direct_beam=0
    atm = 1;

    % Broadband surface-incident solar flux [W/m2]:
    %  (used to scale spectral flux fractions)
    flx_dwn_bb = 1.0;

    % COSINE OF SOLAR ZENITH ANGLE FOR DIRECT-BEAM RT
    coszen = 0.5;
  
    % ICE REFRACTIVE INDEX DATASET TO USE:
    ice_ri = 3;
    
    % REFLECTANCE OF SURFACE UNDERLYING SNOW:
    % (value applied to all wavelengths.  user can also specify
    % spectrally-dependent ground albedo below)
    R_sfc_all_wvl = 0.25;

    % SNOW OR ICE LAYER THICKNESSES [m]:
    %dz = [0.02 0.02 0.05];
    dz = [1000];
    nbr_lyr = length(dz);  % number of layers
    
    % LAYER MEDIUM TYPE [ 1=snow, 2=ice]
    %  Must have same length as dz
    lyr_typ(1:nbr_lyr) = [1];
    
    % SNOW DENSITY FOR EACH LAYER (units: kg/m3)
    rho_snw(1:nbr_lyr) = 150;

    % SNOW GRAIN SIZE FOR EACH SNOW LAYER ?R BUBBLE RADIUS FOR ICE
    % (units: microns):
    rds_snw(1:nbr_lyr) = 1000;
  
    % Options added by Cenlin He for nonspherical ice particles based on
    % the parameterizations described by He et al. (2017,
    % doi:10.1175/JCLI-D-17-0300.1)
    % these options are only relevant to snow grains (ie lyr_typ = 1) 
    
    sno_shp(1:nbr_lyr)  = 2;    % Snow grain shape option
                                % 1=sphere; 2=spheroid; 3=hexagonal plate; 4=koch snowflake

    sno_fs(1:nbr_lyr)   = 0;    % Shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
                                % 0=use recommended default value (He et al. 2017);
                                % others(0<fs<1)= use user-specified value
                                % only activated when sno_shp > 1 (i.e. nonspherical)
                                
    sno_ar(1:nbr_lyr)   = 0;    % Aspect ratio: ratio of grain width to length
                                % 0=use recommended default value (He et al. 2017);
                                % others(0.1<fs<20)= use user-specified value
                                % only activated when sno_shp > 1 (i.e. nonspherical)
    
    % type of dust:
    dust_type = 1;              % 1=Saharan, 2=Colorado, 3=Greenland, 4=Mars
    
    % type of volcanic ash:
    ash_type = 1;               % 1 = Eyjafjallajokull
    
    % PARTICLE MASS MIXING RATIOS (units: ng g-1)
    % NOTE: This is mass of impurity per mass of snow
    %  (i.e., mass of impurity / mass of ice+impurity)
    mss_cnc_sot1(1:nbr_lyr)  = 0.0;    % uncoated black carbon [ng/g]
    mss_cnc_sot2(1:nbr_lyr)  = 0.0;    % coated black carbon [ng/g]
    mss_cnc_brc1(1:nbr_lyr)  = 0.0;    % uncoated brown carbon [ng/g]
    mss_cnc_brc2(1:nbr_lyr)  = 0.0;    % coated brown carbon [ng/g]
    mss_cnc_dst1(1:nbr_lyr)  = 0.0;    % dust species 1 [ug/g]
    mss_cnc_dst2(1:nbr_lyr)  = 0.0;    % dust species 2 [ug/g]
    mss_cnc_dst3(1:nbr_lyr)  = 0.0;    % dust species 3 [ug/g]
    mss_cnc_dst4(1:nbr_lyr)  = 0.0;    % dust species 4 [ug/g]
    mss_cnc_dst5(1:nbr_lyr)  = 0.0;    % dust species 5 [ug/g]
    mss_cnc_ash1(1:nbr_lyr)  = 0.0;    % volcanic ash species 1 [ug/g]
    mss_cnc_ash2(1:nbr_lyr)  = 0.0;    % volcanic ash species 2 [ug/g]
    mss_cnc_ash3(1:nbr_lyr)  = 0.0;    % volcanic ash species 3 [ug/g]
    mss_cnc_ash4(1:nbr_lyr)  = 0.0;    % volcanic ash species 4 [ug/g]
    mss_cnc_ash5(1:nbr_lyr)  = 0.0;    % volcanic ash species 5 [ug/g]

    snw_alg_cell_nbr_conc(1:nbr_lyr) = 0.0;    % algae [cells/mL]
    alg_rds(1:nbr_lyr)               = 10;     % mean cell radius (um)
    dcmf_pig_chla(1:nbr_lyr)         = 0.02;   % dry cell mass fraction of chlorophyll-a
    dcmf_pig_chlb(1:nbr_lyr)         = 0.02;   % dry cell mass fraction of chlorophyll-b
    dcmf_pig_cara(1:nbr_lyr)         = 0.05;   % dry cell mass fraction of photoprotective_carotenoids
    dcmf_pig_carb(1:nbr_lyr)         = 0.00;   % dry cell mass fraction of photosynthetic_carotenoids  
    
    glc_alg_mss_cnc(1:nbr_lyr)       = [0];    % GLACIER algae [UNITS ng/g] 
    glc_alg_rds(1)                   = [4];    % GLACIER algae radius [um]
    glc_alg_len(1)                   = [40];   % GLACIER algae length [um]
    
else
    % USE THE DRIVER INPUT ARGUMENTS
    direct_beam   = input_args.direct_beam;     % direct beam or diffuse
    atm           = input_args.atm;             % atmospheric profile
    flx_dwn_bb    = input_args.flx_dwn_bb;      % broadband surface insolation
    coszen        = input_args.coszen;          % cosine of solar zenith angle
    ice_ri        = input_args.ice_ri;          % ice refractive index data
    R_sfc_all_wvl = input_args.R_sfc_all_wvl;   % albedo of underlying surface
    dz            = input_args.dz;              % layer thicknesses
    lyr_typ       = input_args.lyr_typ;         % snow or ice layer indicator [1 = snow, 2=ice]
    rho_snw       = input_args.rho_snw;         % layer densities
    rds_snw       = input_args.rds_snw;         % layer grain radii
    sno_shp       = input_args.sno_shp;         % snow layer grain shape
    sno_fs        = input_args.sno_fs;          % snow layer grain shape factor
    sno_ar        = input_args.sno_ar;          % snow layer grain aspect ratio
    dust_type     = input_args.dust_type;       % type of dust
    ash_type      = input_args.ash_type;        % type of ash
    mss_cnc_sot1  = input_args.mss_cnc_sot1;    % uncoated black carbon [ng/g]
    mss_cnc_sot2  = input_args.mss_cnc_sot2;    % coated black carbon [ng/g]
    mss_cnc_brc1  = input_args.mss_cnc_brc1;    % uncoated brown carbon [ng/g]
    mss_cnc_brc2  = input_args.mss_cnc_brc2;    % coated brown carbon [ng/g]
    mss_cnc_dst1  = input_args.mss_cnc_dst1;    % dust species 1 [ng/g]
    mss_cnc_dst2  = input_args.mss_cnc_dst2;    % dust species 2 [ng/g]
    mss_cnc_dst3  = input_args.mss_cnc_dst3;    % dust species 3 [ng/g]
    mss_cnc_dst4  = input_args.mss_cnc_dst4;    % dust species 4 [ng/g]
    mss_cnc_dst5  = input_args.mss_cnc_dst5;    % dust species 5 [ng/g]
    mss_cnc_ash1  = input_args.mss_cnc_ash1;    % volcanic ash species 1 [ng/g]
    mss_cnc_ash2  = input_args.mss_cnc_ash2;    % volcanic ash species 2 [ng/g]
    mss_cnc_ash3  = input_args.mss_cnc_ash3;    % volcanic ash species 3 [ng/g]
    mss_cnc_ash4  = input_args.mss_cnc_ash4;    % volcanic ash species 4 [ng/g]
    mss_cnc_ash5  = input_args.mss_cnc_ash5;    % volcanic ash species 5 [ng/g]
    
    snw_alg_cell_nbr_conc = input_args.snw_alg_cell_nbr_conc;   % algae [cells/mL]
    alg_rds               = input_args.alg_rds;         % mean cell radius (um)
    dcmf_pig_chla         = input_args.dcmf_pig_chla;   % dry cell mass fraction of chlorophyll-a
    dcmf_pig_chlb         = input_args.dcmf_pig_chlb;   % dry cell mass fraction of chlorophyll-b
    dcmf_pig_cara         = input_args.dcmf_pig_cara;   % dry cell mass fraction of photoprotective carotenoids
    dcmf_pig_carb         = input_args.dcmf_pig_carb;   % dry cell mass fraction of photosynthetic carotenoids
    
    glc_alg_mss_cnc       = input_args.glc_alg_mss_cnc;  % GLACIER algae [UNITS ng/g] 
    glc_alg_rds           = input_args.glc_alg_rds;      % GLACIER algae radius [um]
    glc_alg_len           = input_args.glc_alg_len;      % GLACIER algae length [um]

    nbr_lyr       = length(dz);  % number of snow layers
end;

% density of pure ice [kg/m3]
rho_ice =  917; 

% density of air [kg / m3]
rho_air = 1.025; 

% Sub-directories of NetCDF files for 
% (1) optical properties of light-absorbing impurities, 
% (2) optical properties of snow algae,
% (3) surface spectral irradiance profiles:
dir_lai     = strcat(dir_op_root,'lai/');
dir_alg     = strcat(dir_op_root,'alg_pig/');
dir_spc     = strcat(dir_op_root,'fsds/');

% Set wavelength grid (um) and wavelength number:
wvl     = [0.205:0.01:4.995];
nbr_wvl = length(wvl);

% file substrings for ice Mie parameters:
if (ice_ri == 1)
    stb1 = 'ice_Wrn84';
    ice_rfidx_re_str = strcat('re_',stb1(5:end));
    ice_rfidx_im_str = strcat('im_',stb1(5:end));
elseif (ice_ri == 2)
    stb1 = 'ice_Wrn08';
    ice_rfidx_re_str = strcat('re_',stb1(5:end));
    ice_rfidx_im_str = strcat('im_',stb1(5:end));
elseif (ice_ri == 3)
    stb1 = 'ice_Pic16';
    ice_rfidx_re_str = strcat('re_',stb1(5:end));
    ice_rfidx_im_str = strcat('im_',stb1(5:end));
elseif (ice_ri == 4)
    stb1 = 'co2ice';
    ice_rfidx_re_str = strcat('re_',stb1);
    ice_rfidx_im_str = strcat('im_',stb1);
end;

% find the first ice layer 
%     occurs between the last snow layer and the first ice layer
%     if the top layer is ice total refelction will occur @ high SZAs
%     we recommend including a SSL to avoid aphysical results 
kfrsnl = find(lyr_typ==2, 1 );
if isempty(kfrsnl) == 1
                  kfrsnl=0;
else
    % read in precalcualted FL diffuse reflection  
    FL_r_dif_a = ncread(strcat(dir_op_root,'FL_reflection_diffuse.nc'),strcat('R_dif_fa_',stb1));
    FL_r_dif_b  = ncread(strcat(dir_op_root,'FL_reflection_diffuse.nc'),strcat('R_dif_fb_',stb1));
end

% subdirectory for ice optical properties:
dir_ice  = strcat(dir_op_root,stb1,'/');

% subdirectory for ice bubble optical properties:
stb2 = 'bbl'; % for log normal bubble size dist
dir_bbl  = strcat(dir_op_root,stb2,'/');

% read file with ice refractive index 
ice_rfidx_file = strcat(dir_op_root,"rfidx_ice.nc");

% real and imaginary ice refractive index 
rfidx_ice_re = ncread(ice_rfidx_file,ice_rfidx_re_str)';
rfidx_ice_im = ncread(ice_rfidx_file,ice_rfidx_im_str)';

% adjusted index of refraction for ice (Liou 2004 Eq. 5.4.18)
temp1 = rfidx_ice_re.^2 - rfidx_ice_im.^2 +sin(acos(coszen)).^2;
temp2 = rfidx_ice_re.^2 - rfidx_ice_im.^2 -sin(acos(coszen)).^2;
Nreal = (sqrt(2)/2) .* ( temp1 + (temp2.^2 + 4*rfidx_ice_re.^2.*rfidx_ice_im.^2).^(0.5) ).^0.5;

% filenames for impurity optical properties:
fl_sot1  = 'bc_ChCB_rn40_dns1270.nc';
fl_sot2  = 'bc_ChCB_rn40_dns1270_slfcot.nc';

fl_brc1  = 'brC_Kirch_BCsd.nc';
fl_brc2  = 'brC_Kirch_BCsd_slfcot.nc';

fl_glc_alg = strcat('Cook2020_glacier_algae_', num2str(glc_alg_rds),'_', num2str(glc_alg_len),'.nc');

if (dust_type==1)
    fl_dst1  = 'dust_balkanski_central_size1.nc';
    fl_dst2  = 'dust_balkanski_central_size2.nc';
    fl_dst3  = 'dust_balkanski_central_size3.nc';
    fl_dst4  = 'dust_balkanski_central_size4.nc';
    fl_dst5  = 'dust_balkanski_central_size5.nc';
elseif (dust_type==2)
    fl_dst1  = 'dust_skiles_size1.nc';
    fl_dst2  = 'dust_skiles_size2.nc';
    fl_dst3  = 'dust_skiles_size3.nc';
    fl_dst4  = 'dust_skiles_size4.nc';
    fl_dst5  = 'dust_skiles_size5.nc';
elseif (dust_type==3)
    fl_dst1  = 'dust_greenland_central_size1.nc';
    fl_dst2  = 'dust_greenland_central_size2.nc';
    fl_dst3  = 'dust_greenland_central_size3.nc';
    fl_dst4  = 'dust_greenland_central_size4.nc';
    fl_dst5  = 'dust_greenland_central_size5.nc';
elseif (dust_type==4)
    fl_dst1  = 'dust_mars_size1.nc';
    fl_dst2  = 'dust_mars_size2.nc';
    fl_dst3  = 'dust_mars_size3.nc';
    fl_dst4  = 'dust_mars_size4.nc';
    fl_dst5  = 'dust_mars_size5.nc';
end;

if (ash_type==1)
    fl_ash1  = 'volc_ash_eyja_central_size1.nc';
    fl_ash2  = 'volc_ash_eyja_central_size2.nc';
    fl_ash3  = 'volc_ash_eyja_central_size3.nc';
    fl_ash4  = 'volc_ash_eyja_central_size4.nc';
    fl_ash5  = 'volc_ash_eyja_central_size5.nc';
end;

% create cell structure of impurity file names
f1  = strcat(dir_lai,fl_sot1);
f2  = strcat(dir_lai,fl_sot2);
f3  = strcat(dir_lai,fl_brc1);
f4  = strcat(dir_lai,fl_brc2);
f5  = strcat(dir_lai,fl_dst1);
f6  = strcat(dir_lai,fl_dst2);
f7  = strcat(dir_lai,fl_dst3);
f8  = strcat(dir_lai,fl_dst4);
f9  = strcat(dir_lai,fl_dst5);
f10 = strcat(dir_lai,fl_ash1);
f11 = strcat(dir_lai,fl_ash2);
f12 = strcat(dir_lai,fl_ash3);
f13 = strcat(dir_lai,fl_ash4);
f14 = strcat(dir_lai,fl_ash5);
f15 = strcat(dir_alg,fl_glc_alg);

tmp_char = strvcat(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15);
files    = cellstr(tmp_char);

% NUMBER OF PARTICLE SPECIES IN SNOW (ICE AND ALGAE EXCLUDED)
tmp_sz  = size(files);
nbr_aer = tmp_sz(1);

% REFLECTANCE OF UNDERLYING SURFACE
% (OPTIONAL: User can set spectrally-dependent surface albedo here):
R_sfc(1:nbr_wvl,1) = R_sfc_all_wvl;

% files containing surface-incident spectral flux fractions:
stb_spc1 = 'swnb_480bnd_';
var_spc  = 'flx_frc_sfc';

if (atm==1)
    stb_atm = 'mlw';
elseif (atm==2)
    stb_atm = 'mls';
elseif (atm==3)
    stb_atm = 'saw';
elseif (atm==4)
    stb_atm = 'sas';
elseif (atm==5)
    stb_atm = 'smm';
elseif (atm==6)
    stb_atm = 'hmn';
elseif (atm==7)
    stb_atm = 'toa';
    var_spc = 'flx_frc_toa';
end;

mu_not  = coszen;
slr_zen = round(acosd(coszen)); % integer value of solar zenith angle
if (slr_zen>89)
    % surface irradiance profiles exist for theta = [0 1 ... 85]
    slr_zen=89;
end;

% Set incident solar flux spectral distribution:
if (direct_beam == 1)
 
    % Option 1: Use zenith-angle-independent (60 degree) surface irradiance profile:
    %fi_spc = strcat(dir_spc,stb_spc1,stb_atm,'_clr.nc');

    % Option 2: Use zenith-angle-dependent surface irradiance profile:
    if (atm < 7)
        fi_spc  = strcat(dir_spc,stb_spc1,stb_atm,'_clr_SZA',sprintf('%02d',slr_zen),'.nc');
    else
        % No zenith-angle dependence for TOA irradiance
        fi_spc  = strcat(dir_spc,stb_spc1,'toa.nc');
    end;
    
    flx_slr                   = ncread(fi_spc,var_spc);
    flx_slr(find(flx_slr==0)) = 1E-30;
    
    % direct-beam incident spectral flux [W/m2/band]
    Fs(1:nbr_wvl,1) = flx_dwn_bb.*flx_slr./(mu_not*pi);
    
    % diffuse incident spectral flux [W/m2/band]:
    Fd(1:nbr_wvl,1) = 0;
    
elseif (direct_beam == 0)
    if (atm < 7)
        fi_spc  = strcat(dir_spc,stb_spc1,stb_atm,'_cld.nc');
    else
        % TOA profile has no clear/cloud distinction
        fi_spc  = strcat(dir_spc,stb_spc1,'toa.nc');
    end;

    flx_slr                   = ncread(fi_spc,var_spc);
    flx_slr(find(flx_slr==0)) = 1E-30;
    
    % direct-beam incident spectral flux [W/m2/band]
    Fd(1:nbr_wvl,1) = flx_dwn_bb.*flx_slr;
    
    % diffuse incident spectral flux [W/m2/band]
    Fs(1:nbr_wvl,1) = 0;
end;

% Set indices for visible (0.2-0.7um) and near-IR (0.7-5.0um) bands
vis_max_idx = 50;
nir_max_idx = length(wvl);


%%% Constants for aspherical ice particles %%%
% g_snw asymmetry factor parameterization coefficients (6 bands) from
% Table 3 & Eqs. 6-7 in He et al. (2017)
% assume same values for 4-5 um band, which leads to very small biases (<3%)
g_wvl = [0.25,0.70,1.41,1.90,2.50,3.50,4.00,5.00]; % wavelength (um) division point
g_wvl_center = g_wvl(2:8)/2 + g_wvl(1:7)/2 ; % center point for wavelength band
g_b0 = [9.76029E-01,9.67798E-01,1.00111E+00,1.00224E+00,9.64295E-01,9.97475E-01,9.97475E-01];
g_b1 = [5.21042E-01,4.96181E-01,1.83711E-01,1.37082E-01,5.50598E-02,8.48743E-02,8.48743E-02];
g_b2 = [-2.66792E-04,1.14088E-03,2.37011E-04,-2.35905E-04,8.40449E-04,-4.71484E-04,-4.71484E-04];
% Tables 1 & 2 and Eqs. 3.1-3.4 from Fu, 2007
g_F07_c2 = [1.349959e-1,1.115697e-1,9.853958e-2,5.557793e-2,-1.233493e-1,0.0,0.0];
g_F07_c1 = [-3.987320e-1,-3.723287e-1,-3.924784e-1,-3.259404e-1,4.429054e-2,-1.726586e-1,-1.726586e-1];
g_F07_c0 = [7.938904e-1,8.030084e-1,8.513932e-1,8.692241e-1,7.085850e-1,6.412701e-1,6.412701e-1];
g_F07_p2 = [3.165543e-3,2.014810e-3,1.780838e-3,6.987734e-4,-1.882932e-2,-2.277872e-2,-2.277872e-2];
g_F07_p1 = [1.140557e-1,1.143152e-1,1.143814e-1,1.071238e-1,1.353873e-1,1.914431e-1,1.914431e-1];
g_F07_p0 = [5.292852e-1,5.425909e-1,5.601598e-1,6.023407e-1,6.473899e-1,4.634944e-1,4.634944e-1];


%%% Constants for snow algae %%%
sigma_alg = alg_rds.*0.1;  % standard deviation of Gaussian size distribution [um] (consistent with Mie)
rho_alg   = 1080;          % algae density [kg/m3] (consistent with Mie calculations)


% Establish layer-specific optical properties for ice and algae:
for n=1:nbr_lyr
    
    if lyr_typ(n) == 1 % SNOW
        
        % ice grains (snow)
        fl_ice = strcat(dir_ice,stb1,'_',sprintf('%04d',rds_snw(n)),'.nc');
        
        % single-scatter albedo, mass extinction coefficient, and
        % asymmetry paramater of ice:
        omega_ice(:,n)       = ncread(fl_ice,'ss_alb');
        ext_cff_mss_ice(:,n) = ncread(fl_ice,'ext_cff_mss');
        
        %%% shape-dependent asymmetry factors (Cenlin He) %%%%
        if (sno_shp(n) == 1)  % spheres from Mie properties (original SNICAR scheme)
            g_ice(:,n)       = ncread(fl_ice,'asm_prm'); % asymmetry paramater
            
        elseif (sno_shp(n) == 2) % 2=spheroid, He et al. (2017) parameterization
            diam_ice = 2.0 .* rds_snw(n); % effective snow grain diameter
            if (sno_fs(n) == 0)
                fs_sphd = 0.929; % default shape factor for spheroid; He et al. (2017), Table 1
            else
                fs_sphd = sno_fs(n);
            end
            fs_hex = 0.788; % shape factor for hexagonal plate (reference);
            if (sno_ar(n) == 0)
                AR_tmp = 0.5; % default aspect ratio for spheroid; He et al. (2017), Table 1
            else
                AR_tmp = sno_ar(n);
            end
            g_ice_Cg_tmp = g_b0 .* (fs_sphd/fs_hex).^g_b1 .* diam_ice.^g_b2; % Eq.7, He et al. (2017)
            gg_ice_F07_tmp = g_F07_c0 + g_F07_c1 .* AR_tmp + g_F07_c2 .* AR_tmp^2;% Eqn. 3.1 in Fu (2007)
            
        elseif (sno_shp(n) == 3) % 3=hexagonal plate, He et al. 2017 parameterization
            diam_ice = 2.0 .* rds_snw(n); % effective snow grain diameter
            if (sno_fs(n) == 0)
                fs_hex0 = 0.788; % default shape factor for hexagonal plates; He et al. (2017), Table 1
            else
                fs_hex0 = sno_fs(n);
            end
            fs_hex = 0.788; % shape factor for hexagonal plate (reference);
            if (sno_ar(n) == 0)
                AR_tmp = 2.5; % default aspect ratio for hexagonal plate; He et al. (2017), Table 1
            else
                AR_tmp = sno_ar(n);
            end
            g_ice_Cg_tmp = g_b0 .* (fs_hex0/fs_hex).^g_b1 .* diam_ice.^g_b2; % Eq.7, He et al. (2017)
            gg_ice_F07_tmp = g_F07_p0 + g_F07_p1 .* log(AR_tmp) + g_F07_p2 .* (log(AR_tmp))^2;   % Eqn. 3.3 in Fu (2007)
            
        elseif (sno_shp(n) == 4) % 4=koch snowflake, He et al. (2017) parameterization
            diam_ice = 2.0 .* rds_snw(n) ./ 0.544; % effective snow grain diameter
            if (sno_fs(n) == 0)
                fs_koch = 0.712; % default shape factor for koch snowflake; He et al. (2017), Table 1
            else
                fs_koch = sno_fs(n);
            end
            fs_hex = 0.788; % shape factor for hexagonal plate (reference);
            if (sno_ar(n) == 0)
                AR_tmp = 2.5; % default aspect ratio for koch snowflake; He et al. (2017), Table 1
            else
                AR_tmp = sno_ar(n);
            end
            g_ice_Cg_tmp = g_b0 .* (fs_koch/fs_hex).^g_b1 .* diam_ice.^g_b2; % Eq.7, He et al. (2017)
            gg_ice_F07_tmp = g_F07_p0 + g_F07_p1 .* log(AR_tmp) + g_F07_p2 .* (log(AR_tmp))^2;   % Eqn. 3.3 in Fu (2007)
            
        end;
        
        if (sno_shp(n) > 1)
            % 6 wavelength bands for g_ice to be interpolated into 480-bands of SNICAR
            % shape-preserving piecewise interpolation into 480-bands
            g_Cg_intp = pchip(g_wvl_center,g_ice_Cg_tmp,wvl) ;
            gg_F07_intp = pchip(g_wvl_center,gg_ice_F07_tmp,wvl) ;
            g_ice_F07 = gg_F07_intp + (1.0 - gg_F07_intp) ./ omega_ice(:,n)' ./ 2; % Eq.2.2 in Fu (2007)
            g_ice(:,n) = g_ice_F07 .* g_Cg_intp; % Eq.6, He et al. (2017)
            g_ice(381:480,n) = g_ice(380); % assume same values for 4-5 um band, with very small biases (<3%)
        end;
        
        g_ice(g_ice > 0.99) = 0.99; % avoid unreasonable values (so far only occur in large-size spheroid cases)
        
        
        % Volume and number concentration are not relevant for snow layers set to NAN
        No(n) = NaN;
        vlm_air(n) = NaN;
        
        % Specific surface area snow [m2/kg]
        ssa(n) = 3/(rho_ice*(rds_snw(n)*(10^-6))); 

    % ice with bubbles
    elseif lyr_typ(n) == 2 % ice 
        fl_ice              = strcat(dir_bbl,stb2,'_',sprintf('%04d',rds_snw(n)),'.nc');
        sca_cff_vlm         = ncread(fl_ice,'sca_cff_vlm'); % scattering cross section unit per volume of bubble
        g_ice(:,n)          = ncread(fl_ice,'asm_prm');     % asymmetry parameter
        abs_cff_mss_ice     = ((4 * pi * rfidx_ice_im) ./ (wvl * 1E-6))/rho_ice; % 
        vlm_frac_air        = (rho_ice - rho_snw(n)) / rho_ice; % (rho_snw(n) - rho_ice) / (-rho_ice + rho_air); 
        ext_cff_mss_ice(:,n)= ((sca_cff_vlm * vlm_frac_air) /rho_snw(n)) + abs_cff_mss_ice';
        omega_ice(:,n)      = ((sca_cff_vlm * vlm_frac_air) /rho_snw(n)) ./ ext_cff_mss_ice(:,n);
        
        % air bubble size distribution properties  
         sigma_tilde_g = log(1.5);
         d_eff_m = (rds_snw(n)*2)*(10^-6);
        
        % number concentration of bubbles [bubbles/m3]    
        No(n) = (6*vlm_frac_air)./(pi*d_eff_m.^3)*exp(3*sigma_tilde_g^2);
        
        % Volume of air in ice [unitless: m3/m-3]
        vlm_air(n) = vlm_frac_air;
        
        % Specific surface area ice [m2/kg]
        ssa(n) = (3*vlm_frac_air)/(rho_snw(n)*(rds_snw(n)*(10^-6)));
        
    end
    
    % algae, based on cell size and pigment concentrations:
    if (snw_alg_cell_nbr_conc(n) > 0)
        fl_alg = strcat(dir_alg,...
            'alg_sph_r',sprintf('%03d',round(alg_rds(n))),'um_',...
            'chla',sprintf('%03d',round(dcmf_pig_chla(n)*1000)),'_',...
            'chlb',sprintf('%03d',round(dcmf_pig_chlb(n)*1000)),'_',...
            'cara',sprintf('%03d',round(dcmf_pig_cara(n)*1000)),'_',...
            'carb',sprintf('%03d',round(dcmf_pig_carb(n)*1000)),...
            '.nc');
        
        % check that derived alg-impurity file exists:
        if (exist(fl_alg,'file')~=2)
            error(strcat('Snow algae impurity file: ',fl_alg,' does not exist.'));
        end;
        
        % single-scatter albedo, mass extinction coefficient, and
        % asymmetry paramater of snow algae:
        omega_snw_alg(:,n)       = ncread(fl_alg,'ss_alb');
        ext_cff_mss_snw_alg(:,n) = ncread(fl_alg,'ext_cff_mss');
        g_snw_alg(:,n)           = ncread(fl_alg,'asm_prm');
    end;
    
end % end loop over column layer

% Read Mie LAI parameters (layer-independent)

% read NetCDF properties for LAC 
for j=1:nbr_aer
    fl_in                = char(files(j));
    omega_aer(:,j)       = ncread(fl_in,'ss_alb');
    g_aer(:,j)           = ncread(fl_in,'asm_prm');
    if ((j==2) | (j==4))
        % unique variable name for mass extinction coefficient of
        % coated aerosols (i.e., coated BC and coated BrC)
        ext_cff_mss_aer(:,j) = ncread(fl_in,'ext_cff_mss_ncl');
    else
        ext_cff_mss_aer(:,j) = ncread(fl_in,'ext_cff_mss');
    end;
end

% Set aerosol mass mixing ratio (units of kg/kg):
mss_cnc_aer(1:nbr_lyr,1)  = mss_cnc_sot1.*1E-9;
mss_cnc_aer(1:nbr_lyr,2)  = mss_cnc_sot2.*1E-9;
mss_cnc_aer(1:nbr_lyr,3)  = mss_cnc_brc1.*1E-9;
mss_cnc_aer(1:nbr_lyr,4)  = mss_cnc_brc2.*1E-9;

mss_cnc_aer(1:nbr_lyr,5)  = mss_cnc_dst1.*1E-6;
mss_cnc_aer(1:nbr_lyr,6)  = mss_cnc_dst2.*1E-6;
mss_cnc_aer(1:nbr_lyr,7)  = mss_cnc_dst3.*1E-6;
mss_cnc_aer(1:nbr_lyr,8)  = mss_cnc_dst4.*1E-6;
mss_cnc_aer(1:nbr_lyr,9)  = mss_cnc_dst5.*1E-6;

mss_cnc_aer(1:nbr_lyr,10)  = mss_cnc_ash1.*1E-6;
mss_cnc_aer(1:nbr_lyr,11)  = mss_cnc_ash2.*1E-6;
mss_cnc_aer(1:nbr_lyr,12)  = mss_cnc_ash3.*1E-6;
mss_cnc_aer(1:nbr_lyr,13)  = mss_cnc_ash4.*1E-6;
mss_cnc_aer(1:nbr_lyr,14)  = mss_cnc_ash5.*1E-6;

mss_cnc_aer(1:nbr_lyr,15)  = glc_alg_mss_cnc.*1E-9;


% Calculate effective tau, omega, g for the (ice+algae+impurity) system
for n=1:nbr_lyr
    
    % Snow column mass [kg/m^2] (array)
    % Mass of snow is ice+impurities
    L_snw(n)     = rho_snw(n)*dz(n);
    
     % burden and optical thickness of algae
     if (snw_alg_cell_nbr_conc(n) > 0)
         % mean algal cell volume (3rd moment of Gaussian distribution) by layer:
         mean_vol_cell  = 4/3*pi * (alg_rds(n)^3 + 3*alg_rds(n)*sigma_alg(n)^2); %[um^3/cell]
         
         % mean mass per cell by layer:
         mass_per_cell  = mean_vol_cell*1E-18*rho_alg; % [kg/cell]
         
         % mass concentration of algae by layer
         mss_cnc_alg(n) = snw_alg_cell_nbr_conc(n)*1000*mass_per_cell; % [kg/kg]
         
         L_alg(n)     = L_snw(n)*mss_cnc_alg(n);
         tau_alg(:,n) = L_alg(n).*ext_cff_mss_snw_alg(:,n);
     else
         L_alg(n)     = 0.0;
         tau_alg(:,n) = 0.0;
     end;
    
    % burdens and optical thicknesses of LAI
    for j=1:nbr_aer
        L_aer(n,j)     = L_snw(n)*mss_cnc_aer(n,j);
        tau_aer(:,n,j) = L_aer(n,j).*ext_cff_mss_aer(:,j);
    end
    
    % ice mass = snow mass - impurity mass (generally a tiny correction)
    L_ice(n)     =  L_snw(n) - L_alg(n) - sum(L_aer(n,:));
    
    if (L_ice(n) < 0)
        error(['Impurity load cannot exceed snow load. Snow mass ' ...
               'is assumed to be that of ice+impurities, so the sum ' ...
               'of impurity mixing ratios cannot exceed 1']);
    end;
    
    % optical thickness due to ice:
    tau_ice(:,n) = L_ice(n).*ext_cff_mss_ice(:,n);
    
    tau_sum(1:nbr_wvl,1)   = 0.0;
    omega_sum(1:nbr_wvl,1) = 0.0;
    g_sum(1:nbr_wvl,1)     = 0.0;
    
    for j=1:nbr_aer
        tau_sum   = tau_sum + tau_aer(:,n,j);
        omega_sum = omega_sum + (tau_aer(:,n,j).*omega_aer(:,j));
        g_sum     = g_sum + (tau_aer(:,n,j).*omega_aer(:,j).*g_aer(:,j));
    end
  
    % add algae contribution to weighted LAI sums:
    if (snw_alg_cell_nbr_conc(n) > 0)
        tau_sum   = tau_sum + tau_alg(:,n);
        omega_sum = omega_sum + (tau_alg(:,n).*omega_snw_alg(:,n));
        g_sum     = g_sum + (tau_alg(:,n).*omega_snw_alg(:,n).*g_snw_alg(:,n));
    end;
    
    % weighted sums, including ice:
    tau(:,n)   = tau_sum + tau_ice(:,n);
    omega(:,n) = (1./tau(:,n)).*(omega_sum+ (omega_ice(:,n).*tau_ice(:,n)));
    g(:,n)     = (1./(tau(:,n).*omega(:,n))) .* (g_sum+ (g_ice(:,n).*omega_ice(:,n).*tau_ice(:,n)));
end


% ----- BEGIN Radiative Solver Adding Doubling Method -----
tau0    = tau;
g0      = g;
omega0  = omega;

epsilon = 1e-5;     %to deal with singularity in alpha and gamma 
exp_min = 1e-5;     %exp(-500); % minimum number that is not zero - zero will raise error
trmin   = 1e-4;     %minimum transmissivity
puny    = 1e-10;    


gauspt = [0.9894009, 0.9445750, 0.8656312, 0.7554044, ... % gaussian angles (radians)
    0.6178762, 0.4580168, 0.2816036, 0.0950125];
gauswt = [0.0271525, 0.0622535, 0.0951585, 0.1246290, ... % gaussian weights
    0.1495960, 0.1691565, 0.1826034, 0.1894506];

for iw = 1: nbr_wvl;
    for k = 1: nbr_lyr+1;   %klevp
        trndir(iw,k) = 0;   %solar beam down transmission from top
        trntdr(iw,k) = 0;   %total transmission from layers above
        trndif(iw,k) = 0;   %diffuse transmission for layers above
        rupdir(iw,k) = 0;   %reflectivity to direct radiation for layers below
        rupdif(iw,k) = 0;   %reflectivity to diffuse radiation for layers below
        rdndif(iw,k) = 0;   %reflectivity to diffuse radiation for layers above
    end
    
    %! initialize top interface of top layer
    trndir(iw,1) =  1;
    trntdr(iw,1) =  1;
    trndif(iw,1) =  1;
    rdndif(iw,1) =  0;
    
end;

% ! proceed down one layer at a time; if the total transmission to
% ! the interface just above a given layer is less than trmin, then no
% ! Delta-Eddington computation for that layer is done.
for iw = 1:nbr_wvl      %wavelengths
    
    %! begin main level loop
    for k = 1:nbr_lyr   %number of layers
        
        %! initialize all layer apparent optical properties to 0
        rdir  (k) = 0;  %layer reflectivity to direct radiation
        rdif_a(k) = 0;  %layer reflectivity to diffuse radiation from above
        rdif_b(k) = 0;  %layer reflectivity to diffuse radiation from below
        tdir  (k) = 0;  %layer transmission to direct radiation (solar beam + diffuse)
        tdif_a(k) = 0;  %layer transmission to diffuse radiation from above
        tdif_b(k) = 0;  %layer transmission to diffuse radiation from below
        trnlay(k) = 0;  %solar beam transm for layer (direct beam only)
        
        %    ! compute next layer Delta-eddington solution only if total transmission
        %    ! of radiation to the interface just above the layer exceeds trmin.
        
        if (trntdr(iw,k) > trmin )
            %    ! initialize current layer properties to zero; only if total
            %    ! transmission to the top interface of the current layer exceeds the
            %    ! minimum, will these values be computed below:
            mu0 = mu_not;
          
            % ice adjusted refractive index (Liou 2002)
            nr = Nreal(iw);

            if(k < kfrsnl || kfrsnl==0)
                % above FL mu0 is unchanged 
                mu0n = mu0;
            elseif(k >= kfrsnl)
                % mu0 under the Fresnel Layer
                % Eq. (5.4.13) Liou 2002
                mu0n = cos(asin(sin(acos(mu0))/nr)); 
            end
            
            %! calculation over layers with penetrating radiation
            tautot = tau0(iw,k);
            wtot   = omega0(iw,k);
            gtot   = g0(iw,k);
            ftot   = g0(iw,k) * g0(iw,k);
            
            % coefficient for delta eddington solution for all layers;
            % Eq. 50: Briegleb and Light 2007
            ts   = (1-wtot.*ftot) .* tautot;         %layer delta-scaled extinction optical depth
            ws   = ((1-ftot).*wtot)./(1-wtot.*ftot); %layer delta-scaled single scattering albedo
            gs   = gtot/(1+gtot);                    %layer delta-scaled asymmetry parameter 
            lm   = sqrt(3*(1-ws) .* (1-ws.*gs));     %lambda
            ue   = 1.5 * (1-ws.*gs)./lm;             %u equation, term in diffuse reflectivity and transmissivity
           
            extins = max(exp_min, exp(-lm*ts));       %extinction, MAX function keeps from getting an error if the exp(-lm*ts) is < 1e-5
            ne     = (ue+1).^2./extins - (ue-1).^2.*extins; %N equation, term in diffuse reflectivity and transmissivity
            
            % ! first calculation of rdif, tdif using Delta-Eddington formulas
            % Eq.: Briegleb 1992; alpha and gamma for direct radiation
            rdif_a(k) = (ue.^2-1)*(1/extins - extins)/ne;   %R BAR = layer reflectivity to DIFFUSE radiation
            tdif_a(k) = 4*ue/ne;                            %T BAR layer transmissivity to DIFFUSE radiation
            
            %   ! evaluate rdir,tdir for direct beam
            trnlay(k) = max(exp_min, exp(-ts/mu0n));        % transmission from TOA to interface
            
            %  Eq. 50: Briegleb and Light 2007; alpha and gamma for direct radiation
            alp = (0.75.*ws.*mu0n) ...
                .* (1 + gs.*(1-ws))...
                ./ (1 - lm.^2 .* mu0n.^2 + epsilon); %alp = alpha(ws,mu0n,gs,lm)
            gam = (0.5 .* ws)...
                .* (1 + 3.*gs.*mu0n.^2.*(1-ws))...
                ./ (1-lm.^2 .* mu0n.^2 + epsilon);   %gam = gamma(ws,mu0n,gs,lm)
            apg = alp + gam;
            amg = alp - gam;
            rdir(k) = apg*rdif_a(k) +  amg*(tdif_a(k)*trnlay(k) - 1);   %layer reflectivity to DIRECT radiation
            tdir(k) = apg*tdif_a(k) + (amg* rdif_a(k)-apg+1)*trnlay(k); %layer transmissivity to DIRECT radiation
            
            %    ! recalculate rdif,tdif using direct angular integration over rdir,tdir,
            %    ! since Delta-Eddington rdif formula is not well-behaved (it is usually
            %    ! biased low and can even be negative); use ngmax angles and gaussian
            %    ! integration for most accuracy:
            
            R1 = rdif_a(k); %! use R1 as temporary
            T1 = tdif_a(k); %! use T1 as temporary
            swt = 0;
            smr = 0;
            smt = 0;
            
            % loop through the gaussian angles for the AD integral
            for ng=1:length(gauspt)     %gaussian angles (radians)
                mu  = gauspt(ng);       %solar zenith angles
                gwt = gauswt(ng);       %gaussian weight
                swt = swt + mu*gwt;     % sum of weights
                trn = max(exp_min, exp(-ts/mu)); %transmission
                
                alp = (0.75.*ws.*mu) ...
                    .* (1 + gs.*(1-ws))...
                    ./ (1 - lm.^2 .* mu.^2 + epsilon); %alp = alpha(ws,mu0n,gs,lm)
                gam = (0.5 .* ws)...
                    .* (1 + 3.*gs.*mu.^2.*(1-ws))...
                    ./ (1-lm.^2 .* mu.^2 + epsilon);%gam = gamma(ws,mu0n,gs,lm)
                
                apg = alp + gam;
                amg = alp - gam;
                rdr = apg*R1 + amg*T1*trn - amg;
                tdr = apg*T1 + amg*R1*trn - apg*trn + trn;
                smr = smr + mu*rdr*gwt; %accumulator for rdif gaussian integration
                smt = smt + mu*tdr*gwt; %accumulator for tdif gaussian integration
            end      %! ng; gaussian angles for the AD integral
            
            rdif_a(k) = smr/swt; 
            tdif_a(k) = smt/swt;
            
            %! homogeneous layer
            rdif_b(k) = rdif_a(k);
            tdif_b(k) = tdif_a(k);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fresnel layer
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (k == kfrsnl)
                
                % ice complex index of refraction 
                refindx = complex(rfidx_ice_re(iw),rfidx_ice_im(iw));
                
                % critical angle where total internal reflection occurs 
                critical_angle = asin(refindx);
                
                % compare incoming angle to critical angle 
                if acos(mu_not) < critical_angle
                    
                %! compute fresnel reflection and transmission amplitudes
                %! for two polarizations: 1=perpendicular and 2=parallel to
                %! the plane containing incident, reflected and refracted rays.
                    
                    %! Eq. (5.4.18a-b); Liou 2002
                    R1 = (mu0-nr*mu0n) / (mu0 + nr*mu0n);      %reflection amplitude factor for perpendicular polarization
                    R2 = (nr*mu0 - mu0n) / (nr*mu0 + mu0n);    %reflection amplitude factor for parallel polarization
                    T1 = 2*mu0 / (mu0 + nr*mu0n);              %transmission amplitude factor for perpendicular polarization
                    T2 = 2*mu0 / (nr*mu0 + mu0n);              %transmission amplitude factor for parallel polarization
                    
                    %! unpolarized light for direct beam
                    %! Eq. 21; Brigleb and light 2007
                    Rf_dir_a = 0.5 * ((R1^2) + (R2^2));
                    Tf_dir_a = 0.5 * (T1*T1 + T2*T2)*nr*mu0n/mu0;
               
                else % total internal reflection 
                    Tf_dir_a = 0;
                    Rf_dir_a = 1;
                    
                end % if critical angle check 
                
                %      ! precalculated diffuse reflectivities and transmissivities
                %      ! for incident radiation above and below fresnel layer, using
                %      ! the direct albedos and accounting for complete internal
                %      ! reflection from below; precalculated because high order
                %      ! number of gaussian points is required for convergence:
                
                Rf_dif_a = FL_r_dif_a(iw);
                Tf_dif_a = 1 - Rf_dif_a;
                Rf_dif_b = FL_r_dif_b(iw);
                Tf_dif_b = 1 - Rf_dif_b;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %! the k = kfrsnl layer properties are updated to combined
                %! the fresnel (refractive) layer, always taken to be above
                %! the present layer k (i.e. be the top interface):
                
                rintfc   = 1 / (1-Rf_dif_b*rdif_a(k)); % denom interface scattering
                
                tdir(k)   = Tf_dir_a*tdir(k) + ...    % layer transmissivity to DIRECT radiation
                    Tf_dir_a*rdir(k) * ...            % Eq. B7; Briegleb & Light 2007
                    Rf_dif_b*rintfc*tdif_a(k);
                
                rdir(k)   = Rf_dir_a + ...            % layer reflectivity to DIRECT radiation
                    Tf_dir_a*rdir(k) * ...            % Eq. B7; Briegleb & Light 2007
                    rintfc*Tf_dif_b;
                
                rdif_a(k) = Rf_dif_a + ...            % R BAR = layer reflectivity to DIFFUSE radiation (above)
                    Tf_dif_a*rdif_a(k) * ...          % Eq. B9; Briegleb & Light 2007
                    rintfc*Tf_dif_b;
                
                rdif_b(k) = rdif_b(k) + ...           % R BAR = layer reflectivity to DIFFUSE radiation (below)
                    tdif_b(k)*Rf_dif_b * ...          % Eq. B10; Briegleb & Light 2007
                    rintfc*tdif_a(k);
                
                tdif_a(k) = tdif_a(k)*rintfc*Tf_dif_a;  % T BAR layer transmissivity to DIFFUSE radiation (above), Eq. B9; Briegleb & Light 2007
                tdif_b(k) = tdif_b(k)*rintfc*Tf_dif_b;  % Eq. B10; Briegleb & Light 2007
                
                %! update trnlay to include fresnel transmission
                trnlay(k) = Tf_dir_a*trnlay(k);
                
            end     %k = kfrsnl
           
            
        end %! trntdr(k,iw) > trmin
        
        %       ! Calculate the solar beam transmission, total transmission, and
        %       ! reflectivity for diffuse radiation from below at interface k,
        %       ! the top of the current layer k:
        %       !
        %       !              layers       interface
        %       !
        %       !       ---------------------  k-1
        %       !                k-1
        %       !       ---------------------  k
        %       !                 k
        %       !       ---------------------
        %       ! note that we ignore refraction between sea ice and underlying ocean:
        %       !
        %       !              layers       interface
        %       !
        %       !       ---------------------  k-1
        %       !                k-1
        %       !       ---------------------  k
        %       !       \\\\\\\ ocean \\\\\\\
        
        % Eq. 51; Briegleb and Light 2007
        trndir(iw,k+1) = trndir(iw,k)*trnlay(k); % solar beam transmission from top
                                    % trnlay = exp(-ts/mu_not) = direct solar beam transmission
        refkm1         = 1/(1 - rdndif(iw,k)*rdif_a(k)); %interface multiple scattering for k-1
        
        tdrrdir        = trndir(iw,k)*rdir(k);           %direct tran times layer direct ref
        
        tdndif         = trntdr(iw,k) - trndir(iw,k);    %total down diffuse = tot tran - direct tran
        
        trntdr(iw,k+1) = trndir(iw,k)*tdir(k) + ...
            (tdndif + tdrrdir*rdndif(iw,k))*refkm1*tdif_a(k); %total transmission for layers above
        
        %Eq. 51; Briegleb and Light 2007
        rdndif(iw,k+1) = rdif_b(k) + ...
            (tdif_b(k)*rdndif(iw,k)*refkm1*tdif_a(k));   %reflectivity to diffuse radiation for layers above
        trndif(iw,k+1) = trndif(iw,k)*refkm1*tdif_a(k);  %diffuse transmission to diffuse beam for layers above
        
    end       %! k   end main level loop; number of layers
    
    % ! compute reflectivity to direct and diffuse radiation for layers
    % ! below by adding succesive layers starting from the underlying
    % ! ocean and working upwards:
    % !
    % !              layers       interface
    % !
    % !       ---------------------  k
    % !                 k
    % !       ---------------------  k+1
    % !                k+1
    % !       ---------------------
    
    % set the underlying ground albedo
    rupdir(iw,nbr_lyr+1) = R_sfc(iw);   % reflectivity to direct radiation for layers below
    rupdif(iw,nbr_lyr+1) = R_sfc(iw);   % reflectivity to diffuse radiation for layers below
    
    for k=nbr_lyr:-1:1  % starts at the bottom and works its way up to the top layer
        
        %Eq. B2; Briegleb and Light 2007
        %! interface scattering
        refkp1        = 1/( 1 - rdif_b(k)*rupdif(iw,k+1));
        
        %   ! dir from top layer plus exp tran ref from lower layer, interface
        %   ! scattered and tran thru top layer from below, plus diff tran ref
        %   ! from lower layer with interface scattering tran thru top from below
        rupdir(iw,k) = rdir(k) ...
            + (        trnlay(k)  * rupdir(iw,k+1) ...
            +  (tdir(k)-trnlay(k))* rupdif(iw,k+1))*refkp1*tdif_b(k);
        %   ! dif from top layer from above, plus dif tran upwards reflected and
        %   ! interface scattered which tran top from below
        rupdif(iw,k) = rdif_a(k) + tdif_a(k)*rupdif(iw,k+1)*refkp1*tdif_b(k);
    end      %! k
    
end      %! iw; number of wavelengths

% fluxes at interface

for iw = 1:nbr_wvl
    for k=1:nbr_lyr+1
        
        % Eq. 52; Briegleb and Light 2007
        %! interface scattering
        refk          = 1/(1 - rdndif(iw,k)*rupdif(iw,k));
        %     ! dir tran ref from below times interface scattering, plus diff
        %     ! tran and ref from below times interface scattering
        fdirup(iw,k) = (trndir(iw,k)*rupdir(iw,k) + ...
            (trntdr(iw,k)-trndir(iw,k))  ...
            *rupdif(iw,k))*refk;
        %     ! dir tran plus total diff trans times interface scattering plus
        %     ! dir tran with up dir ref and down dif ref times interface scattering
        fdirdn(iw,k) = trndir(iw,k) + (trntdr(iw,k) ...
            - trndir(iw,k) + trndir(iw,k)  ...
            *rupdir(iw,k)*rdndif(iw,k))*refk;
        %     ! diffuse tran ref from below times interface scattering
        fdifup(iw,k) = trndif(iw,k)*rupdif(iw,k)*refk;
        %     ! diffuse tran times interface scattering
        fdifdn(iw,k) = trndif(iw,k)*refk;
        %
        %     ! dfdir = fdirdn - fdirup
        dfdir(iw,k) = trndir(iw,k) ...
            + (trntdr(iw,k)-trndir(iw,k)) * (1 - rupdif(iw,k)) * refk ...
            -  trndir(iw,k)*rupdir(iw,k)  * (1 - rdndif(iw,k)) * refk;
        if (dfdir(iw,k) < puny)
            dfdir(iw,k) = 0; %!echmod necessary?
            %! dfdif = fdifdn - fdifup
        end
        dfdif(iw,k) = trndif(iw,k) * (1 - rupdif(iw,k)) * refk;
        if (dfdif(iw,k) < puny)
            dfdif(iw,k) = 0; %!echmod necessary?
        end
    end      %! k
end

% ----- END Radiative Solver Adding Doubling Method -----


% ----- Calculate and package radiative terms ----

% Upward and downward fluxes at each layer interface:
for n = 1:nbr_lyr+1
    F_up(:,n)  = (fdirup(:,n).*(Fs*mu_not*pi) + fdifup(:,n).*Fd);
    F_dwn(:,n) = (fdirdn(:,n).*(Fs*mu_not*pi) + fdifdn(:,n).*Fd);
end

% Net flux at each layer interface:
F_net = F_up - F_dwn;

% Absorbed flux in each layer:
F_abs = F_net(:,2:end) - F_net(:,1:end-1);

% Upward flux at upper model boundary
F_top_pls = F_up(:,1);

% Net flux at lower model boundary = bulk transmission through entire
% media = absorbed radiation by underlying surface:
F_btm_net = -F_net(:,nbr_lyr+1);

% Spectrally-integrated absorption in each layer:
F_abs_slr = sum(F_abs);
for n=1:nbr_lyr
    F_abs_vis(n) = sum(F_abs(1:vis_max_idx,n));
    F_abs_nir(n) = sum(F_abs(vis_max_idx+1:nir_max_idx,n));
end

% Spectrally-integrated absorption by underlying surface:
F_abs_slr_btm = sum(F_btm_net);
F_abs_vis_btm = sum(F_btm_net(1:vis_max_idx));
F_abs_nir_btm = sum(F_btm_net(vis_max_idx+1:nir_max_idx));


% Radiative heating rate:
heating_rate = F_abs_slr./(L_snw.*2117);   %[K/s], 2117 = specific heat ice (J kg-1 K-1)
heating_rate = heating_rate.*3600;         %[K/hr]

% Energy conservation check:
% Incident direct+diffuse radiation equals (absorbed+transmitted+bulk_reflected)
energy_sum = (mu_not*pi*Fs)+Fd - (sum(F_abs,2) + F_btm_net + F_top_pls);

energy_conservation_error = sum(abs(energy_sum));
if (energy_conservation_error > 1e-10)
    disp(['energy conservation error: ' num2str(energy_conservation_error)])
end

% Spectral albedo:
albedo  = F_up(:,1)./F_dwn(:,1);

% Alt spectral albedo (identical to previous within ~machine precision)
%if (direct_beam==1)
%    albedo = rupdir(:,1)'; % F_top_pls./((mu_not*pi*Fs)+Fd);
%else
%    albedo = rupdif(:,1)';
%end

% Spectrally-integrated solar, visible, and NIR albedos:
alb_bb  = sum(flx_slr.*albedo)./sum(flx_slr);

alb_vis = sum(flx_slr(1:vis_max_idx).*albedo(1:vis_max_idx))/...
          sum(flx_slr(1:vis_max_idx));

alb_nir = sum(flx_slr(vis_max_idx+1:nir_max_idx).*albedo(vis_max_idx+1:nir_max_idx))/...
          sum(flx_slr(vis_max_idx+1:nir_max_idx));


% Spectrally-integrated VIS and NIR total snowpack absorption:
abs_vis = sum(flx_slr(1:vis_max_idx).*(1-albedo(1:vis_max_idx)));
abs_nir = sum(flx_slr(vis_max_idx+1:nir_max_idx).*(1-albedo(vis_max_idx+1:nir_max_idx)));

% Spectrally-integrated downwelling solar fluxes at top of snowpack [W/m2]:
flx_dwn_top_slr  = sum((mu_not*pi*Fs))+sum(Fd);
flx_dwn_top_vis  = sum((mu_not*pi*Fs(1:vis_max_idx)))+sum(Fd(1:vis_max_idx));
flx_dwn_top_nir  = sum((mu_not*pi*Fs(vis_max_idx+1:nir_max_idx)))+sum(Fd(vis_max_idx+1:nir_max_idx));


%%%%%%%%%%%%%%%%%%%%%%%%%  OUTPUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_out.wvl             = wvl;                % spectral wavelength bands (um)
data_out.flx_dwn_spc     = (mu_not*pi*Fs+Fd)'; % spectral downwelling flux at model top [W/m2/band]
data_out.albedo          = albedo';            % spectral hemispheric albedo
data_out.alb_slr         = alb_bb;             % solar broadband albedo
data_out.alb_vis         = alb_vis;            % visible (0.2-0.7um) albedo
data_out.alb_nir         = alb_nir;            % near-IR (0.7-5.0um) albedo

data_out.abs_snw_slr     = sum(F_abs_slr);     % total solar absorption by entire snow column (not including underlying substrate) [W/m2]
data_out.abs_snw_vis     = sum(F_abs_vis);     % visible solar absorption by entire snow column (not including underlying substrate) [W/m2]
data_out.abs_snw_nir     = sum(F_abs_nir);     % near-IR solar absorption by entire snow column (not including underlying substrate) [W/m2]
data_out.abs_spc         = sum(F_abs,2)';      % spectral absorption by entire snow column [W/m2/band]

data_out.abs_snw_top_slr = F_abs_slr(1);       % top snow layer solar absorption [W/m2]
data_out.abs_snw_top_vis = F_abs_vis(1);       % top snow layer VIS absorption [W/m2]
data_out.abs_snw_top_nir = F_abs_nir(1);       % top snow layer NIR absorption [W/m2]

data_out.abs_ground_slr  = F_abs_slr_btm;      % total solar absorption by underlying substrate [W/m2]
data_out.abs_ground_vis  = F_abs_vis_btm;      % visible absorption by underlying substrate [W/m2]
data_out.abs_ground_nir  = F_abs_nir_btm;      % near-IR absorption by underlying substrate [W/m2]

data_out.flx_dwn_top_slr = flx_dwn_top_slr;    % downwelling solar broadband flux on upper boundary [W/m2]
data_out.flx_dwn_top_vis = flx_dwn_top_vis;    % downwelling visible broadband flux on upper boundary [W/m2]
data_out.flx_dwn_top_nir = flx_dwn_top_nir;    % downwelling visible broadband flux on upper boundary [W/m2]

data_out.vlm_frc_air    = vlm_air;            % volume fraction of air in each ice layer [unitless m3/m3]
data_out.No_bbl_cnc     = No;                 % number concentration of air bubbles in ice [bbls/m3]
data_out.ssa            = ssa;                % specific surface area of snow or ice [m2/kg]


if (0==1)
    figure(1)
    hold on
    plot(data_out.wvl,data_out.albedo,'linewidth',3,'DisplayName','SNICAR-ADv3');
    axis([0.2 1.8 0 1]);
    set(gca,'xtick',0.2:0.2:1.8,'fontsize',14)
    set(gca,'ytick',0:0.1:1,'fontsize',14);
    xlabel('Wavelength (\mum)','fontsize',20);
    ylabel('Hemispheric Albedo','fontsize',20);
    grid on;
    legend();
end;


end % end of function


%%%%%%%%%%%%%%%%%%%%%%%%%  REFERENCES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Balkanski, Y., M. Schulz, T. Claquin, and S. Guibert (2007),
% Reevaluation of mineral aerosol radiative forcings suggests a better
% agreement with satellite and AERONET data, Atmos. Chem. Phys., 7
% (1), 81-95, doi:10.5194/acp-7-81-2007.
%
% Briegleb, P. and Light, B. (2007) A Delta-Eddington mutiple
% scattering parameterization for solar radiation in the sea ice
% component of the community climate system model, NCAR Technical Note
% NCAR/TN-472+STR, https://doi.org/10.5065/D6B27S71.
%
% Cook, J. M., Hodson, A. J., Gardner, A. S., Flanner, M., Tedstone,
% A. J., Williamson, C., Irvine-Fynn, T. D. L., Nilsson, J., Bryant,
% R., and Tranter, M. (2017), Quantifying bioalbedo: a new physically
% based model and discussion of empirical methods for characterising
% biological influence on ice and snow albedo, The Cryosphere, 11,
% 2611-2632, doi: 10.5194/tc-11-2611-2017.
% 
% Cook, J. M., Tedstone, A. J., Williamson, C., McCutcheon, J., Hodson, 
% A. J., Dayal, A., Skiles, M., Hofer, S., Bryant, R., McAree, 
% O., McGonigle, A., Ryan, J., Anesio, A. M., Irvine-Fynn, 
% T. D. L., Hubbard, A., Hanna, E., Flanner, M., Mayanna, 
% S., Benning, L. G., ? Tranter, M. (2020). 
% Glacier algae accelerate melt rates on the south-western Greenland 
% Ice Sheet. The Cryosphere, 14(1), 309?330. 
% https://doi.org/10.5194/tc-14-309-2020
%
% Dang, C., Zender, C. S., and Flanner, M. G. (2019), Intercomparison
% and improvement of two-stream shortwave radiative transfer schemes
% in Earth system models for a unified treatment of cryospheric
% surfaces, The Cryosphere, 13, 2325-2343,
% doi:10.5194/tc-13-2325-2019.
% 
% Flanner, M. G., Arnheim, J., Cook, J. M., Dang, C., He, C., Huang, X., 
% Singh, D., Skiles, S. M., Whicker, C. A., and Zender, C. S.: SNICAR-AD v3: 
% A Community Tool for Modeling Spectral Snow Albedo, 1?49, 
% https://doi.org/10.5194/gmd-2021-182, 2021.
%
% Flanner, M. G., A. S. Gardner, S. Eckhardt, A. Stohl, J. Perket
% (2014), Aerosol radiative forcing from the 2010 Eyjafjallajokull
% volcanic eruptions, J. Geophys. Res. Atmos., 119, 9481-9491,
% doi:10.1002/2014JD021977.
%
% Flanner, M. G., C. S. Zender, J. T. Randerson, and P. J. Rasch
% (2007), Present day climate forcing and response from black carbon
% in snow, J. Geophys. Res., 112, D11202, doi:10.1029/2006JD008003.
%
% Fu, Q. (2007), A new parameterization of an asymmetry factor of
% cirrus clouds for climate models, J. Atmos. Sci., 64, 4140-4150,
% doi:10.1175/2007JAS2289.1.
%
% Hansen, G. B. (2005), Ultraviolet to near-infrared absorption
% spectrum of carbon dioxide ice from 0.174 to 1.8 um,
% J. Geophys. Res., 110, E11003, doi:10.1029/2005JE002531.
%
% He, C., Y. Takano, K.-N. Liou, et al. (2017), Impact of snow grain
% shape and black carbon-snow internal mixing on snow optical
% properties: Parameterizations for climate models, J. Climate,
% 30(24), 10019-10036, doi:10.1175/JCLI-D-17-0300.1.
%
% Kirchstetter, T. W., T. Novakov, and P. V. Hobbs (2004), Evidence
% that the spectral dependence of light absorption by aerosols is
% affected by organic carbon, J. Geophys. Res., 109, D21208,
% doi:10.1029/2004JD004999.
%
% Lawrence, D. and others (2018), Technical Description of version 5.0
% of the Community Land Model (CLM),
% http://www.cesm.ucar.edu/models/cesm2/land/CLM50_Tech_Note.pdf
%
% Liou, K. N.: An Introduction to Atmospheric Radiation, Academic Press, 
% Amsterdam, 2002.
%
% Picard, G., Libois, Q., and Arnaud, L. (2016), Refinement of the ice
% absorption spectrum in the visible using radiance profile
% measurements in Antarctic snow, The Cryosphere, 10, 2655-2672,
% doi:10.5194/tc-10-2655-2016.
%
% Polashenski, C. M., J. E. Dibb, M. G. Flanner, J. Y. Chen,
% Z. R. Courville, A. M. Lai, J. J. Schauer, M. M. Shafer, and
% M. Bergin (2015), Neither dust nor black carbon causing apparent
% albedo decline in Greenland's dry snow zone: Implications for MODIS
% C5 surface reflectance, Geophys. Res. Lett., 42, 9319-9327,
% doi:10.1002/2015GL065912.
%
% Singh, D., and M. G. Flanner (2016), An improved carbon dioxide snow
% spectral albedo model: Application to Martian conditions,
% J. Geophys. Res. Planets, 121, 2037-2054, doi:10.1002/2016JE005040.
%
% Skiles, S. M., Painter, T., and Okin, G. (2017), A method to
% retrieve the spectral complex refractive index and single scattering
% optical properties of dust deposited in mountain snow, Journal of
% Glaciology. Vol. 63, 133-147, doi:10.1017/jog.2016.126.
%
% Warren, S. G. (1984), Optical constants of ice from the ultraviolet
% to the microwave, Appl. Opt., 23, 1206-1225.
%
% Warren, S. G., and R. E. Brandt (2008), Optical constants of ice
% from the ultraviolet to the microwave: A revised compilation,
% J. Geophys. Res., 113, D14220, doi:10.1029/2007JD009744.
%
% Wolff, M. J., M. D. Smith, R. T. Clancy, R. Arvidson, M. Kahre,
% F. Seelos IV, S. Murchie, and H. Savijarvi (2009), Wavelength
% dependence of dust aerosol single scattering albedo as observed by
% the Compact Reconnaissance Imaging Spectrometer, J. Geophys. Res.,
% 114, E00D04, doi:10.1029/2009JE003350.
%
% Wolff, M. J., R. T. Clancy, J. D. Goguen, M. C. Malin, and
% B. A. Cantor (2010), Ultraviolet dust aerosol properties as observed
% by MARCI, Icarus, 208(1), 143-155.



% Copyright (c) 2020 SNICAR contributors
%
% Permission is hereby granted, free of charge, to any person
% obtaining a copy of this software and associated documentation files
% (the "Software"), to deal in the Software without restriction,
% including without limitation the rights to use, copy, modify, merge,
% publish, distribute, sublicense, and/or sell copies of the Software,
% and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
%
% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
