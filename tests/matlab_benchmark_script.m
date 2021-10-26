% Test-suite of runs for comparison with Joe Cook's Python code.
% See email sent on 20211022

clear;

lyr_array     = [1,2];
density_array = [400, 500, 600, 700, 800];
reff_array    = [200, 400, 600, 80, 1000];
zen_array     = [30, 40, 50, 60, 70];
bc_array      = [500, 1000, 1500, 2000];
dz_array      = [[0.02,0.04,0.06, 0.08, 0.1],
    [0.04,0.06,0.08, 0.10, 0.15],
    [0.05,0.10,0.15, 0.2, 0.5],
    [0.15, 0.2, 0.25, 0.3, 0.5],
    [0.5,0.5,0.5,1,10]];


nbr_lyr = 5;  % number of snow layers (true in all test cases)

% 1= Direct-beam incident flux, 0= Diffuse incident flux
% NOTE that cloudy-sky spectral fluxes are loaded when direct_beam=0
input_args.direct_beam   = 1;


% ICE REFRACTIVE INDEX DATASET TO USE:
% 1=Warren (1984); 2=Warren and Brandt (2008); 3=Picard et al. (2016); 4=CO2 ice
input_args.ice_ri = 3;

% Snow grain shape option
% 1=sphere; 2=spheroid; 3=hexagonal plate; 4=koch snowflake
input_args.sno_shp(1:nbr_lyr)  = 1;

% Shape factor: ratio of nonspherical grain effective radii to that of equal-volume sphere
% 0=use recommended default value (He et al. 2017);
% others(0<fs<1)= use user-specified value
% only activated when sno_shp > 1 (i.e. nonspherical)
input_args.sno_fs(1:nbr_lyr)   = 0;

% Aspect ratio: ratio of grain width to length
% 0=use recommended default value (He et al. 2017);
% others(0.1<fs<20)= use user-specified value
% only activated when sno_shp > 1 (i.e. nonspherical)
input_args.sno_ar(1:nbr_lyr)   = 0;

% type of dust: 1=Sahara, 2=Colorado, 3=Greenland, 4=Mars
input_args.dust_type = 1;

% type of volcanic ash: 1 = Eyjafjallajokull
input_args.ash_type = 1;

% PARTICLE MASS MIXING RATIOS (ng/g or ug/g, depending on species)
% NOTE: This is mass of impurity per mass of snow
%  (i.e., mass of impurity / mass of ice+impurity)

input_args.mss_cnc_sot2(1:nbr_lyr)  = 0.0;    % sulfate-coated black carbon [ng/g]
input_args.mss_cnc_brc1(1:nbr_lyr)  = 0.0;    % uncoated brown carbon [ng/g]
input_args.mss_cnc_brc2(1:nbr_lyr)  = 0.0;    % sulfate-coated brown carbon [ng/g]
input_args.mss_cnc_dst1(1:nbr_lyr)  = 0.0;    % dust size 1 (r=0.05-0.5um) [ug/g]
input_args.mss_cnc_dst2(1:nbr_lyr)  = 0.0;    % dust size 2 (r=0.5-1.25um) [ug/g]
input_args.mss_cnc_dst3(1:nbr_lyr)  = 0.0;    % dust size 3 (r=1.25-2.5um) [ug/g]
input_args.mss_cnc_dst4(1:nbr_lyr)  = 0.0;    % dust size 4 (r=2.5-5.0um)  [ug/g]
input_args.mss_cnc_dst5(1:nbr_lyr)  = 0.0;    % dust size 5 (r=5.0-50um)   [ug/g]
input_args.mss_cnc_ash1(1:nbr_lyr)  = 0.0;    % volcanic ash size 1 (r=0.05-0.5um) [ug/g]
input_args.mss_cnc_ash2(1:nbr_lyr)  = 0.0;    % volcanic ash size 2 (r=0.5-1.25um) [ug/g]
input_args.mss_cnc_ash3(1:nbr_lyr)  = 0.0;    % volcanic ash size 3 (r=1.25-2.5um) [ug/g]
input_args.mss_cnc_ash4(1:nbr_lyr)  = 0.0;    % volcanic ash size 4 (r=2.5-5.0um)  [ug/g]
input_args.mss_cnc_ash5(1:nbr_lyr)  = 0.0;    % volcanic ash size 5 (r=5.0-50um)   [ug/g]

input_args.snw_alg_cell_nbr_conc(1:nbr_lyr) = 0.0;    % SNOW algae [UNITS: cells/mL]
input_args.alg_rds(1:nbr_lyr)               = 10;     % mean algae cell radius (um)
input_args.dcmf_pig_chla(1:nbr_lyr)         = 0.015;  % dry cell mass fraction of chlorophyll-a
input_args.dcmf_pig_chlb(1:nbr_lyr)         = 0.005;  % dry cell mass fraction of chlorophyll-b
input_args.dcmf_pig_cara(1:nbr_lyr)         = 0.05;   % dry cell mass fraction of photoprotective_carotenoids
input_args.dcmf_pig_carb(1:nbr_lyr)         = 0.00;   % dry cell mass fraction of photosynthetic_carotenoids  

input_args.glc_alg_mss_cnc(1:nbr_lyr)       = 0.0;  % GLACIER algae [UNITS ng/g] = ppb 
input_args.glc_alg_rds                      = 4;    % GLACIER algae radius [um]
input_args.glc_alg_len                      = 40;   % GLACIER algae length [um]


% REFLECTANCE OF SURFACE UNDERLYING SNOW: (Value is applied to all
% wavelengths. User can alternatively specify spectrally-dependent
% ground albedo in snicar_v3.m)
input_args.R_sfc_all_wvl = 0.25;

% ATMOSPHERIC PROFILE for surface-incident flux:
%     1 = mid-latitude winter
%     2 = mid-latitude summer
%     3 = sub-Arctic winter
%     4 = sub-Arctic summer
%     5 = Summit,Greenland (sub-Arctic summer, surface pressure of 796hPa)
%     6 = High Mountain (summer, surface pressure of 556 hPa)
%     7 = Top-of-atmosphere
% NOTE that clear-sky spectral fluxes are loaded when direct_beam=1,
% and cloudy-sky spectral fluxes are loaded when direct_beam=0
input_args.atm = 5;

% Broadband surface-incident solar flux [W/m2]:
%  (used to scale spectral flux fractions)
input_args.flx_dwn_bb = 1.0;


% TEST SUITE LOOP:
idx = 1;
for h = 1:length(lyr_array)
    for i=1:length(density_array)
        for j=1:length(reff_array)
            for k=1:length(zen_array)
                for l=1:length(bc_array)
                    for m=1:length(dz_array)
                        
                        % LAYER MEDIUM TYPE [ 1=snow, 2=ice]
                        input_args.lyr_typ(1:nbr_lyr) = lyr_array(h);
                        
                        % SNOW DENSITY FOR EACH LAYER (units: kg/m3)
                        input_args.rho_snw(1:nbr_lyr) = density_array(i);
                        
                        % SNOW GRAIN SIZE FOR EACH LAYER (units: microns):
                        input_args.rds_snw(1:nbr_lyr) = reff_array(j);
                        
                        % COSINE OF SOLAR ZENITH ANGLE FOR DIRECT-BEAM RT
                        input_args.coszen = cosd(zen_array(k));
                        
                        % uncoated black carbon [ng/g]
                        input_args.mss_cnc_sot1(1:nbr_lyr)  = bc_array(l);
                        
                        % SNOW LAYER THICKNESSES [m]:
                        input_args.dz = dz_array(m,:); % multi-layer snowpack
                        
                        
                        % CALL SNICAR WITH INPUT ARGUMENTS
                        di = snicarAD_v4(input_args);
                        
                        albedo_out(:,idx) = di.albedo;
                        alb_bb_out(idx)   = di.alb_slr;
                        
                        idx = idx+1;
                        
                    end
                end
            end
        end
    end
end

% set last row equal to bb albedo and save as csv file:
albedo_out(481,:) = alb_bb_out;

writematrix(albedo_out,'testsuite_20211022b.csv');

