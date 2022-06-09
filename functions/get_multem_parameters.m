function input_multislice = get_multem_parameters(nx, conv_angle, E0, gpu)

%     nx = tfm_pn_fact(g_max * 2 * lx, 3);

    %%%%%%%%%%%%%%%%%% Load multem default parameter %%%%%%%%$$%%%%%%%%%
    input_multislice = multem_input.parameters();          % Load default values;

    %%%%%%%%%%%%%%%%%%%%% Set system configuration %%%%%%%%%%%%%%%%%%%%%
    input_multislice.system_conf.precision = 1;                           % eP_Float = 1, eP_double = 2
    if gpu >= 0
        input_multislice.system_conf.device = 2;                          % eD_CPU = 1, eD_GPU = 2
        input_multislice.system_conf.cpu_nthread = 1; 
        input_multislice.system_conf.gpu_device = gpu;
    else
        input_multislice.system_conf.device = 1;
        input_multislice.system_conf.cpu_nthread = abs(gpu); 
    end  

    %%%%%%%%%%%%%%%%%%%% Set simulation experiment %%%%%%%%%%%%%%%%%%%%%
    % eTEMST_STEM=11, eTEMST_ISTEM=12, eTEMST_CBED=21, eTEMST_CBEI=22, eTEMST_ED=31, eTEMST_HRTEM=32, eTEMST_PED=41, eTEMST_HCTEM=42, eTEMST_EWFS=51, eTEMST_EWRS=52, 
    % eTEMST_EELS=61, eTEMST_EFTEM=62, eTEMST_ProbeFS=71, eTEMST_ProbeRS=72, eTEMST_PPFS=81, eTEMST_PPRS=82,eTEMST_TFFS=91, eTEMST_TFRS=92
    input_multislice.simulation_type = 51;

    %%%%%%%%%%%%%% Electron-Specimen interaction model %%%%%%%%%%%%%%%%%
    input_multislice.interaction_model = 1;              % eESIM_Multislice = 1, eESIM_Phase_Object = 2, eESIM_Weak_Phase_Object = 3
    input_multislice.potential_type = 6;                 % ePT_Doyle_0_4 = 1, ePT_Peng_0_4 = 2, ePT_Peng_0_12 = 3, ePT_Kirkland_0_12 = 4, ePT_Weickenmeier_0_12 = 5, ePT_Lobato_0_12 = 6

    %%%%%%%%%%%%%%%%%%%%%%% Potential slicing %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.potential_slicing = 2;              % ePS_Planes = 1, ePS_dz_Proj = 2, ePS_dz_Sub = 3, ePS_Auto = 4
    input_multislice.spec_dz = 1;
    %%%%%%%%%%%%%%% Electron-Phonon interaction model %%%%%%%%%%%%%%%%%%
    input_multislice.pn_model = 3;                       % ePM_Still_Atom = 1, ePM_Absorptive = 2, ePM_Frozen_Phonon = 3
    input_multislice.pn_coh_contrib = 0;                 % elastic / inelastic contributions
    input_multislice.pn_single_conf = 0;                 % 1: true, 0:false (extract single configuration)
    input_multislice.pn_nconf = 30;                      % true: specific phonon configuration, false: number of frozen phonon configurations
    input_multislice.pn_dim = 110;                       % phonon dimensions (xyz)
    input_multislice.pn_seed = 353183;                   % Random seed(frozen phonon)

    %%%%%%%%%%%%%%%%%%%%%% Specimen thickness %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.thick_type = 1;                     % eTT_Whole_Spec = 1, eTT_Through_Thick = 2, eTT_Through_Slices = 3
    % input_multislice.thick = c:c:1000;                   % Array of thickes (�)

    %%%%%%%%%%%%%%%%%%%%%% x-y sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.nx = nx;
    input_multislice.ny = nx;
    input_multislice.bwl = 0;                            % Band-width limit, 1: true, 0:false

    %%%%%%%%%%%%%%%%%%%% Microscope parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.E_0 = E0;                          % Acceleration Voltage (keV)
    input_multislice.theta = 0.0;                        % Till ilumination (�)
    input_multislice.phi = 0.0;                          % Till ilumination (�)

    %%%%%%%%%%%%%%%%%%%%%% Illumination model %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.illumination_model = 1;             % 1: coherente mode, 4: Numerical integration
    input_multislice.temporal_spatial_incoh = 1;         % 1: Temporal and Spatial, 2: Temporal, 3: Spatial

    %%%%%%%%%%%%%%%%%%%%%%%% condenser lens %%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.cond_lens_m = 0;                  % Vortex momentum
    input_multislice.cond_lens_c_30 = 0.001;            % Third order spherical aberration (mm)
    input_multislice.cond_lens_c_10 = ...
        ilc_scherzer_defocus(input_multislice.E_0,input_multislice.cond_lens_c_30);            % Defocus (�)

    input_multislice.cond_lens_c_50 = 0;                        % Fifth order spherical aberration (mm)
    % input_multislice.cond_lens_c_12 = 0.0;                    % Twofold astigmatism (�)
    % input_multislice.cond_lens_phi_12 = 0.0;                  % Azimuthal angle of the twofold astigmatism (�)
    % input_multislice.cond_lens_c_23 = 0.0;                    % Threefold astigmatism (�)
    % input_multislice.cond_lens_phi_23 = 0.0;                  % Azimuthal angle of the threefold astigmatism (�)
    input_multislice.cond_lens_inner_aper_ang = 0.0;            % Inner aperture (mrad) 
    input_multislice.cond_lens_outer_aper_ang = conv_angle;     % Outer aperture (mrad)


    %%%%%%%%% defocus spread function %%%%%%%%%%%%
%     dsf_sigma = ilc_iehwgd_2_sigma(32);                     % from defocus spread to standard deviation
%     input_multislice.cond_lens_ti_sigma = dsf_sigma;        % standard deviation (�)


    %%%%%%%%%% source spread function %%%%%%%%%%%%
    % ssf_sigma = il_fwhm_2_sigma(0.7);                  % half width at half maximum to standard deviation
%     input_multislice.cond_lens_si_sigma = 0.5;           % standard deviation: For parallel ilumination(�^-1); otherwise (�)

    %%%%%%%%% zero defocus reference %%%%%%%%%%%%
    input_multislice.cond_lens_zero_defocus_type = 2;   % eZDT_First = 1, eZDT_Middle = 2, eZDT_Last = 3, eZDT_User_Define = 4

    %%%%%%%%%%%%%%%%%%%%%%%%%%% Incident wave %%%%%%%%%%%%%%%%%%%%%%%%%%
    input_multislice.iw_type = 2;                                   % 1: Plane_Wave, 2: Convergent_wave, 3:User_Define, 4: auto

end
