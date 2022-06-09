function generate_training_data_64(db_n, gpu_n, b_val, target_dir)
    clc;
    
    % activate / deactivate tests and plotting stuff to check each function
    global b_test
    b_test = False;
    
    if nargin < 4
        target_dir = ['.' filesep 'Data_64'];
    end
    if nargin < 3
        b_val = 0;
    end
    if nargin < 2
        gpu_n = 0;
    end
    if nargin < 1
        db_n = 0;
    end
    
    
    s.gpu = gpu_n+1;
    warning('off','all')
    
    addpath('functions');
    addpath(['atomic_specimen_creation' filesep 'src']);

    s.n_batch = 8;
    s.NUM_DATA = [8192 32 inf]; 
    s.DS = {'Validation', 'Test', 'Training'};
    
    %Simulation Parametes
    s.E0 = [30 40 50 60 80 100 120 140 160 180 200 300];        % keV
    s.SCAN_STEPS = 3;                   % points along each dimension
    s.OBJ_SIZE = 80;                    % Effective simulation box size
    s.THICKNESS = [3 30];               % Sample Thickness range
    s.CONV_ANGLE = [5 30];             % mrad

    % Detector parameters
    s.NX_REAL = 64;                     % Real space pixels
    s.NX_DET = 64;                      % Detector pixels
    s.NX_POT = 1440;                    % Potential sampling

    % Specimen information
    cif_folder = ['cif_files' filesep 'pymatgen'];
    cif_files = dir(cif_folder); 
    b_hid = startsWith({cif_files.name}, '.');
    b_dir = [cif_files.isdir];
    b_cif = contains({cif_files.name}, '.cif');
    cif_files(b_dir | ~b_cif | b_hid) = [];
    rm_fields = {'folder','date','bytes','datenum','isdir'};
    cif_files = rmfield(cif_files,rm_fields);
    clear b_dir b_cif b_hid
    
    % Log
    header = {"mp-id","Composition","Space Group","Zone Axis","Rotation","Thickness","E0","Aperture","Detector Size","Step Size"}; %#ok<CLARRSTR>
    
    T_hkls = [  1 1 0;  0 1 1; 1 0 1; ...
                0 0 1; 1 0 0; 0 1 0; 1 1 1]; 
    
    s.FILENAME = ['db_h5_b_' num2str(db_n)];

    % load standard parameters for Multem input
    input_multislice = get_multem_parameters(s.NX_POT, s.CONV_ANGLE(1), s.E0(1), gpu_n);

    features = zeros(s.NX_DET,s.NX_DET,s.SCAN_STEPS^2,s.n_batch,'uint16');
    labels_k = zeros(s.NX_DET,s.NX_DET, 2, s.n_batch,'uint16');
    labels_r = zeros(s.NX_REAL,s.NX_REAL, 2, s.n_batch,'uint16');
    probe_r = zeros(s.NX_REAL,s.NX_REAL, 2, s.n_batch,'uint16');
    meta = zeros(s.SCAN_STEPS^2 + 10 ,s.n_batch,'single');
    log_prms = cell(s.n_batch,10);
    
    % Go!
    if b_val; ds1 = 1; else ds1 = 3; end
    for i_dat = ds1:numel(s.DS)
        s.FOLDER = [target_dir filesep s.DS{i_dat}];
        hdf_file = [s.FOLDER filesep s.FILENAME '_' s.DS{i_dat} '.h5']; 
        [~,~] = mkdir(s.FOLDER);
        [idx_s, ~] = get_or_make_h5(hdf_file, s);
        log_path = [s.FOLDER filesep 'data_log_' num2str(db_n) '_' s.DS{i_dat} '.xlsx'];
        writecell(header,log_path,'FileType','spreadsheet','Sheet',1,'Range','A1')
        
        tic;
        msg = fprintf(' ');
        for idx = idx_s:s.n_batch:s.NUM_DATA(i_dat) 
            n_batch = min([s.n_batch,(s.NUM_DATA(i_dat) - idx)]);
            for b = 1:n_batch
               [features(:,:,:,b), labels_k(:,:,:,b), labels_r(:,:,:,b), probe_r(:,:,:,b), meta(:,b), spec] = run(input_multislice, cif_folder, cif_files, T_hkls, s);
                log_prms(b,1:10) = {spec.trial_mat, spec.trial_mat_par.formula, spec.trial_mat_par.sgn, num2str(spec.trial_ori), spec.trial_rot, spec.trial_thick, spec.E_0, spec.alpha, spec.gmax, spec.step_size};
            end
            
            % Write to log
            writecell(log_prms(1:n_batch,:),log_path,'FileType','spreadsheet','Sheet',1,'Range',['A' num2str(idx+2)])
            
            % Write to db
            h5write(hdf_file, '/features', features,[1 1 1 idx+1],[s.NX_DET,s.NX_DET,s.SCAN_STEPS^2,s.n_batch]);
            h5write(hdf_file, '/labels_k', labels_k, [1 1 1 idx+1],[s.NX_DET,s.NX_DET, 2, s.n_batch]);
            h5write(hdf_file, '/labels_r', labels_r, [1 1 1 idx+1],[s.NX_REAL,s.NX_REAL, 2, s.n_batch]);
            h5write(hdf_file, '/probe_r', probe_r, [1 1 1 idx+1],[s.NX_REAL,s.NX_REAL, 2, s.n_batch]);
            h5write(hdf_file, '/meta', meta, [1, idx+1],[s.SCAN_STEPS^2+10,s.n_batch]);
            post_batch_idx = idx  +s.n_batch;
            h5writeatt(hdf_file,'/','idx', post_batch_idx);
            h5writeatt(hdf_file,'/','State', rng().State);           
            
            % Print progress to CLI
            fprintf(repmat('\b', 1, msg));
            str = [s.FILENAME , '   Progress: ' num2str(post_batch_idx) ,' / ' num2str(s.NUM_DATA(i_dat)) '   avg time =' num2str(toc/(post_batch_idx-idx_s),3) 's per sample \n'];
            msg = fprintf(str);
        end    
    end
end

%% Functions

function [feature, label_k, label_r, probe, meta, spec] = run(input_multislice, cif_folder, cif_files, T_hkls, s)
    % Run Simulation
    [X, Y, inp, nx, g_max, spec] = rnd_prms(input_multislice, cif_folder, cif_files, T_hkls,s);
    [feature, label_k, label_r, probe, max_sc] = fcn_run_simulation_over_grid(inp, X, Y, nx, s);
    spec.E_0 = inp.E_0;
    spec.alpha = inp.cond_lens_outer_aper_ang;
    spec.gmax = g_max;
    meta = single([inp.E_0, inp.cond_lens_outer_aper_ang, g_max, spec.step_size, max_sc]');

end

function [idx, rng_gen] = get_or_make_h5(hdf_file, s)
    if ~exist(hdf_file, 'file')
        rng_gen = fcn_new_rng();

        h5create(hdf_file,'/features',[s.NX_DET,s.NX_DET,s.SCAN_STEPS^2, inf],'ChunkSize',[s.NX_DET,s.NX_DET,s.SCAN_STEPS^2,1],'Datatype','uint16','Deflate',9)
        h5create(hdf_file,'/labels_k',[s.NX_DET,s.NX_DET, 2, inf],'ChunkSize',[s.NX_DET,s.NX_DET, 2, 1],'Datatype','uint16','Deflate',9)
        h5create(hdf_file,'/labels_r',[s.NX_REAL,s.NX_REAL, 2, inf],'ChunkSize',[s.NX_REAL,s.NX_REAL, 2, 1],'Datatype','uint16','Deflate',9)
        h5create(hdf_file,'/probe_r',[s.NX_REAL,s.NX_REAL, 2, inf],'ChunkSize',[s.NX_REAL,s.NX_REAL, 2, 1],'Datatype','uint16','Deflate',9)
        h5create(hdf_file,'/meta',[s.SCAN_STEPS^2+10, inf],'ChunkSize',[s.SCAN_STEPS^2+10, 1],'Datatype','single','Deflate',9)

        idx = 0;
        h5writeatt(hdf_file,'/','Seed',rng_gen.Seed)
        h5writeatt(hdf_file,'/','State',rng_gen.State)
        h5writeatt(hdf_file,'/','Type',rng_gen.Type)
        h5writeatt(hdf_file,'/','idx',idx)
        h5writeatt(hdf_file,'/','arch',computer('arch'))
        h5writeatt(hdf_file,'/','gpu',gpuDevice(s.gpu).Name)
        h5writeatt(hdf_file,'/','matlab_ver',version('-release'))
    else
        idx = h5readatt(hdf_file,'/','idx');
        rng_struct.Seed = h5readatt(hdf_file,'/','Seed');
        rng_struct.Type = h5readatt(hdf_file,'/','Type');
        rng_struct.State = h5readatt(hdf_file,'/','State');
        rng(rng_struct);
        rng_gen = rng();
    end    
end


function [X,Y, input_multislice, nx, g_max, spec] = rnd_prms(input_multislice, cif_folder, cif_files, T_hkls, s)

    while 1
        
        % Convergence / Collection angle ratio
        bf_r = tfm_rnd_unif_lim(0.1,0.95);
        
        % Random Energy
        rnd_e0 = tfm_rnd_norm_lim(s.E0(1),s.E0(end)+30,100,180);
        while 220 < rnd_e0 && rnd_e0 < 280
            rnd_e0 = tfm_rnd_norm_lim(s.E0(1),s.E0(end)+30,100,180);
        end
        [~,e0_id] = min(abs(s.E0-rnd_e0));
        e0 = s.E0(e0_id);
        
        while 1
            % Random Probe size
            ap = tfm_rnd_norm_lim(s.CONV_ANGLE(1),s.CONV_ANGLE(2),10,20);

            % Collection angle and g-max follows from Energy and Convergence / Collection angle ratio
            col_ang = ap / bf_r;
            g_max = ilm_mrad_2_rAng(e0,col_ang);
            nx = ceil(g_max*2*s.OBJ_SIZE);
            if nx >= s.NX_DET
                break;
            end
        end
        
        % Quasi Probe size in real space
        probe_s_r = fcn_first_zero_radius(e0, ap);
        
        % Step size (Ensure some probe overlap)
        spec.step_size = tfm_rnd_unif_lim(probe_s_r*0.1,probe_s_r*0.9);

        % Set scan grid and simulation parameters
        [X, Y] = fcn_scan_points(spec.step_size);
        input_multislice.E_0 = e0;
        input_multislice.cond_lens_outer_aper_ang = ap;
        
        % Specimen   
        mat_x = ceil(rand(1,1)*numel(cif_files));
        spec.trial_mat = cif_files(mat_x).name;
        spec.trial_mat_par = tfm_get_uc_from_cif([cif_folder filesep spec.trial_mat]);
        spec.trial_ori = T_hkls(ceil(rand()*length(T_hkls)),:);
        spec.trial_thick = tfm_rnd_unif_lim(s.THICKNESS(1),s.THICKNESS(2));
        spec.trial_rot = rand(1,1) * 360;
        
        input_multislice = fcn_gen_specimen(spec.trial_mat_par, spec.trial_ori, spec.trial_rot, spec.trial_thick, s.OBJ_SIZE, input_multislice);
        
        % Ensure specimen creation was successful
        if  ~isempty(input_multislice.spec_atoms)
            break;
        end
    end
end                                         

function [X, Y] = fcn_scan_points(step_size)
    % Returns X and Y coordinates for a grid of scan points around [0, 0]
    sr2 = (step_size*3)/2;
    pts = linspace(-sr2,sr2,3);
    [X, Y] = meshgrid(pts,pts);

end

function input_multislice =  fcn_gen_specimen(trial_mat_par, trial_ori, trial_rot, trial_thick, lx, input_multislice)
    global b_test
    if b_test
        figure(1); clf;
        h = gca;
    else
        h = [];
    end
    % Generate Supercell
    trial_mat_par.atoms = tfm_align_duplicate_cut(trial_mat_par, trial_ori, trial_rot, lx-10, lx-10, trial_thick, b_test, false, h);         
    
    if ~isempty(trial_mat_par.atoms)
        % Center specimen
        trial_mat_par.atoms = ilm_center_spec(trial_mat_par.atoms,  lx,  lx, trial_thick);

        % Convert Supercell to Multem spec_atoms input
        trial_mat_par.atoms(:,5) = 0.08;
        trial_mat_par.atoms(:,6) = 1;
        input_multislice.spec_atoms = double(gather(trial_mat_par.atoms));
        input_multislice.spec_lx = lx;
        input_multislice.spec_ly = lx;
        input_multislice.spec_lz = trial_thick;

        if b_test
           [~, slices] = ilc_spec_slicing(input_multislice.toStruct);
           [nslice, ~] = size(slices);
           V = zeros(input_multislice.ny,input_multislice.nx);
           for ix=1:nslice
               input_multislice.islice = ix;
               pr = input_multislice.ilc_projected_potential;
               V = V + pr.V;
           end
           figure(2); clf; imagesc(pr.x,pr.y,V);
           set(gca,'YDir','normal')
           axis image
        end
    end
end

function [x, x_max] = parse_uint16(x)
    x_max = max(x(:));
    x = uint16(x./x_max*65536);
end

function [A, P, Max_A, Max_P] = fcn_wave_parts(wave, nx, s ,domain)
    if strcmp(domain, 'real')
        nxc = ceil(nx/2) * 2 + 1;
        val_1 = sum_wave(wave);
        wave(end+1,:) = 0;
        wave(:,end+1) = 0;
        wave = imresize(wave, [nxc,nxc],'bilinear');
        val_2 = sum_wave(wave);
        wave = wave*val_1/val_2/s.NX_POT;
        n_crop = ceil((nxc - s.NX_REAL)/2);
        wave = wave(n_crop:end-n_crop,n_crop:end-n_crop);
    elseif strcmp(domain, 'fourier')
        nxc = ceil(nx/2) * 2;
        n_crop = round((s.NX_POT - nxc)/2);
        wave = wave(n_crop+1:end-n_crop,n_crop+1:end-n_crop);
%         nxc = nxc + 1;
        wave(end+1,:) = 0;
        wave(:,end+1) = 0;
        val_1 = sum_wave(wave);
        wave = imresize(wave, [s.NX_DET+1 s.NX_DET+1],'bilinear');
        wave(end,:) = [];
        wave(:,end) = [];
        val_2 = sum_wave(wave);
        wave = wave*val_1/val_2;
    end

    if isreal(wave)
        [A, Max_A] = parse_uint16(wave);
        P = []; Max_P = [];
    else
        [A, Max_A] = parse_uint16(abs(wave));
        [P, Max_P] = parse_uint16(angle(wave)+pi);
    end
end

function [feature, label_k, label_r, probe_r, max_sc] = fcn_run_simulation_over_grid(input_multislice, X, Y, nx, s)
    global b_test
    if b_test
        n_pat = sqrt(numel(Y(:)));
        plot_order = reshape(fliplr(reshape(1:n_pat^2,n_pat,n_pat))',n_pat^2,1);
        figure(3); clf();
        tiledlayout(n_pat,n_pat,'TileSpacing','none');
    end

    feature = zeros(s.NX_DET, s.NX_DET, numel(Y(:)),'uint16');
    label_k = zeros(s.NX_DET, s.NX_DET, 2, 'uint16');
    label_r = zeros(s.NX_REAL, s.NX_REAL, 2, 'uint16');
    probe_r = zeros(s.NX_DET, s.NX_DET, 2, 'uint16');
    max_sc = zeros(1,numel(Y(:))+6, 'single');

    for py = 1:numel(Y(:))
        b_center = py == (numel(Y(:))+1)/2;
        input_multislice.iw_x = X(py) + input_multislice.spec_lx/2;      % x position 
        input_multislice.iw_y = Y(py) + input_multislice.spec_ly/2;      % y position
        input_multislice.pn_coh_contrib = b_center;                      % Output also coherent wave for central cbed
        
        wave_r = input_multislice.ilc_multem;
        [cbed_r, ~, Max_I, ~] = fcn_wave_parts(wave_r.data.m2psi_tot, nx, s, 'fourier');
        feature(:,:,py) = cbed_r;
        max_sc(py) = Max_I;
        
        if b_center
            wave_real = fftshift(ifft2(ifftshift(wave_r.data.psi_coh)))*numel(wave_r.data.psi_coh);
            probe = input_multislice.ilc_incident_wave.psi_0;
            
            [I_k, phase_k, Max_I_k, Max_P_k] = fcn_wave_parts(wave_r.data.psi_coh, nx, s, 'fourier');
            [I_r, phase_r, Max_I_r, Max_P_r] = fcn_wave_parts(wave_real, nx, s, 'real');
            [P_I_r, P_phase_r, P_Max_I_r, P_Max_P_r] = fcn_wave_parts(probe, nx, s, 'real');
            
            label_k(:,:,1) = phase_k;
            label_k(:,:,2) = I_k;
            label_r(:,:,1) = phase_r;
            label_r(:,:,2) = I_r;
            probe_r(:,:,1) = P_phase_r;
            probe_r(:,:,2) = P_I_r;
            max_sc(end-5:end) = [Max_P_k Max_I_k Max_P_r Max_I_r, P_Max_P_r, P_Max_I_r];
        end

        if b_test
            figure(3);
            nexttile(plot_order(py))  
            plot_im = double(cbed_r) * Max_I / 65536;
            imagesc(plot_im); set(gca,'ColorScale','log'); axis image off;
            drawnow()
        end
    end

    if b_test  
        figure(4); clf();
        tiledlayout(3,2);
        nexttile
        imagesc(double(label_k(:,:,1)) * max_sc(end-5) / 65536); axis image off;
        nexttile
        imagesc(double(label_k(:,:,2)) * max_sc(end-4) / 65536); set(gca,'ColorScale','log'); axis image off;
        nexttile
        imagesc(double(label_r(:,:,1)) * max_sc(end-3) / 65536); axis image off;
        nexttile
        imagesc(double(label_r(:,:,2)) * max_sc(end-2) / 65536); set(gca,'ColorScale','log'); axis image off;
        
        nexttile
        p = (double(probe_r(:,:,2))*P_Max_I_r / 65536) .* exp(1i*double(probe_r(:,:,1))*P_Max_P_r / 65536 - pi);
        ew = (double(label_r(:,:,2))* max_sc(end-2) / 65536) .* exp(1i*double(label_r(:,:,1)) * max_sc(end-3) / 65536 - pi);
        
        imagesc(angle(ew./p)); axis image off;

    end
end

function val = sum_wave(wave)
    if isreal(wave)
        val = sum(sum(wave));
    else
        val = sqrt(sum(sum(abs(wave).^2)));
    end
end

function r = fcn_first_zero_radius(E0, conv_angle)
    emass = 510.99906;      % electron rest mass in keV
    hc = 12.3984244;        % Planck's const x speed of light
    lambda = hc/sqrt(E0*(2.0*emass + E0)); % A
    rAng = 2 * conv_angle / lambda*1e-3;
    % 3.8317 is the fisrt zero of the bessel function of the first type
    r = 3.8317/rAng/pi;
end