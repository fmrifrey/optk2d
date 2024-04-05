function spgr(kspace,varargin)
% spgr sequence with arbitrary readout trajectory
% written in toppe v6 by David Frey
% kspace = kspace trajectory ([N x Nshots x 3] cm^-1)

%% Set defaults
if nargin < 1 || isempty(kspace)
    kspace = zeros(500,1,3); % zeros - FID
end
arg.fov = 22; % XY FOV (cm)
arg.N = 320; % image matrix size
arg.slthick = 0.2; % slice thickness (cm)
arg.slspace = 0; % slice spacing (cm)
arg.nslices = 1; % number of slices to acquire
arg.nprephs = 500; % number of pre-phaser points to add to beginning of traj
arg.nrewind = 500; % number of rewinder points to add to end of traj
arg.dofatsat = 1; % option to do play fat sat pulse before readout
arg.rffa = 90; % tipdown flip angle (deg)
arg.rftbw = 8; % tipdown time-bandwidth product
arg.rfdur = 3; % tipdown pulse duration (ms)
arg.fsfa = 90; % fatsat flip angle (deg)
arg.fstbw = 2; % fatsat time-bandwidth product
arg.fsdur = 1; % fatsat pulse duration (ms)
arg.fsdf = -440; % fatsat frequency offset (Hz)
arg.rfspoil = 1; % option to do rf spoiling
arg.spoilcyc = 12; % number of phase cycles per voxel for gradient spoiling
arg.padtime = 100; % extra deadtime in each core (us)
arg.deadtime = 5; % extra deadtime at the end of each shot (ms)
arg.ndisdaqs = 0; % number of disdaq shots at beginning of sequence

% Parse variable inputs
arg = toppe.utils.vararg_pair(arg, varargin);

%% Generate pulse sequence modules
sys = toppe.systemspecs('maxRF',0.3);

% Create a gradient spoiler
spoilgz = toppe.utils.makecrusher(arg.spoilcyc, ...
    arg.nslices*(arg.slthick + arg.slspace), ...
    sys, ...
    0, ...
    sys.maxSlew, ...
    sys.maxGrad);

% Create tipdown SLR pulse
[rf, rfgz, rffreq] = toppe.utils.rf.makeslr( ...
    arg.rffa, ...
    arg.slthick, ...
    arg.rftbw, ...
    arg.rfdur, ...
    0, ... % no spoil, we'll add the crusher manually
    sys, ...
    'sliceOffset', arg.nslices*(arg.slthick + arg.slspace), ...
    'writeModFile', false ...
    );

% Append the crusher to tipdown
rf = [zeros(size(spoilgz)); rf];
rfgz = [spoilgz; rfgz];

% Write tipdown out to mod file
toppe.writemod(sys, ...
    'rf', rf, ...
    'gz', rfgz, ...
    'ofname', 'tipdown.mod');

% Create the fatsat SLR pulse
[fs,fsgz] = toppe.utils.rf.makeslr( ...
    arg.fsfa, ...
    arg.nslices*(arg.slthick + arg.slspace), ...
    arg.fstbw, ...
    arg.fsdur, ...
    0, ...
    sys, ...
    'ftype', 'ls', ...
    'type', 'ex', ...
    'writeModFile', false ...
    );
fsgz = 0*fsgz; % zero out slice selective gradient

% Append the crusher to fatsat
fs = [fs; zeros(size(spoilgz))];
fsgz = [fsgz; spoilgz];

% Write fatsat out to mod file
toppe.writemod(sys, ...
    'rf', fs, ...
    'gz', fsgz, ...
    'ofname', 'fatsat.mod');

% Generate the readout gradients
nshots = size(kspace,2);
g = gengrads(sys,kspace);
gx = g(:,:,1);
gy = g(:,:,2);
gz = g(:,:,3);
toppe.writemod(sys, ...
    'gx', gx, ...
    'gy', gy, ...
    'gz', gz, ...
    'ofname', 'readout.mod');

% Calculate all module durations
dur_tipdown = length(toppe.readmod('tipdown.mod'))*sys.raster + arg.padtime; % us
dur_fatsat = length(toppe.readmod('fatsat.mod'))*sys.raster + arg.padtime; % us
dur_readout = length(toppe.readmod('readout.mod'))*sys.raster + arg.padtime; % us

% Write entry file
toppe.writeentryfile('toppeN.entry');

% Set cores file text
coresfiletext = [ ...
    "1 1" % tipdown
    "1 2" % fatsat
    "1 3" % readout
    ];

% Set modules file text
modulesfiletext = [ ...
    sprintf("tipdown.mod %d 1 0 -1", dur_tipdown)
    sprintf("fatsat.mod %d 1 0 -1", dur_fatsat)
    sprintf("readout.mod %d 0 1 -1", dur_readout)
    ];

% Write out cores and modules files
writemodulesfile(modulesfiletext);
writecoresfile(coresfiletext);

%% Write the scan loop
toppe.write2loop('setup',sys,'version',6);
sloffset = round(linspace(-rffreq,rffreq,arg.nslices));
for shotn = 1-arg.ndisdaqs:nshots % loop through shots
    for slicen = 1:arg.nslices % loop through slices

        % Write fatsat to loop
        if arg.dofatsat
            toppe.write2loop('fatsat.mod', sys, ...
                'RFoffset', arg.fsdf, ...
                'RFphase', 0, ...
                'RFspoil', 0, ...
                'version', 6, ...
                'trigout', 0, ...
                'core', 2);
        end

        % Write tipdown to loop
        toppe.write2loop('tipdown.mod', sys, ...
            'RFoffset', sloffset(slicen), ...
            'RFphase', 0, ...
            'RFspoil', arg.rfspoil, ...
            'version', 6, ...
            'trigout', 0, ...
            'core', 1);

        % Write readout to loop
        if shotn > 0
            toppe.write2loop('readout.mod', sys, ...
                'echo', 1, ...
                'slice', slicen, ...
                'view', shotn, ...
                'waveform', shotn, ...
                'RFspoil', arg.rfspoil, ...
                'version', 6, ...
                'trigout', 0, ...
                'dabmode', 'on', ...
                'textra', arg.deadtime, ...
                'core', 3);
        else
            toppe.write2loop('readout.mod', sys, ...
                'RFspoil', arg.rfspoil, ...
                'version', 6, ...
                'trigout', 0, ...
                'dabmode', 'off', ...
                'textra', arg.deadtime, ...
                'core', 3);
        end

    end
end

% Finish loop
toppe.write2loop('finish',sys);
toppe.preflightcheck('toppeN.entry','seqstamp.txt',sys);

end