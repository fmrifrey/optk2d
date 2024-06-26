function spgr(kspace,varargin)
% spgr sequence with arbitrary readout trajectory
% written in toppe v6 by David Frey
% kspace = kspace trajectory ([N x Nshots x 3] cm^-1)

import optk2d.utl.*

%% Set defaults
if nargin < 1 || isempty(kspace)
    kspace = zeros(500,1,3); % zeros - FID
end

arg.loadmat = 0; % option to load last used args and kspace from
                  % 'arg.mat' and 'kspace.mat', respectively

arg.ecdel = 'auto'; % eff echo time delay (# of samples from beginning of readout)
arg.te = 'min'; % echo time (ms)
arg.tr = 'min'; % repitition time (ms)

arg.slabthk = 2; % slice thickness (cm)

arg.rffa = 10; % tipdown flip angle (deg)
arg.rftbw = 8; % tipdown time-bandwidth product
arg.rfdur = 3.2; % tipdown pulse duration (ms)
arg.rfspoil = 1; % option to do rf spoiling
arg.spoilcyc = 4; % number of phase cycles per voxel for gradient spoiling

arg.dofatsat = 1; % option to dplay fat sat pulse before readout
arg.fsfa = 90; % fatsat flip angle (deg)
arg.fstbw = 2; % fatsat time-bandwidth product
arg.fsdur = 2.4; % fatsat pulse duration (ms)
arg.fsdf = -440; % fatsat frequency offset (Hz)

arg.padtime = 500; % extra deadtime in each core (us)
arg.ndisdaqs = 0; % number of disdaq shots at beginning of sequence

% load in variable input arguments
if arg.loadmat
    load arg.mat arg
    load arg.kspace kspace
end

% Parse variable inputs (will overwrite loaded arguments!)
arg = toppe.utils.vararg_pair(arg, varargin);

%% Generate pulse sequence modules
sys = toppe.systemspecs;

% Create tipdown SLR pulse
[rf,gzrf,rffreq] = toppe.utils.rf.makeslr( ...
    arg.rffa, ...
    arg.slabthk, ...
    arg.rftbw, ...
    arg.rfdur, ...
    eps(), ...
    sys, ...
    'writeModFile', false ...
    );
toppe.writemod(sys, ...
    'rf', rf, ...
    'gz', gzrf, ...
    'nChop', [10,10], ...
    'ofname', 'tipdown.mod');

% Create the fatsat SLR pulse
fs = toppe.utils.rf.makeslr( ...
    arg.fsfa, ...
    5, ... % dummy variable - we won't use the slicesel gradient
    arg.fstbw, ...
    arg.fsdur, ...
    0, ...
    sys, ...
    'writeModFile', false ...
    );
toppe.writemod(sys, ...
    'rf', fs, ...
    'nChop', [10,10], ...
    'ofname', 'fatsat.mod');

% Generate the readout gradients
nshots = size(kspace,2);
[g,nramp] = gengrads(sys,kspace);
gx = g(:,:,1);
gy = g(:,:,2);
gz = g(:,:,3);
toppe.writemod(sys, ...
    'gx', gx, ...
    'gy', gy, ...
    'gz', gz, ...
    'nChop', nramp, ...
    'ofname', 'readout.mod');

% Create a crusher gradient
gspoil = toppe.utils.makecrusher(arg.spoilcyc, ...
    arg.slabthk, ...
    sys, ...
    0, ...
    sys.maxSlew, ...
    sys.maxGrad);
toppe.writemod(sys, ...
    'gz', gspoil, ...
    'ofname', 'crusher.mod');

% Calculate all module durations
dur_tipdown = length(toppe.readmod('tipdown.mod'))*sys.raster + arg.padtime; % us
dur_fatsat = length(toppe.readmod('fatsat.mod'))*sys.raster + arg.padtime; % us
dur_readout = length(toppe.readmod('readout.mod'))*sys.raster + arg.padtime; % us
dur_crusher = length(toppe.readmod('crusher.mod'))*sys.raster + arg.padtime; % us

% Determine effective echo delay in readout
if strcmpi(arg.ecdel,'auto')
    [~,arg.ecdel] = min(vecnorm(kspace(:,1,:),2,3));
    fprintf('spgr auto ecdel: sample delay of %d detected from first shot\n', arg.ecdel);
end

% Calculate minimum echo time and tgap1
[~,peakidx] = max(abs(rf)); % get rf peak sample
minte = dur_tipdown - peakidx*sys.raster + ... % rf pulse duration (post-peak)
    arg.padtime + ... % some extra fluffiness
    (nramp(1) + arg.ecdel)*sys.raster; % effective echo delay (wrt start of readout)
minte = minte*1e-3; % convert to ms
if strcmpi(arg.te,'min')
    arg.te = minte;
    fprintf('spgr auto min te: effective echo time = %.3fms\n', arg.te);
elseif (arg.te < minte)
    error('te (%.3fms) < minte (%.3fms)', arg.te, minte);
end
tgap1 = arg.te - minte;
fprintf('spgr gap: gap time 1 = %.3fms\n', tgap1);

% Calculate minimum repitition time and tgap2
mintr = arg.dofatsat*dur_fatsat + ... % fatsat (us)
    dur_tipdown + ... % tipdown (us)
    tgap1*1e3 + ... % gap 1 (us)
    dur_readout + ... % readout (us)
    dur_crusher; % crusher (us)
mintr = mintr*1e-3; % convert to ms
if strcmpi(arg.tr,'min')
    arg.tr = mintr;
    fprintf('spgr auto tr: repitition time = %.3fms\n', arg.tr);
elseif (arg.tr < mintr)
    error('tr (%.3fms) < mintr (%.3fms)', arg.tr, mintr);
end
tgap2 = arg.tr - mintr;
fprintf('spgr gap: gap time 2 = %.3fms\n', tgap2);

% Write entry file
toppe.writeentryfile('toppeN.entry');

% Set cores file entries
toppe.writecoresfile( {...
    [4,1,0], ... % tipdown
    2, ... % fatsat
    [3,0] ... % readout
    })

% Set modules file text
modulesfiletext = [ ...
    sprintf("tipdown.mod %d 1 0 -1", dur_tipdown)
    sprintf("fatsat.mod %d 1 0 -1", dur_fatsat)
    sprintf("readout.mod %d 0 1 -1", dur_readout)
    sprintf("crusher.mod %d 0 0 -1", dur_crusher)
    ];

% Write out cores and modules files
writemodulesfile(modulesfiletext);

%% Write the scan loop
toppe.write2loop('setup',sys,'version',6);
for shotn = 1-arg.ndisdaqs:nshots % loop through shots

    % Write fatsat to loop
    if arg.dofatsat
        toppe.write2loop('fatsat.mod', sys, ...
            'RFoffset', arg.fsdf, ...
            'RFphase', 0, ...
            'Gamplitude', [0;0;0], ...
            'version', 6, ...
            'trigout', 0, ...
            'core', 2);
    end

    % Write tipdown to loop
    toppe.write2loop('crusher.mod', sys, ...
        'version', 6, ...
        'trigout', 0, ...
        'core', 1);
    toppe.write2loop('tipdown.mod', sys, ...
        'RFspoil', arg.rfspoil, ...
        'RFoffset', rffreq, ...
        'RFphase', 0, ...
        'version', 6, ...
        'trigout', 0, ...
        'core', 1);

    % Write tgap1 to loop
    toppe.write2loop('delay', sys, ...
        'textra', tgap1, ...
        'core', 1);

    % Write readout to loop
    if shotn > 0
        toppe.write2loop('readout.mod', sys, ...
            'RFspoil', arg.rfspoil, ...
            'echo', 1, ...
            'slice', 1, ...
            'view', shotn, ...
            'waveform', shotn, ...
            'RFphase', 0, ...
            'version', 6, ...
            'trigout', 0, ...
            'dabmode', 'on', ...
            'core', 3);
    else
        toppe.write2loop('readout.mod', sys, ...
            'RFspoil', arg.rfspoil, ...
            'version', 6, ...
            'RFphase', 0, ...
            'trigout', 0, ...
            'dabmode', 'off', ...
            'core', 3);
    end

    % Write tgap2 to loop
    toppe.write2loop('delay', sys, ...
        'textra', tgap2, ...
        'core', 3);

end

% Finish loop
toppe.write2loop('finish',sys);
toppe.preflightcheck('toppeN.entry','seqstamp.txt',sys);

% Save arguments to file
save arg.mat arg
save kspace.mat kspace

end