function fse(kspace,varargin)
% fse sequence with arbitrary readout trajectory
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
arg.etl = 1; % echo train length (number of echoes per TR)

arg.slabthk = 2; % slice thickness (cm)

arg.rffa = 90; % tipdown flip angle (deg)
arg.rftbw = 8; % tipdown time-bandwidth product
arg.rfdur = 4.8; % tipdown pulse duration (ms)

arg.rf2fa = 180; % refocuser flip angle (deg)
arg.rf2tbw = 2; % refocuser time-bandwidth product
arg.rf2dur = 1.6; % refocuser pulse duration (ms)

arg.spoilcyc = 4; % number of phase cycles per voxel for gradient spoiling

arg.dofatsat = 1; % option to play fat sat pulse before echo train
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

% Create refocuser SLR pulse
[rf2,gzrf2,rffreq2] = toppe.utils.rf.makeslr( ...
    arg.rf2fa, ...
    arg.slabthk, ...
    arg.rf2tbw, ...
    arg.rf2dur, ...
    0, ...
    sys, ...
    'writeModFile', false ...
    );
toppe.writemod(sys, ...
    'rf', rf2, ...
    'gz', gzrf2, ...
    'nChop', [10,10], ...
    'ofname', 'refocuser.mod');

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
if mod(nshots,arg.etl)
    error('nshots must be divisible by etl')
end
ntrains = nshots/arg.etl;
ept = nshots/ntrains; % echos per train
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
dur_refocuser = length(toppe.readmod('refocuser.mod'))*sys.raster + arg.padtime; % us
dur_fatsat = length(toppe.readmod('fatsat.mod'))*sys.raster + arg.padtime; % us
dur_readout = length(toppe.readmod('readout.mod'))*sys.raster + arg.padtime; % us
dur_crusher = length(toppe.readmod('crusher.mod'))*sys.raster + arg.padtime; % us

% Determine effective echo delay in readout
if strcmpi(arg.ecdel,'auto')
    [~,arg.ecdel] = min(vecnorm(kspace(:,1,:),2,3));
    fprintf('spgr auto ecdel: sample delay of %d detected from first shot\n', arg.ecdel);
end

% Calculate minte1 = minimum time from peak tipdown to echo center
minte1 = dur_tipdown/2 + ... % tipdown pulse duration (post-peak)
    dur_refocuser + ... % refocuser pulse duration
    arg.padtime + ... % some extra fluffiness
    (nramp(1) + arg.ecdel)*sys.raster; % effective echo delay (wrt start of readout)
minte1 = minte1 + 1-mod(minte1,sys.raster);
minte1 = minte1*1e-3; % convert to ms

% Calculate minte2 = minimum time from echo center to next echo center
if arg.etl == 1
    minte2 = 0;
else
    minte2 = dur_refocuser + ... % refocuser pulse duration
        dur_readout + ... % readout length
        arg.padtime; % some extra fluffiness
end
minte2 = minte2 + 1-mod(minte2,sys.raster);
minte2 = minte2*1e-3; % convert to ms

% Set min te
minte = max(minte1,minte2);
if strcmpi(arg.te,'min')
    arg.te = minte;
    fprintf('fse auto min te: effective echo time = %.3fms\n', arg.te);
elseif (arg.te < minte)
    error('te (%.3fms) < minte (%.3fms)', arg.te, minte);
end

% Calculate gap times 1-3
tgap1 = arg.te/2 - dur_tipdown*1e-3/2 - dur_refocuser*1e-3/2; % time between end of tipdown & start of refocuser
tgap2 = arg.te/2 - dur_refocuser*1e-3/2 - (nramp(1) + arg.ecdel)*sys.raster*1e-3; % time between end of refocuser and beginning of readout
tgap3 = arg.te/2 - (nramp(2) + size(kspace,1)-arg.ecdel)*sys.raster*1e-3 - ...
    dur_refocuser*1e-3/2; % time between end of readout and beginning of next refocuser
fprintf('fse gap: gap time 1/2/3 = %.3f/%.3f/%.3fms\n', tgap1, tgap2, tgap3);

% Calculate minimum repitition time and tgap4
mintr = arg.dofatsat*dur_fatsat + ... % fatsat (us)
    arg.te*1e3/2 + ... % tipdown + gap time (us)
    arg.etl*arg.te*1e3; % refocuser + readout train (us)
mintr = mintr*1e-3; % convert to ms
if strcmpi(arg.tr,'min')
    arg.tr = mintr;
    fprintf('fse auto min tr: repitition time = %.3fms\n', arg.tr);
elseif (arg.tr < mintr)
    error('tr (%.3fms) < mintr (%.3fms)', arg.tr, mintr);
end
tgap4 = arg.tr - mintr;
fprintf('fse gap: gap time 4 = %.3fms\n', tgap4);

% Write entry file
toppe.writeentryfile('toppeN.entry');

% Set cores file entries
toppe.writecoresfile( {...
    [1,0], ... % tipdown
    [2,0], ... % refocuser
    [3,5], ... % fatsat
    [4,0], ... % readout
    [1,0] ... % TR deadtime
    })

% Set modules file text
modulesfiletext = [ ...
    sprintf("tipdown.mod %d 1 0 -1", dur_tipdown)
    sprintf("refocuser.mod %d 1 0 -1", dur_refocuser)
    sprintf("fatsat.mod %d 1 0 -1", dur_fatsat)
    sprintf("readout.mod %d 0 1 -1", dur_readout)
    sprintf("crusher.mod %d 0 0 -1", dur_crusher)
    ];

% Write out cores and modules files
writemodulesfile(modulesfiletext);

%% Write the scan loop
toppe.write2loop('setup',sys,'version',6);
for trainn = 1-arg.ndisdaqs:ntrains % loop through trains

    % Write fatsat to loop
    if arg.dofatsat
        toppe.write2loop('fatsat.mod', sys, ...
            'RFoffset', arg.fsdf, ...
            'RFphase', 0, ...
            'Gamplitude', [0;0;0], ...
            'version', 6, ...
            'trigout', 0, ...
            'core', 3);

        toppe.write2loop('crusher.mod', sys, ...
            'version', 6, ...
            'trigout', 0, ...
            'core', 3);
    end

    % Write tipdown to loop
    toppe.write2loop('tipdown.mod', sys, ...
        'RFoffset', rffreq, ...
        'RFphase', 0, ...
        'version', 6, ...
        'trigout', 0, ...
        'core', 1);

    % Write tgap1 to loop
    toppe.write2loop('delay', sys, ...
        'textra', tgap1, ...
        'core', 1);

    for echon = 1:ept % loop through echoes

        % Calculate shot index
        shotn = (trainn-1)*ept + echon;

        % Write refocuser to loop
        toppe.write2loop('refocuser.mod', sys, ...
            'RFoffset', rffreq2, ...
            'RFphase', pi/2, ...
            'version', 6, ...
            'trigout', 0, ...
            'core', 2);

        % Write tgap2 to loop
        toppe.write2loop('delay', sys, ...
            'textra', tgap2, ...
            'core', 2);

        % Write readout to loop
        if trainn > 0
            toppe.write2loop('readout.mod', sys, ...
                'RFphase', pi/2, ...
                'echo', 1, ...
                'slice', 1, ...
                'view', shotn, ...
                'waveform', shotn, ...
                'version', 6, ...
                'trigout', 0, ...
                'dabmode', 'on', ...
                'core', 4);
        else
            toppe.write2loop('readout.mod', sys, ...
                'version', 6, ...
                'trigout', 0, ...
                'dabmode', 'off', ...
                'core', 4);
        end

        % Write tgap3 to loop
        toppe.write2loop('delay', sys, ...
            'textra', tgap3, ...
            'core', 4);
    end

    % Write tgap4 to loop
    toppe.write2loop('delay', sys, ...
        'textra', tgap4, ...
        'core', 5);

end

% Finish loop
toppe.write2loop('finish',sys);
toppe.preflightcheck('toppeN.entry','seqstamp.txt',sys);

% Save arguments to file
save arg.mat arg
save kspace.mat kspace

end