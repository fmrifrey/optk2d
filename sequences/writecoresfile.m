function writecoresfile(cores)

fid = fopen('cores.txt', 'wt');
fprintf(fid, 'Total number of cores\n');
fprintf(fid, '%d\n', length(cores));
fprintf(fid, 'nmodules modIds... \n');

% Recall that module id's are defined by their row number in modules.txt
for icore = 1:length(cores)
    fprintf(fid, '%s\n', cores(icore));
end

fclose(fid);