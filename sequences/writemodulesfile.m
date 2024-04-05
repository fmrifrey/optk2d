function writemodulesfile(mods)

fid = fopen("modules.txt","w");
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', length(mods));
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ? hastrigout?\n');

for imod = 1:length(mods)
    fprintf(fid, '%s\n', mods(imod));
end

fclose(fid);