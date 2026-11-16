function capture_water_processing(output_dir, spatial_points, points)
% Instrument a private copy; leave bundled MRSCloud source untouched.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'simulation','basissets','MRSCloud-main')));
mkdir(output_dir);
code=fileread(which('sim_signals'));
needle="outA = op_dccorr_FID_A(outA,'p');";
assert(numel(strfind(code,needle))==3); % UnEdited, MEGA, and multiplexed branches
code=strrep(code,'= sim_signals(MRS_opt)','= captured_sim_signals(MRS_opt)');
replacement="save('-mat7-binary',fullfile(save_dir,'before_first_dc.mat'),'outA'); outA = op_dccorr_FID_A(outA,'p');";
code=strrep(code,needle,replacement);
file=fullfile(output_dir,'captured_sim_signals.m');
fid=fopen(file,'w');fprintf(fid,'%s',code);fclose(fid);
addpath(output_dir);
p=struct('flipAngle',180,'centreFreq',4.7,'edit_flipAngle',180,'excflipAngle',90,...
    'nX',spatial_points,'nY',spatial_points,'nZ',spatial_points,...
    'localization',{{'PRESS'}},'FieldStr',{{'3T'}},'vendor',{{'Universal_Philips'}},...
    'seq',{{'UnEdited'}},'metab','H2O','Nmetab',1,'TEs',{{30}},'save_dir',output_dir);
p=load_parameters(p);p.Npts=points;
[p,after]=captured_sim_signals(p);
before=load(fullfile(output_dir,'before_first_dc.mat'));
raw_fid=before.outA.fids;
raw_spec=before.outA.specs;
ppm=before.outA.ppm;
after_first_dc=after.fids;
sw=p.sw;lw=p.lw;
save('-mat7-binary',fullfile(output_dir,'processing_trace.mat'),...
     'raw_fid','raw_spec','ppm','after_first_dc','sw','lw','spatial_points','points');
end
