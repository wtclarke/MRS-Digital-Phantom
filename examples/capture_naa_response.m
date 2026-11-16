function capture_naa_response(output_dir, spatial_points, points)
% Capture spatially integrated pre-DC responses using a private source copy.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'simulation','basissets','MRSCloud-main')));
mkdir(output_dir);
code=fileread(which('sim_signals'));
needle="outA = op_dccorr_FID_A(outA,'p');";
assert(numel(strfind(code,needle))==3);
code=strrep(code,'= sim_signals(MRS_opt)','= naa_capture_sim_signals(MRS_opt)');
replacement="save('-mat7-binary',fullfile(save_dir,'before_first_dc.mat'),'outA'); outA = op_dccorr_FID_A(outA,'p');";
code=strrep(code,needle,replacement);
fid=fopen(fullfile(output_dir,'naa_capture_sim_signals.m'),'w');
fprintf(fid,'%s',code);fclose(fid);addpath(output_dir);
metabs={'NAA','NAA_Acetyl_only','NAA_Acetyl_only','NAA_Acetyl_only','NAA_Aspartate_only','NAA_Acetyl_only'};
labels={'full','methyl','methyl_water_carrier','methyl_on_carrier','aspartyl','methyl_repeat'};
carriers=[3.,3.,4.7,2.008,3.,3.];
raw=struct();first_dc=struct();raw_spec=struct();ppm=struct();configs=struct();
for k=1:length(labels)
    folder=fullfile(output_dir,labels{k});mkdir(folder);
    p=struct('flipAngle',180,'centreFreq',carriers(k),'edit_flipAngle',180,'excflipAngle',90,...
        'nX',spatial_points,'nY',spatial_points,'nZ',spatial_points,...
        'localization',{{'PRESS'}},'FieldStr',{{'3T'}},'vendor',{{'Universal_Philips'}},...
        'seq',{{'UnEdited'}},'metab',metabs{k},'Nmetab',1,'TEs',{{30}},'save_dir',folder);
    p=load_parameters(p);p.Npts=points;
    configs.(labels{k})=struct('metabolite',metabs{k},'carrier_ppm',carriers(k),...
        'spin_system_carrier_relative',p.sys,'TE_ms',30,'TE1_ms',p.TE1,'TE2_ms',p.TE2,...
        'refocusing_duration_ms',p.refTp,'refocusing_waveform',p.refRF,...
        'bandwidth_hz',p.sw,'points',points,'spatial_points',spatial_points,...
        'readout_linewidth_hz',p.lw,'B0_T',p.Bfield);
    [p,after]=naa_capture_sim_signals(p);
    before=load(fullfile(folder,'before_first_dc.mat'));
    raw.(labels{k})=before.outA.fids;
    raw_spec.(labels{k})=before.outA.specs;
    ppm.(labels{k})=before.outA.ppm;
    first_dc.(labels{k})=after.fids;
end
save('-mat7-binary',fullfile(output_dir,'naa_trace.mat'),...
    'raw','first_dc','raw_spec','ppm','configs','carriers','labels');
end
