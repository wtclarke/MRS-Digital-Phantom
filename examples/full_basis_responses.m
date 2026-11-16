function full_basis_responses(output_dir, entries_file, spatial_points, points, signature)
% Independent experiment; per-response caches are bound to configuration/source hashes.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'simulation','basissets','MRSCloud-main')));
mkdir(output_dir);
code=fileread(which('sim_signals'));
needle="outA = op_dccorr_FID_A(outA,'p');";
assert(numel(strfind(code,needle))==3);
code=strrep(code,'= sim_signals(MRS_opt)','= full_library_sim_signals(MRS_opt)');
code=strrep(code,needle,"save('-mat7-binary',fullfile(save_dir,'before_dc.mat'),'outA'); outA = op_dccorr_FID_A(outA,'p');");
fid=fopen(fullfile(output_dir,'full_library_sim_signals.m'),'w');fprintf(fid,'%s',code);fclose(fid);
addpath(output_dir);
fid=fopen(entries_file);entries=textscan(fid,'%s %s %f');fclose(fid);
for k=1:length(entries{1})
    label=entries{1}{k};metabolite=entries{2}{k};carrier=entries{3}(k);
    folder=fullfile(output_dir,label);cache=fullfile(folder,'response.mat');
    if exist(cache,'file')
        saved=load(cache,'signature');
        assert(strcmp(saved.signature,signature),'Cached raw response configuration/source mismatch');
        fprintf('CACHE %s\n',label);continue;
    end
    mkdir(folder);
    p=struct('flipAngle',180,'centreFreq',carrier,'edit_flipAngle',180,'excflipAngle',90,...
        'nX',spatial_points,'nY',spatial_points,'nZ',spatial_points,...
        'localization',{{'PRESS'}},'FieldStr',{{'3T'}},'vendor',{{'Universal_Philips'}},...
        'seq',{{'UnEdited'}},'metab',metabolite,'Nmetab',1,'TEs',{{30}},'save_dir',folder);
    p=load_parameters(p);p.Npts=points;
    spin_system=p.sys;
    [p,after]=full_library_sim_signals(p);
    before=load(fullfile(folder,'before_dc.mat'));
    raw_fid=before.outA.fids;raw_ppm=before.outA.ppm;after_first_dc=after.fids;
    sequence=struct('TE_ms',30,'TE1_ms',p.TE1,'TE2_ms',p.TE2,'RF_carrier_ppm',carrier,...
        'refocusing_duration_ms',p.refTp,'refocusing_waveform',p.refRF,'sw',p.sw,...
        'points',points,'spatial_points',spatial_points,'readout_linewidth_hz',p.lw,'B0_T',p.Bfield);
    save('-mat7-binary',cache,'raw_fid','raw_ppm','after_first_dc','sequence',...
        'spin_system','metabolite','carrier','signature');
    fprintf('COMPLETE %s\n',label);fflush(stdout);
end
end
