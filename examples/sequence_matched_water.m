function sequence_matched_water(output_dir, spatial_points, points)
% Isolated MRSCloud PRESS experiment. All sequence amplitudes remain unscaled.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'simulation','basissets','MRSCloud-main')));
mkdir(output_dir);
metabs = {'NAA','Ref0ppm','H2O','H2O','H2O','H2O'};
carriers = [3.0,3.0,4.7,3.0,4.7,4.7];
labels = {'NAA','Ref0ppm','water_on','water_off','water_single_proton','water_repeat'};
signals = struct();
for k=1:length(metabs)
    p = struct('flipAngle',180,'centreFreq',carriers(k),'edit_flipAngle',180,...
        'excflipAngle',90,'nX',spatial_points,'nY',spatial_points,'nZ',spatial_points,...
        'localization',{{'PRESS'}},'FieldStr',{{'3T'}},'vendor',{{'Universal_Philips'}},...
        'seq',{{'UnEdited'}},'metab',metabs{k},'Nmetab',1,'TEs',{{30}},...
        'save_dir',fullfile(output_dir,labels{k}));
    mkdir(p.save_dir);
    p = load_parameters(p);
    p.Npts = points;
    if k==5
        % Control only: halve the already two-proton initial density matrix.
        for j=1:length(p.d), p.d{j}=p.d{j}/2; end
    end
    [p,signal] = sim_signals(p);
    signals.(labels{k}) = signal.fids;
    if k==1
        TE1=p.TE1; TE2=p.TE2; refTp=p.refTp; refRF=p.refRF;
        sw=p.sw; Bfield=p.Bfield; lw=p.lw;
    end
end
save('-mat7-binary',fullfile(output_dir,'raw_sequence.mat'),'signals','carriers',...
     'TE1','TE2','refTp','refRF','sw','Bfield','lw','spatial_points','points');
end
