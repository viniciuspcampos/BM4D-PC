
function dwi_prepareFiles(file_magn,file_phase,temp_files_folder)



fprintf('Preparing files ... \n')

file_dir = dir(file_magn);


file_folder = file_dir.folder;

if isempty(file_phase)
    fprintf('Magnitude-only data ...  \n')
    I_noisy_to_denoise = (niftiread(file_magn));

else
    fprintf('Complex-valued data ... performing phase stabilization \n')
    I_noisy_magn = double(niftiread(file_magn));
    I_noisy_phase = double(niftiread(file_phase));

    I_noisy_to_denoise = perform_phase_stabilization(I_noisy_magn,I_noisy_phase);
   
end

nii_file_info = niftiinfo(file_magn);
nii_file_info.Datatype = 'single';


fname_tosave = 'dwi_temp';
file_temp_path = fullfile(file_folder , temp_files_folder , fname_tosave);
niftiwrite(cast(real(I_noisy_to_denoise),nii_file_info.Datatype),file_temp_path, nii_file_info, 'Compressed', false);




fprintf('Finished preparing files. \n\n');


