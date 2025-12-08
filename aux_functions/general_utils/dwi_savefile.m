
function dwi_savefile(file_magn,I_denoised)

fprintf('Saving file ... \n');

save_file_prefix = 'BM4D_PC_';

file_dir = dir(file_magn);
file_folder = file_dir.folder;


fname_tosave = [save_file_prefix file_dir.name];
fname_tosave = erase(fname_tosave,'.nii.gz');
fname_tosave = erase(fname_tosave,'.nii');


nii_info = niftiinfo([file_folder filesep file_dir.name]);

if contains(file_dir.name,'nii.gz')
    niftiwrite(cast(I_denoised,nii_info.Datatype),[file_folder filesep fname_tosave], nii_info, 'Compressed', true);
else
    niftiwrite(cast(I_denoised,nii_info.Datatype),[file_folder filesep fname_tosave], nii_info, 'Compressed', false);
end
