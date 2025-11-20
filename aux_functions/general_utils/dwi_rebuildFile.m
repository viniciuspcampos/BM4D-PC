
function dwi_rebuildFile(file_magn,temp_files_folder,I_denoised)

fprintf('Saving file and deleting temps... \n');

save_file_prefix = 'BM4D_PC_';

file_dir = dir(file_magn);
file_folder = file_dir.folder;


fname_tosave = [save_file_prefix file_dir.name];
fname_tosave = erase(fname_tosave,'.nii.gz');
fname_tosave = erase(fname_tosave,'.nii');


nii_info = niftiinfo([file_folder filesep file_dir.name]);

niftiwrite(cast(I_denoised,nii_info.Datatype),[file_folder filesep fname_tosave], nii_info, 'Compressed', false);



files_to_delete = dir([file_folder filesep temp_files_folder filesep '*.nii']);

for i=1:size(files_to_delete,1)
   delete([files_to_delete(i).folder filesep files_to_delete(i).name])
end

rmdir([file_folder filesep temp_files_folder])


