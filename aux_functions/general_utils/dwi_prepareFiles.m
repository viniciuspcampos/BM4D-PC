
function I_noisy_phase_stab = dwi_prepareFiles(file_magn,file_phase)



fprintf('Preparing files ... \n')

if isempty(file_phase)
    fprintf('Magnitude-only data ...  \n')
    I_noisy_phase_stab = (niftiread(file_magn));

else
    fprintf('Complex-valued data ... performing phase stabilization \n')
    I_noisy_magn = double(niftiread(file_magn));
    I_noisy_phase = double(niftiread(file_phase));

    I_noisy_phase_stab = perform_phase_stabilization(I_noisy_magn,I_noisy_phase);
    I_noisy_phase_stab = real(I_noisy_phase_stab);
   
end


fprintf('Finished preparing files. \n\n');


