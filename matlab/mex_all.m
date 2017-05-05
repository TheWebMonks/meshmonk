disp('Mexing "compute_correspondences"...')
mex compute_correspondences.cpp -lmeshmonk
disp('Mexing "compute_inlier_weights"...')
mex compute_inlier_weights.cpp -lmeshmonk
disp('Mexing "compute_nonrigid_transformation"...')
mex compute_nonrigid_transformation.cpp -lmeshmonk
disp('Mexing "compute_rigid_transformation"...')
mex compute_rigid_transformation.cpp -lmeshmonk
disp('Mexing "downsample_mesh"...')
mex downsample_mesh.cpp -lmeshmonk
disp('Mexing "nonrigid_registration"...')
mex nonrigid_registration.cpp -lmeshmonk
disp('Mexing "pyramid_registration"...')
mex pyramid_registration.cpp -lmeshmonk
disp('Mexing "rigid_registration"...')
mex rigid_registration.cpp -lmeshmonk
disp('Mexing "scaleshift_mesh"...')
mex scaleshift_mesh.cpp -lmeshmonk
disp('Mexing "compute_normals"...')
mex compute_normals.cpp -lmeshmonk