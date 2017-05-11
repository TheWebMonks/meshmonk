disp('Mexing "compute_correspondences"...')
mex -I/usr/local/include/ compute_correspondences.cpp -lmeshmonk
disp('Mexing "compute_inlier_weights"...')
mex -I/usr/local/include/ compute_inlier_weights.cpp -lmeshmonk
disp('Mexing "compute_nonrigid_transformation"...')
mex -I/usr/local/include/ compute_nonrigid_transformation.cpp -lmeshmonk
disp('Mexing "compute_rigid_transformation"...')
mex -I/usr/local/include/ compute_rigid_transformation.cpp -lmeshmonk
disp('Mexing "downsample_mesh"...')
mex -I/usr/local/include/ downsample_mesh.cpp -lmeshmonk
disp('Mexing "nonrigid_registration"...')
mex -I/usr/local/include/ nonrigid_registration.cpp -lmeshmonk
disp('Mexing "pyramid_registration"...')
mex -I/usr/local/include/ pyramid_registration.cpp -lmeshmonk
disp('Mexing "rigid_registration"...')
mex -I/usr/local/include/ rigid_registration.cpp -lmeshmonk
disp('Mexing "scaleshift_mesh"...')
mex -I/usr/local/include/ scaleshift_mesh.cpp -lmeshmonk
disp('Mexing "compute_normals"...')
mex -I/usr/local/include/ compute_normals.cpp -lmeshmonk