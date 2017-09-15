disp('Mexing "compute_correspondences"...')
mex -I/usr/local/include/ mex/compute_correspondences.cpp -lmeshmonk
disp('Mexing "compute_inlier_weights"...')
mex -I/usr/local/include/ mex/compute_inlier_weights.cpp -lmeshmonk
disp('Mexing "compute_nonrigid_transformation"...')
mex -I/usr/local/include/ mex/compute_nonrigid_transformation.cpp -lmeshmonk
disp('Mexing "compute_rigid_transformation"...')
mex -I/usr/local/include/ mex/compute_rigid_transformation.cpp -lmeshmonk
disp('Mexing "downsample_mesh"...')
mex -I/usr/local/include/ mex/downsample_mesh.cpp -lmeshmonk
disp('Mexing "nonrigid_registration"...')
mex -I/usr/local/include/ mex/nonrigid_registration.cpp -lmeshmonk
disp('Mexing "pyramid_registration"...')
mex -I/usr/local/include/ mex/pyramid_registration.cpp -lmeshmonk
disp('Mexing "rigid_registration"...')
mex -I/usr/local/include/ mex/rigid_registration.cpp -lmeshmonk
disp('Mexing "scaleshift_mesh"...')
mex -I/usr/local/include/ mex/scaleshift_mesh.cpp -lmeshmonk
disp('Mexing "compute_normals"...')
mex -I/usr/local/include/ mex/compute_normals.cpp -lmeshmonk