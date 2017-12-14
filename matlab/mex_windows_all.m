%Defining compiler flags
c_om = ['-I', '%ProgramFiles%\OpenMesh 6.3\include'];
c_mesh = ['-I', '%USERPROFILE%\Documents\GitHub\meshmonk'];
c_nano = ['-I', '%USERPROFILE%\Documents\GitHub\meshmonk\vendor'];
c_eigen = ['-I', '%USERPROFILE%\Documents\GitHub\eigen-eigen-3.3.4'];
c_math = ['-D', '_USE_MATH_DEFINES'];
c_lib = ['-l', 'meshmonk'];

disp('Mexing "compute_correspondences"...')
mex(c_om, c_mesh, c_nano, c_eigen, c_math, c_lib, 'matlab/mex/compute_correspondences.cpp')

disp('Mexing "compute_inlier_weights"...')
mex(c_om, c_mesh, c_nano, c_eigen, c_math, c_lib, 'matlab/mex/compute_inlier_weights.cpp')

disp('Mexing "compute_nonrigid_transformation"...')
mex(c_om, c_mesh, c_nano, c_eigen, c_math, c_lib, 'matlab/mex/compute_nonrigid_transformation.cpp')

disp('Mexing "compute_rigid_transformation"...')
mex(c_om, c_mesh, c_nano, c_eigen, c_math, c_lib, 'matlab/mex/compute_rigid_transformation.cpp')

disp('Mexing "downsample_mesh"...')
mex(c_om, c_mesh, c_nano, c_eigen, c_math, c_lib, 'matlab/mex/downsample_mesh.cpp')

disp('Mexing "nonrigid_registration"...')
mex(c_om, c_mesh, c_nano, c_eigen, c_math, c_lib, 'matlab/mex/nonrigid_registration.cpp')

disp('Mexing "pyramid_registration"...')
mex(c_om, c_mesh, c_nano, c_eigen, c_math, c_lib, 'matlab/mex/pyramid_registration.cpp')

disp('Mexing "rigid_registration"...')
mex(c_om, c_mesh, c_nano, c_eigen, c_math, c_lib, 'matlab/mex/rigid_registration.cpp')

disp('Mexing "scaleshift_mesh"...')
mex(c_om, c_mesh, c_nano, c_eigen, c_math, c_lib, 'matlab/mex/scaleshift_mesh.cpp')

disp('Mexing "compute_normals"...')
mex(c_om, c_mesh, c_nano, c_eigen, c_math, c_lib, 'matlab/mex/compute_normals.cpp')
