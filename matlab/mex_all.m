% Define base paths
meshmonkRoot = '..'; 

% Include directories
eigenIncludeDir = fullfile(meshmonkRoot, 'build', '_deps', 'eigen-src');
openMeshIncludeDir = fullfile(meshmonkRoot, 'build', 'vendor', 'OpenMesh-11.0.0', '_build', 'src');
meshmonkIncludeDir = fullfile(meshmonkRoot, 'library', 'include');
meshmonkSrcDir = fullfile(meshmonkRoot, 'library', 'src');
meshmonkVendorDir = fullfile(meshmonkRoot, 'vendor'); % For nanoflann.hpp

% Library directory
libDir = fullfile(meshmonkRoot, 'build', 'library');

% Common MEX flags
% Using C++17 standard, consistent with the main CMake build
% Adding -std=c++17. Add CXXFLAGS='-std=c++17' if MEX complains.
% For OpenMP, if used by MeshMonk: CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"
% For now, assuming no explicit OpenMP flags needed directly in mex command here if not in main build.
commonFlags = { ...
    ['-I', meshmonkIncludeDir], ...
    ['-I', fullfile(meshmonkIncludeDir, 'meshmonk')], ... % to access headers like meshmonk/global.hpp
    ['-I', meshmonkSrcDir], ... % if internal headers from library/src are directly included by MEX files
    ['-I', meshmonkVendorDir], ... % for nanoflann.hpp
    ['-I', eigenIncludeDir], ...
    ['-I', openMeshIncludeDir], ...
    ['-I', fullfile(openMeshIncludeDir, 'OpenMesh', 'Core')], ... % OpenMesh headers are often in subdirectories
    ['-I', fullfile(openMeshIncludeDir, 'OpenMesh', 'Tools')], ...
    ['-L', libDir], ...
    '-lmeshmonk_shared', ...
    '-lstdc++' ... % Common on Linux, might need adjustment for other OS
};

% Source files directory
mexSrcDir = fullfile(meshmonkRoot, 'matlab', 'mex');

% List of MEX files to compile (base names)
mexFiles = {
    'compute_correspondences', ...
    'compute_inlier_weights', ...
    'compute_nonrigid_transformation', ...
    'compute_rigid_transformation', ...
    'downsample_mesh', ...
    'nonrigid_registration', ...
    'pyramid_registration', ...
    'rigid_registration', ...
    'scaleshift_mesh', ...
    'compute_normals' ...
};

disp('Starting MEX compilation for MeshMonk...');

for i = 1:length(mexFiles)
    baseName = mexFiles{i};
    cppFile = fullfile(mexSrcDir, [baseName, '.cpp']);
    outputName = baseName; % Default output name, can be explicit with -output

    disp(['Mexing "', baseName, '" from "', cppFile, '"...']);
    
    % Construct the mex command
    % Note: Using '-output' explicitly ensures the output name matches the baseName,
    % which is good practice and avoids potential issues with platform-specific prefixes/extensions
    % if the .cpp file was named, e.g., pyramid_registration_mex.cpp but we want pyramid_registration as output.
    % The original script implies output name is derived from cpp filename.
    % The .cpp files are e.g. 'pyramid_registration.cpp', so output will be 'pyramid_registration.mexa64'
    
    mexCommand = ['mex ', sprintf('%s ', commonFlags{:}), ...
                  % '-output ', outputName, ' ', ... % Let MATLAB derive output name for now
                  cppFile];
    
    disp(['Executing: ', mexCommand]);
    
    try
        eval(mexCommand);
        disp(['Successfully mexed "', outputName, '".']);
    catch e
        disp(['Error mexing "', outputName, '":']);
        disp(e.message);
        disp('Skipping this file and continuing...');
    end
    disp('---');
end

disp('MEX compilation finished.');
