% Define base paths
meshmonkRoot = '..'; 

% Include directories
eigenIncludeDir = fullfile(meshmonkRoot, 'build', '_deps', 'eigen-src'); % Assumes Eigen is available via CMake's FetchContent. User might need to change this if their Eigen is elsewhere.
openMeshIncludeDir = fullfile(meshmonkRoot, 'build', 'vendor', 'OpenMesh-11.0.0', '_build', 'src'); % This path is based on the CMake build structure.
meshmonkIncludeDir = fullfile(meshmonkRoot, 'library', 'include');
meshmonkSrcDir = fullfile(meshmonkRoot, 'library', 'src');
meshmonkVendorDir = fullfile(meshmonkRoot, 'vendor'); % For nanoflann.hpp

% Library directory
libDir = fullfile(meshmonkRoot, 'build', 'library'); % Expects meshmonk_shared.lib here. The corresponding .dll should be in system's PATH or MATLAB's path at runtime.

% Common MEX flags
% For C++ standard, using COMPFLAGS. MSVC uses /std:c++17
% Using CXXFLAGS here as a placeholder if COMPFLAGS doesn't work as expected directly in the list.
% The typical way is: mex COMPFLAGS="\$COMPFLAGS /std:c++17" ...
% For this script, we'll add it to the list of flags directly.
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
    '-lmeshmonk_shared', ... % mex handles .lib extension automatically
    '-D_USE_MATH_DEFINES', ... % Common define for Windows
    'COMPFLAGS="/std:c++17"' % For MSVC compiler, ensure C++17 standard
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

disp('Starting MEX compilation for MeshMonk (Windows)...');

for i = 1:length(mexFiles)
    baseName = mexFiles{i};
    cppFile = fullfile(mexSrcDir, [baseName, '.cpp']);
    outputName = baseName; % Output name will be 'baseName.mexw64' or similar

    disp(['Mexing "', baseName, '" from "', cppFile, '"...']);
    
    % Construct the mex command
    % The flags in commonFlags are passed as separate arguments to mex.
    % COMPFLAGS is a special argument to mex.
    
    % Separate COMPFLAGS from other flags for clarity in command construction
    compileFlags = {};
    linkFlags = {};
    otherMexArgs = {};
    compFlagsValue = '';

    for k = 1:length(commonFlags)
        flag = commonFlags{k};
        if startsWith(flag, 'COMPFLAGS=')
            compFlagsValue = flag; % Keep it as is, e.g., 'COMPFLAGS="/std:c++17"'
        elseif startsWith(flag, '-L') || startsWith(flag, '-l')
            linkFlags{end+1} = flag; %#ok<AGROW>
        else
            compileFlags{end+1} = flag; %#ok<AGROW>
        end
    end
    
    mexCommandParts = {'mex'};
    if ~isempty(compFlagsValue)
        mexCommandParts{end+1} = compFlagsValue;
    end
    mexCommandParts = [mexCommandParts, compileFlags, linkFlags, {cppFile}];
    
    % Create the command string for display and evaluation
    mexCommand = strjoin(mexCommandParts, ' ');
    
    disp(['Executing: ', mexCommand]);
    
    try
        % Use feval for safer execution of mex with dynamic arguments
        argsToMex = [compileFlags, linkFlags, {cppFile}];
        if ~isempty(compFlagsValue)
             % Find the actual value for COMPFLAGS, e.g. "/std:c++17"
            compVal = extractBetween(compFlagsValue, '"', '"');
            if isempty(compVal) % Fallback if quotes were not used as expected
                compVal = extractAfter(compFlagsValue, '=');
            end
            feval('mex', 'COMPFLAGS', compVal{1}, argsToMex{:});
        else
            feval('mex', argsToMex{:});
        end
        disp(['Successfully mexed "', outputName, '".']);
    catch e
        disp(['Error mexing "', outputName, '":']);
        disp(e.message);
        disp(e.getReport('extended', 'on', 'hyperlinks', 'off'));
        disp('Skipping this file and continuing...');
    end
    disp('---');
end

disp('MEX compilation finished.');
