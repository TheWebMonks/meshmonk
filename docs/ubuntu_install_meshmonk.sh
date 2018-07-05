#!/bin/sh
#https://github.com/TheWebMonks/meshmonk/blob/master/docs/ubuntu.md

############# CHANGE CONFIG SETTINGS IF NEEDED #####################

# Set download directory (relative to users home directory)
DOWNLOADS_DIR=Downloads

# Change download URL for all necessary software packages
EIGEN_DOWNLOAD_FILE=3.3.4.zip
EIGEN_URL=http://bitbucket.org/eigen/eigen/get/$EIGEN_DOWNLOAD_FILE

OPENMESH_DOWNLOAD_FILE=OpenMesh-6.3.zip
OPENMESH_URL=https://www.openmesh.org/media/Releases/6.3/OpenMesh-6.3.zip

# Needed compiler version
CC=gcc-4.9
CXX=g++-4.9

# Specify log directory, relative to home directory
LOGDIR=logs
LOGFILE=meshmonk_installation.log
LOGFILE_EIGEN=eigen.log
LOGFILE_CODE_BLOCKS=code_blocks.log
LOGFILE_COMPILER=compiler.log
LOGFILE_OPENMESH=openmesh.log

##################### DO NOT EDIT BELOW #############################

# log installation process
log() {
    echo $1
    if [ -f $log_file ]; then
        echo `date` ' | ' $1 >> $log_file
    fi
}

info() {
    /bin/echo ""
    /bin/echo "  Usage: "
    /bin/echo " "
    /bin/echo "       Installing/downloading packages/libraries:"
    /bin/echo "           sudo $0 prepare"
    /bin/echo "       "
    /bin/echo "       Compiling MeshMonk:"
    /bin/echo "           sudo $0 compile"
    /bin/echo "       "
    /bin/echo "       Installing MeshMonk:"
    /bin/echo "           sudo $0 install"
    /bin/echo "       "

}

# Store action given as program argument
action=$1
if [ -z "$action" ]; then
    info
    exit 1
fi



# Check if scripts has been run as sudo...
# ref: https://askubuntu.com/a/30157/8698
checkRunAsRoot() {
    if ! [ $(id -u) = 0 ]; then
        echo "The script need to be run as root." >&2
        exit 1
    fi
}

# Set ownership of file to real user
setOwner() {
    file=$1
    log "$file: owner($real_user:$real_user_group)"
    chown $real_user:$real_user_group $file
}

# Set ownership of directory/file to real user
# Do this recursively for directories
# Preserve execution permission of user for group/other if requested
setOwnerAndFilePermissions() {
    dir=$1
    mod=$2
    recursive="-R"
    if [ -f $dir ]; then
        recursive=""
    fi

    log "$dir: owner ($real_user:$real_user_group)"
    chown $recursive $real_user:$real_user_group $dir

    log "$dir: permissions ($mod)"
    chmod $recursive $mod $dir

    if [ "$3" = "preserve_exec" ]; then
        log "$dir: preserve exec permissions"
        find $dir -perm /u+x -exec chmod a+x {} +
    fi

}

# Create directory if it does not exist
# Set ownership to real user and permissions as requested
createRealUserDirectory() {
    dir=$1
    mod=$2
    absdir=$real_user_home/$dir
    if [ ! -d $absdir ]; then
        log "$absdir: Creating directory"
        mkdir $absdir
        setOwnerAndFilePermissions $absdir $mod
    fi
}

# Initialize some stuff
#  - check if root runs this script
#  - get real user and group
#  - create log directory and file
#  - create download directory
initialize() {
    checkRunAsRoot

    date=`date +%Y-%m-%d`
    hostname=`hostname`

    # Get username who started initialization
    if [ $SUDO_USER ]; then
        real_user=$SUDO_USER
    else
        real_user=$(whoami)
    fi

    # Get primary group of the real user
    real_user_group=`id -Gn $real_user`
    set -- $real_user_group
    real_user_group=$1
    
    # Store home directory of the user who started installation
    real_user_home=`su $real_user -c 'echo $HOME'`

    # Prepare logging 
    log_dir=$real_user_home/$LOGDIR
    log_file=$log_dir/$LOGFILE

    createRealUserDirectory $LOGDIR 755

    # Rename last log file if exists
    if [ -f $log_file ]; then
        mv $log_file $log_file.$date
    fi

    # Create new log file
    touch $log_file
    chown $real_user:$real_user_group $log_file

    log "$real_user started MeshMonk installation on $hostname"

    # Create download directory
    createRealUserDirectory $DOWNLOADS_DIR 755
}

# Check if the needed compiler is installed
# install if necessary
checkCompilerVersions() {
    if [ -f /usr/bin/$CC ] && [ -f /usr/bin/$CXX ]; then
        log "$CC and $CXX already installed"
    else
        log "$CC or $CXX not found, installing now..."
        # install gcc/g++ 4.9
        DEBIAN_FRONTEND=noninteractive apt-get -yq install $CC $CXX > $log_dir/$LOGFILE_COMPILER
        setOwner $log_dir/$LOGFILE_COMPILER
    fi
}

# install Install code::blocks
installCodeBlocks() {
    log "installing code::blocks"
    add-apt-repository ppa:damien-moore/codeblocks-stable -y > $log_dir/$LOGFILE_CODE_BLOCKS 2>&1  
    DEBIAN_FRONTEND=noninteractive apt-get -yq update >> $log_dir/$LOGFILE_CODE_BLOCKS  
    DEBIAN_FRONTEND=noninteractive apt-get -yq install codeblocks codeblocks-contrib >> $log_dir/$LOGFILE_CODE_BLOCKS 
    setOwner $log_dir/$LOGFILE_CODE_BLOCKS
}

# Clone MeshMonk
#
# remove previous version, can be optimized by git pull...
# run a command as the normal user, use 'EOF' instead of EOF to avoid variable expansion
# 
# If any character in word is quoted, the delimiter shall be formed by performing quote 
# removal on word, and the here-document lines shall not be expanded. Otherwise, the 
# delimiter shall be the word itself.
#
# So the <<EOF version has the shell expand all variables before running the here doc 
# contents and the <<\EOF (or <<'EOF' or <<EO'F' etc.) versions don't expand the contents
#  (which lets bash in this case do that work).
cloneMeshMonk() {
    log "Clone MeshMonk"
    su $real_user <<EOF
        mkdir -p $real_user_home/projects
        cd $real_user_home/projects
        rm -rf meshmonk
        git clone https://github.com/TheWebMonks/meshmonk.git
EOF
}

# Install required libraries: Eigen, nanoflann and OpenMesh

# Download Eigen from EIGEN_URL to DOWNLOAD_DIR
installEigenLibrary() {
    log "Download Eigen"
    su $real_user <<EOF
        cd $real_user_home/$DOWNLOADS_DIR
        rm -rf EIGEN $EIGEN_DOWNLOAD_FILE
        wget $EIGEN_URL > $log_dir/$LOGFILE_EIGEN 2>&1
        unzip -o $EIGEN_DOWNLOAD_FILE -d EIGEN >> $log_dir/$LOGFILE_EIGEN 2>&1
EOF

    log "Copy include files, remove old files first"
    cd $real_user_home/$DOWNLOADS_DIR/EIGEN/eigen*
    eigendir=`pwd`
    rm -rf /usr/local/include/Eigen
    cp -r $eigendir/Eigen /usr/local/include/
    setOwnerAndFilePermissions /usr/local/include/Eigen "a+r" "preserve_exec"
}

# Nanoflann only needs header and 
#    is provided from meshmonk github in vendor directory
installNanoflann() {
    log "Copy Nanoflann include file, remove old file first"
    nanoflann_file=/usr/local/include/nanoflann.hpp
    rm -f $nanoflann_file
    cp $real_user_home/projects/meshmonk/vendor/nanoflann.hpp /usr/local/include/
    setOwnerAndFilePermissions $nanoflann_file "a+r" 
}


# Download OpenMesh from
installOpenMesh() {
    log "Download OpenMesh"
    su $real_user <<EOF
        cd $real_user_home/$DOWNLOADS_DIR
        rm -rf OpenMesh $OPENMESH_DOWNLOAD_FILE
        wget $OPENMESH_URL > $log_dir/$LOGFILE_OPENMESH 2>&1
        unzip -o $OPENMESH_DOWNLOAD_FILE -d OpenMesh >> $log_dir/$LOGFILE_OPENMESH 2>&1
        cd $real_user_home/$DOWNLOADS_DIR/OpenMesh/OpenMesh*
        openmesh_dir=`pwd`

        mkdir -p build
        cd build

        CC=$CC CXX=$CXX cmake .. >> $log_dir/$LOGFILE_OPENMESH 2>&1
        CC=$CC CXX=$CXX make >> $log_dir/$LOGFILE_OPENMESH 2>&1
EOF

    cd $real_user_home/$DOWNLOADS_DIR/OpenMesh/OpenMesh*
    openmesh_dir=`pwd`

    rm -rf /usr/local/include/OpenMesh
    cp -r $openmesh_dir/src/OpenMesh /usr/local/include/
    setOwnerAndFilePermissions /usr/local/include/OpenMesh "a+r" "preserve_exec"

    rm -f /usr/local/lib/libOpenMesh*
    for f in $openmesh_dir/build/Build/lib/*; do
        cp $f /usr/local/lib/
        setOwnerAndFilePermissions /usr/local/lib/$(basename $f) "a+r" "preserve_exec"
    done

    # Update library loader
    log "Update library loader - ldconfig"
    ldconfig
}

installMeshMonk() {
    cp $real_user_home/projects/meshmonk/bin/Release/libmeshmonk.so /usr/local/lib/
    setOwnerAndFilePermissions /usr/local/lib/libmeshmonk.so "a+r" "preserve_exec"
    (cd $real_user_home/projects/meshmonk/ && find . -name '*.hpp' -print | tar --create --files-from -) | (cd /usr/local/include/ && tar xvfp -)
    setOwnerAndFilePermissions /usr/local/include/src "a+r" "preserve_exec"
    setOwnerAndFilePermissions /usr/local/include/vendor "a+r" "preserve_exec"
    setOwnerAndFilePermissions /usr/local/include/meshmonk.hpp "a+r" 
    setOwnerAndFilePermissions /usr/local/include/global.hpp "a+r" 
        
    # Update library loader
    log "Update library loader - ldconfig"
    ldconfig
    
}

case "$action" in 
    prepare)
        initialize

        log "Preparing MeshMonk installation"
        checkCompilerVersions
        installCodeBlocks
        cloneMeshMonk
        installEigenLibrary
        installNanoflann
        installOpenMesh
        log " >>> Preparation finished, run program with 'compile' option now <<<"
        log ""
        ;;

    compile)
        echo "Compiling MeshMonk"
        echo ""
        echo "  run codeblocks command and follow the steps from"
        echo ""
        echo "  https://github.com/TheWebMonks/meshmonk/blob/master/docs/ubuntu.md"
        echo ""
        echo " >>> When compilation is finished, run program with 'install' option <<<"
        echo ""
        exit 0
        ;;

    install)
        initialize

        log "Installing MeshMonk"
        installMeshMonk

        echo ""
        echo " >>> MeshMonk installation complete <<<"
        echo ""
        ;;

    *)
        info
        exit 1
esac

log "Finished installation"
echo "Check your log file $log_file"
exit 0
