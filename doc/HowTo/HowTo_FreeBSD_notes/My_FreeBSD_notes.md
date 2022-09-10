## Notebook of what I had done to get FreeBSD up and running
* I started with the virtual machine image downloaded from [here](https://download.freebsd.org/ftp/releases/VM-IMAGES/12.2-RELEASE/amd64/Latest/FreeBSD-12.2-RELEASE-amd64.vmdk.xz)
* booted the 1st machine on this external virtual disk
* login: `root`
* adduser: `adduser wdklotz`
* adding a new member to a group using [__pw__](https://docs.freebsd.org/en/books/handbook/basics/#users-synopsis)
to allow **su** for user **wdklotz**: `pw groupmod wheel -M wdklotz`

# How to add a second hard disk on FreeBSD
* on VirtualBox add a new virtual disk 
![Virtual Media Manager](VBx1.gif)
* reboot; login as root
* use **sade** to format and create file system on virtual disk. [Instructions](https://www.cyberciti.biz/faq/freebsd-adding-second-hard-disk-howto/)...
* virtual disk will be **/dev/ada1p1** mounted on **/data**
* copy /usr to /data: `copy -R /usr/* /data`
* unmount /data: `umount /data`
* mount virtual disk again on mountpoint /usr: `mount /dev/ada1p1 /usr`

# Getting/Updating the Kernel Sources
* [Instructions](https://docs.freebsd.org/en_US.ISO8859-1/books/handbook/makeworld.html)...
* check the release: `uname -r`
* get the kernel sources: `svnlite checkout https://svn.freebsd.org/base/releng/12.2 /usr/src`
* after (very long!) download disk usage is: `du -sh /usr/src` **2.7G!** 

# Using the Ports Collection
To compile and install the port, change to the directory of the port to be installed, 
then type make install at the prompt. **lsof** needs kernel sources for compilation. Therefore I had to 
install the kernel sources as described above.
```
# cd /usr/ports/sysutils/lsof
# make install
# make clean
```

# Give Up!!!!
Ended with _disk full_ because I started from **VMDK** virtual disk 
which [cannot be resized](https://docs.oracle.com/en/virtualization/virtualbox/6.0/user/vboxmanage-modifymedium.html).

-----------------------------------------------
# Restart whole installation with VHD image from FreeBSD
* resized the virtual disk image from 4G to 100G with the **Virtual Media Manager**
![Virtual disk VHD resized from 4G to 100G](VBx2.gif)
* helpful reading about VirtualBox's disk images [here](https://docs.oracle.com/en/virtualization/virtualbox/6.0/user/vdidetails.html) and [here](https://docs.oracle.com/en/virtualization/virtualbox/6.0/user/vboxmanage-modifymedium.html)
* running in severe **problems** with `python -m pip install numpy scipy matplotlib`
* these python modules need BLAS/LAPACK libraries, which are FORTRAN libraries.
* FORTRAN comes with gcc. The make 
```
cd /usr/ports/math/lapack95
make install
```
therefore compiles and installs gcc10 and its dependencies (a very __long__ process).
* finally also this did not work. `python -m pip --no-cache-dir install numpy scipy` ended with unsatisdied library dependencies.

# Numpy, Scipy installation
* started from **python3.7** that I installed with ports (if I remember well) and made a 
virtual environment:
```
cd ~wdklotz
python -m venv py37v
alias py37='source /home/wdklotz/py37v/bin/activate'
py37
which python
```
* installed BLAS and LAPACK (which installs **gcc9** as dependency).
```
pkg install blas
pkg install lapack95
```
* linked **gcc** and **gfortran** to gcc9:
```
cd /usr/local/bin
ln -s gcc9 gcc
ln -s gfortran9 gfortran
```
* installed (with success this time) **numpy** and **scipy**.
```
python -m pip --no-cache-dir install numpy scipy
```

# Matplotlib
* But many **ERRORS** when trying `python -m pip --no-cache-dir matplotlib` after manupulations above.
* __impossible to install python modules this way!__
-----------------------------------------------
# The new Approach that worked!
* deactivate virtual python environment and remove it.
* as root install the missing modules with **pkg**:
```
cd
pkg install py37-matplotlib
pkg install py37-scipy
pkg install py37-pyaml
```
* as wdklotz check imports (all python imports are working!) and run SIMULINAC:
```
cd SIMULINAC
python3.7 simu.py
```
* create a virtual environment with access to the site-installation:
```
python3.7 -m venv --system-site-packages --symlinks py37v
```
* activate the virtual environment and run SIMULINAC:
```
export DISPLAY=192.168.1.52:0
source ~/py37v/bin/activate
cd SIMULINAC
python simu.py
```
* __bingo! -- SIMULINAC on a FreeBSD virtual machine__ without graphical desktop.

# How to run *BSD OS on Docker?
Docker doesn't actually run a full OS. Because it uses the host's kernel to run the container contents, it's not able to run a different kernel than the one used by its host OS. Further, as far as I understand, Docker relies on Linux-specific features for its fundamental operation. So it's not possible to run it with a BSD or another non-Linux kernel, including the XNU kernel used by MacOS, as its host environment. On a Mac, Docker actually runs within a virtualized Linux environment, so its host environment is Linux.

Now, in theory, if someone got a BSD userland to run on a Linux kernel, it might be possible to have a nearly-BSD Docker container. However, some research suggests that no project doing this has succeeded.

__All of that means that there's no way to run a true BSD as a Docker image, which is why there is no BSD image for Docker.__

# Adding [VirtualBoxGuestAdditions](https://docs.freebsd.org/en_US.ISO8859-1/books/handbook/virtualization-guest-virtualbox.html)
* had to install missing kernel-souces before running ports make (see above).