# SIMULINAC as docker container
All work done in SIMULINAC's base directory
```
$cd $HOME/SIMULINAC
```
Make a Dockerfile:
```
FROM python:3

WORKDIR /Simulinac
ENV DISPLAY=192.168.1.52:0
VOLUME /Simulinac

RUN python -m pip install --no-cache-dir matplotlib numpy PyQt5 PyYAML scipy

CMD ["python", "simu.py"]
```
1. WORKDIR sets the containers base directory to `/Simulinac`.
2. ENV sets the `DISPLAY` environment variable for graphics output. Here to the host display.
3. VOLUME sets the `path to the persistent data`. Here the WORKDIR.
4. RUN `python -m pip` installs python dependencies.
5. CMD defines the command in `JSON syntax` the container runs when started. 

Build the docker image with `-t name:tag`:
```
$docker build -t simu:1.2 .
```

When build is done the new image should be visible:
```
$docker images
REPOSITORY               TAG       IMAGE ID       CREATED        SIZE
simu                     1.2       106c17d0c60d   5 hours ago    1.3GB
debian                   latest    5890f8ba95f6   4 weeks ago    114MB
portainer/portainer-ce   latest    980323c8eb3f   2 months ago   196MB
jekyll/jekyll            latest    76e17ded11d1   2 months ago   656MB
eclipse-mosquitto        latest    131df074d5c8   2 months ago   9.54MB
docker/getting-started   latest    021a1b85e641   3 months ago   27.6MB
```

Start the image (1st time) using the default command given during build:
```
$docker run -it -v ${PWD}:/Simulinac --name simu simu:1.2
```
1. `-it` attach a terminal for input/output (optional)
2. `-v ${PWD}:/Simulinac` bind the working directory to the volume defined during build.
3. `--name simu` give the container a name.
4. `simu:1.2` the image to run.

When the container has executed `CMD ["python", "simu.py"]` it will exit:
```
$ docker container ls -l
CONTAINER ID   IMAGE      COMMAND            CREATED         STATUS                     PORTS     NAMES
2ff2e236a4df   simu:1.2   "python simu.py"   4 minutes ago   Exited (0) 3 minutes ago             simu
$
```

You can run the image to create a second container with a shell and verify the `VOLUME` is mounted:
```
wdklotz@xps8500:/mnt/c/Users/wdklotz/SIMULINAC$ docker run -it -v ${PWD}:/Simulinac simu:1.2 sh
# pwd
/Simulinac
# ls
DYNAC_lattice_generator.py  TTFG.py           dynac                 lattice.py                 roadmap.txt
Dockerfile                  __init__.py       dynacEzTab            lattice_generator.py       setutil.py
DynacG.py                   __pycache__       dynacIN               marker_actions.py          sigma.py
Ez0.py                      bucket_size.py    dynacv6_0_wdk_modifs  nb-configuration.xml       simu.py
OXAL.py                     bunch.py          elements.py           notes.txt                  trackPlot.py
README.md                   changelog.txt.py  formulae.py           pyOR_lattice_generator.py  tracker.py
README_docker.md            doc               frames                py_testing                 xml_utils
SF_WDK2g44.TBL              docksimu          jupyter.ipynb         requirements.txt           yml
#
```

When you exit the shell with `Ctrl+d` the container will exit. You can remove it now.

You can now edit your project files on the host and then start the container again with your modifications:
```
$docker start -a simu
simu.py v8.0.6a2 on python 3.9.2 on linux
run Version 20.02.2019_nlat
   macros=macros_20.02.2019_nlat
   template=tmpl_20.02.2019_nlat
   input=yml/simuIN.yml
READING SF-DATA from "SF_WDK2g44.TBL"
12 gap-slices, 192 SF-intervals, 16 SF-intervals/gap-slice
DEBUG: Ez_poly: (il,i0,ir) (  0,  8, 16),       (zl,zr) (-3.84,-3.2)
DEBUG: Ez_poly: (il,i0,ir) ( 16, 24, 32),       (zl,zr) (-3.2,-2.56)
DEBUG: Ez_poly: (il,i0,ir) ( 32, 40, 48),       (zl,zr) (-2.56,-1.92)
DEBUG: Ez_poly: (il,i0,ir) ( 48, 56, 64),       (zl,zr) (-1.92,-1.28)
DEBUG: Ez_poly: (il,i0,ir) ( 64, 72, 80),       (zl,zr) (-1.28,-0.64)
DEBUG: Ez_poly: (il,i0,ir) ( 80, 88, 96),       (zl,zr) (-0.64,0.0)
DEBUG: Ez_poly: (il,i0,ir) ( 96,104,112),       (zl,zr) (0.0,0.64)
DEBUG: Ez_poly: (il,i0,ir) (112,120,128),       (zl,zr) (0.64,1.28)
DEBUG: Ez_poly: (il,i0,ir) (128,136,144),       (zl,zr) (1.28,1.92)
                          :
                       use aperture : False
               use emittance growth : False
                        use express : False
                   use ring lattice : False
                 use sigma tracking : False
CALCULATE C+S TRAJECTORIES
CALCULATE TWISS ENVELOPES
PREPARE DISPLAY
$
```

That's it. You have a docker container that runs the simu.py simulator with graphics, which you can develop/edit on the host.