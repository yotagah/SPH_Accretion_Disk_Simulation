# SPH High Speed Collision Simulation

This Smoothed Particle Hydrodynamics C program simulates an accretion disk around a white dwarf in a binary dwarf novae system. It is a variation of the code in https://github.com/yotagah/SPH_Accretion_Disk_Simulation. The code and all parameters are based in the references below:

Liu G. R. and Liu M. B., Smoothed particle hydrodynamics: a meshfree particle method, World Scientific, 3rd printing, 2003.

Simpson J. C., Numerical techniques for three-dimensional smoothed particle hydrodynamics simulations: applications to accretion disks, The Astrophysical Journal, 448: 822-831, 1995 August 1.

Libersky L. D. at al., High strain lagrangian hydrodynamics: a three-dimensional SPH code for dynamic material response, Journal of Computational Physics 109, 67-75, 1993.

Warner B., Cataclysmic Variable Stars, Cambridge University Press, Cambridge, UK, 1995

# Usage

Download and extract or clone the project, then run make in the terminal to compile.

```sh
$ make
```

After the compile finishes, run the sph executable directing the output to a file data.dat for instance.

```sh
$ ./sph > data.dat
```

The simulation will take a while and when it finishes, run the load 'plot.gnu' command into GNUPLOT to generate the "frames" of the simulation in the pict folder (which has to be created first) or directly in the screen.

```sh
gnuplot> load 'plot.gnu'
```

A video with the result of a similar simulation (with more particles and time than the code in this repository) can be found in Youtube: https://www.youtube.com/watch?v=T3VKhmvu-gQ