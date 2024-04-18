Optimal 2d kspace project
by David Frey and Tao Hong

This code requires the following packages to run:
  toppe v6 (git@github.com:toppeMRI/toppe.git)
  mirt (git@github.com:JeffFessler/mirt.git)
  minTimeGradient (https://people.eecs.berkeley.edu/~mlustig/Software.html)

sequences.spgr() will generate toppe v6 files for an SPGR sequence using a user-created 2d/3d kspace trajectory

test.m will run the reconstruction pipeline on phantom data for a user-created 2d kspace trajectory

kspace trajectories can be created outside of this package, and saved as kspace.mat in the current directory
the variable "kspace" in the .mat file should be of dimensions [Npoints x Nshots x Ndims]
