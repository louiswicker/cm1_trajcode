#!/usr/bin/env python
#
#  Traj.py is a trajectory program meant to replace my original IDL code
#     written in the 1990's for storm dynamics
#
#  Created by Lou Wicker March 2012
#

import numpy as N
import netCDF4 as ncdf
import sys
from time import time as timer
from optparse import OptionParser
import matplotlib.pylab as P
import ctables

# Fortran code imports
from trajsubs import trilinear, traj_integrate, read_cm1_raw


missing = -999.
debug = True

#---------------------------------------------------------------------------------------------------------------------------------
# Default parameters for the initialization of trajectories

default_xcntr  = 60450.
default_ycntr  = 48850.
default_zcntr  = 300.
default_radius = 1000.

#---------------------------------------------------------------------------------------------------------------------------------
#
#
class myTimer():

  def __init__(self, start = 0.0, name = 'nothing', minutes = False, args=[], kwargs={}):
    self.name    = name
    self.last    = timer()
    self.total   = 0.0
    self.minutes = minutes

  def start(self):
    self.last = timer()
      
  def stop(self):
    self.total = self.total + timer() - self.last

  def printit(self,string=None):
    print "\n  ---------------------------------------\n"
    if string == None: string = ""
    total = self.total + timer() - self.last
    if self.minutes:  
      total = round(total/60., 3) 
      print "%s ... Time in minutes for %s: %f" % (string, self.name, total)
    else:
      total = round(total, 3) 
      print "%s ... Time in seconds for %s: %f " % (string, self.name, total)

    print "\n  ---------------------------------------\n"

#---------------------------------------------------------------------------------------------------------------------------------
#
#
def trajinit(tstart, ntraj, trajfile = "test_traj.nc", xcntr = -1.0, ycntr = -1.0, zcntr= -1.0, radius=5000.):

  print "\nTRAJINIT:  CREATING A Horizontal ring of trajectories at %d \n" % tstart
  
  trajdata = N.zeros((ntraj,8))  # 0=time, 1=x, 2=y, 3=z, 4=u, 5=v, 6=w, 7=b, 8=dbz
  
  x0  = xcntr
  y0  = ycntr
  z0  = zcntr
  
  for n in N.arange(ntraj):
    trajdata[n,1] = x0 + radius*N.cos(2.0*N.pi*n/(ntraj))
    trajdata[n,2] = y0 + radius*N.sin(2.0*N.pi*n/(ntraj))
    trajdata[n,3] = z0

# Initialize a trajectory file in netCDF4

  rootgrp = ncdf.Dataset(trajfile, 'w', format='NETCDF4')
  time    = rootgrp.createDimension('time',  None)
  ntraj   = rootgrp.createDimension('ntraj', ntraj)
  
  time    = rootgrp.createVariable('time','f8',('time',))
  time.units = 'seconds since 001-01-01 00:00:00'
  x     = rootgrp.createVariable('x','f4',('time','ntraj',))            
  y     = rootgrp.createVariable('y','f4',('time','ntraj',))            
  z     = rootgrp.createVariable('z','f4',('time','ntraj',))            
  u     = rootgrp.createVariable('u','f4',('time','ntraj',))            
  v     = rootgrp.createVariable('v','f4',('time','ntraj',))            
  w     = rootgrp.createVariable('w','f4',('time','ntraj',))            
  b     = rootgrp.createVariable('b','f4',('time','ntraj',))            
  dbz   = rootgrp.createVariable('dbz','f4',('time','ntraj',))
  circ  = rootgrp.createVariable('circ','f4',('time',))
  buoy  = rootgrp.createVariable('buoy','f4',('time',))
  
# Write info into netcdf file....

  time[0]  = tstart
  x[0,:]   = trajdata[:,1] 
  y[0,:]   = trajdata[:,2] 
  z[0,:]   = trajdata[:,3] 
  u[0,:]   = -999.
  v[0,:]   = -999. 
  w[0,:]   = -999.
  b[0,:]   = -999.
  dbz[0,:] = -999.
  circ[0]  = -999.
  buoy[0]  = -999.
  
  rootgrp.close()
#---------------------------------------------------------------------------------------------------------------------------------
#
#
def trajinit_JD(tstart, trajfile = "test_traj.nc", input_traj_file="w96_ics_4980s.txt"):

  print "\nTRAJINIT_JD:  reading a trajectory input file from J. Dahl....", input_traj_file
  
  ic_data = N.loadtxt(input_traj_file, skiprows=2)
  
  ntraj_in = ic_data.shape[0]
  
  print "INPUT NTRAJ:  ", ntraj_in
  
  trajdata = N.zeros((ntraj_in,8))  # 0=time, 1=x, 2=y, 3=z, 4=u, 5=v, 6=w, 7=b, 8=dbz
  
# Initialize a trajectory file in netCDF4

  rootgrp = ncdf.Dataset(trajfile, 'w', format='NETCDF4')
  time    = rootgrp.createDimension('time',  None)
  ntraj   = rootgrp.createDimension('ntraj', ntraj_in)
  
  time    = rootgrp.createVariable('time','f8',('time',))
  time.units = 'seconds since 001-01-01 00:00:00'
  x     = rootgrp.createVariable('x','f4',('time','ntraj',))            
  y     = rootgrp.createVariable('y','f4',('time','ntraj',))            
  z     = rootgrp.createVariable('z','f4',('time','ntraj',))            
  u     = rootgrp.createVariable('u','f4',('time','ntraj',))            
  v     = rootgrp.createVariable('v','f4',('time','ntraj',))            
  w     = rootgrp.createVariable('w','f4',('time','ntraj',))            
#   b     = rootgrp.createVariable('b','f4',('time','ntraj',))            
#   dbz   = rootgrp.createVariable('dbz','f4',('time','ntraj',))
#   circ  = rootgrp.createVariable('circ','f4',('time',))
#   buoy  = rootgrp.createVariable('buoy','f4',('time',))
  
# Write info into netcdf file....

  for n in N.arange(ntraj_in):
    trajdata[n,1] = ic_data[n,2]
    trajdata[n,2] = ic_data[n,3]
    trajdata[n,3] = ic_data[n,4]
    time[0]       = ic_data[n,1]
    x[0,n]        = ic_data[n,2]
    y[0,n]        = ic_data[n,3] 
    z[0,n]        = ic_data[n,4]
    
  u[0,:]   = -999.
  v[0,:]   = -999. 
  w[0,:]   = -999.
#   b[0,:]   = -999.
#   dbz[0,:] = -999.
#   circ[0]  = -999.
#   buoy[0]  = -999.
  
  rootgrp.close()
  
#---------------------------------------------------------------------------------------------------------------------------------
#
#  
def read_raw_cm1_state(time, nstep, dt, um=0.0, vm=0.0):

# Read in needed data from CM1 Files....
  
  t0 = time/60
  t1 = (time+nstep*dt)/60
  nx = 512
  ny = 512
  nz = 86
  dx = 250.
  dy = 250.

  xc = dx*(0.5 + N.arange(nx))
  yc = dy*(0.5 + N.arange(ny))
  xe = dx*(N.arange(nx+1))
  ye = dy*(N.arange(ny+1))
  
  z  = N.loadtxt("z_input.data")
  zc = z[0:nz]
  ze = z[nz:]
  
  print "XE-min:  %f km  XE-max:  %f km"   % (xe[0], xe[nx])
  print "YE-min:  %f km  YE-max:  %f km"   % (ye[0], ye[ny])
  print "ZE-min:  %f km  ZE-max:  %f km\n" % (ze[0], ze[nz])

# Read first file

  file_time = 542 + (time - 4980)/2
  
  u0, v0, w0 = read_cm1_raw(file_time, nx, ny, nz)
  u0 = u0 - um
  v0 = v0 - vm
  
  print "\nRead RAW CM1 time: %d (min) File: %.6d  W-min: %f m/s  w-max:  %f m/s" % (t0,file_time, w0.max(), w0.min())

# Read second file

  file_time = 542 + (time+nstep*dt - 4980)/2
  
  u1, v1, w1 = read_cm1_raw(file_time, nx, ny, nz)
  u1 = u1 - um
  v1 = v1 - vm
  
  print "Read RAW CM1 time: %d (min) File: %.6d  W-min: %f m/s  w-max:  %f m/s\n" % (t1,file_time, w1.max(), w1.min())
  
  return {'t0': t0*60, 't1': t1*60, 'nx': nx, 'ny': ny, 'nz': nz, 
          'xc': xc, 'yc': yc, 'zc': zc*1000., 'xe': xe, 'ye': ye, 'ze': ze*1000.,
          'u0': u0, 'v0': v0, 'w0': w0, 'u1': u1, 'v1': v1, 'w1': w1,
          'b0': None, 'b1': None}
#---------------------------------------------------------------------------------------------------------------------------------
#
#  
def read_cm1_state(prefix, time, nstep, dt, um=0.0, vm=0.0):

# Read in needed data from CM1 Files....

  t0 = time/60
  t1 = (time + nstep*dt)/60
  
  file = "%s_%.6d.nc" % (prefix, t0)
  f0   = ncdf.Dataset(file, "r")
  file = "%s_%.6d.nc" % (prefix, t1)
  f1   = ncdf.Dataset(file, "r")
  
  xc = f1.variables['xh'][:]
  yc = f1.variables['yh'][:]
  zc = f1.variables['z'][:]
  xe = f1.variables['xf'][:]
  ye = f1.variables['yf'][:]
  ze = f1.variables['zf'][:]
  
  nx = len(f1.dimensions['ni'])
  ny = len(f1.dimensions['nj'])
  nz = len(f1.dimensions['nk'])
  
  print "XE-min:  %f km  XE-max:  %f km"   % (xe[0], xe[nx])
  print "YE-min:  %f km  YE-max:  %f km"   % (ye[0], ye[ny])
  print "ZE-min:  %f km  ZE-max:  %f km\n" % (ze[0], ze[nz])

  u0 = f0.variables['u'][:] - um
  v0 = f0.variables['v'][:] - vm
  w0 = f0.variables['w'][:]
  b0 = f0.variables['b'][:]
  
  print "\nRead CM1 time: %d (min) File: %s_%.6d.nc  W-min: %f m/s  w-max:  %f m/s" % (t0, prefix, t0, w0.max(), w0.min())

  u1 = f1.variables['u'][:] - um
  v1 = f1.variables['v'][:] - vm
  w1 = f1.variables['w'][:]
  b1 = f1.variables['b'][:]
  
  print "Read CM1 time: %d (min) File: %s_%.6d.nc  W-min: %f m/s  w-max:  %f m/s\n" % (t1, prefix, t1, w1.max(), w1.min())
  
  return {'t0': t0*60, 't1': t1*60, 'nx': nx, 'ny': ny, 'nz': nz, 
          'xc': xc, 'yc': yc, 'zc': zc, 'xe': xe, 'ye': ye, 'ze': ze,
          'u0': u0, 'v0': v0, 'w0': w0, 'u1': u1, 'v1': v1, 'w1': w1,
          'b0': b0, 'b1': b1}
          
          
#---------------------------------------------------------------------------------------------------------------------------------
#
#
def traj_interpolate(x, y, z, times, state, ivars=None):

  ntraj = x.shape[1]
  nstep = times.shape[0]
  nx = state['nx']
  ny = state['ny']
  nz = state['nz']
  
  xc = state['xc'][:]
  xe = state['xe'][:]
  yc = state['yc'][:]
  ye = state['ye'][:]
  zc = state['zc'][:]
  ze = state['ze'][:]
  
  if debug:
    print "\nTraj_Interp\n"
    print nx, ny, nz, ntraj, nstep
    print xe[0], xc[0], xc[-1], xe[-1]
    print ye[0], yc[0], yc[-1], ye[-1]
    print ze[0], zc[0], zc[-1], ze[-1]
    
# create arrays to stuff data into

  data = {'u':N.zeros((nstep,ntraj)), 'v':N.zeros((nstep, ntraj)), 'w':N.zeros((nstep,ntraj)), 'b':N.zeros((nstep,ntraj))}
  
  if ivars == None:
    ivars = ["u", "v", "w"]
  
# Compute the time weights

  t10 = (state['t1'] - times) / (state['t1'] - state['t0'])
  t11 = 1.0 - t10
  
  for m, key in enumerate( ivars ):
  
  # time interp of state variable
        
    if key == 'u':    
      tmp1 = trilinear(x.ravel(), y.ravel(), z.ravel(), state['u0'][:], xe, yc, zc, 0.0, 0.0, missing)
      tmp2 = trilinear(x.ravel(), y.ravel(), z.ravel(), state['u1'][:], xe, yc, zc, 0.0, 0.0, missing)

    if key == 'v':
      tmp1 = trilinear(x.ravel(), y.ravel(), z.ravel(), state['v0'][:], xc, ye, zc, 0.0, 0.0, missing)
      tmp2 = trilinear(x.ravel(), y.ravel(), z.ravel(), state['v1'][:], xc, ye, zc, 0.0, 0.0, missing)

    if key == 'w':
      tmp1 = trilinear(x.ravel(), y.ravel(), z.ravel(), state['w0'][:], xc, yc, ze, 0.0, 0.0, missing)
      tmp2 = trilinear(x.ravel(), y.ravel(), z.ravel(), state['w1'][:], xc, yc, ze, 0.0, 0.0, missing)
      
    if key == 'b':
      tmp1 = trilinear(x.ravel(), y.ravel(), z.ravel(), state['b0'][:], xc, yc, zc, 0.0, 0.0, missing)
      tmp2 = trilinear(x.ravel(), y.ravel(), z.ravel(), state['b1'][:], xc, yc, zc, 0.0, 0.0, missing)
        
    for n in N.arange(nstep): 
      data[key][n,:] = t10[n]*tmp1.reshape(nstep,ntraj)[n,:] + t11[n]*tmp2.reshape(nstep,ntraj)[n,:]
      
# Now compute a Simpson's Rule Integration around the circut for circulation and buoyancy

  xm = N.zeros((nstep,ntraj-1))
  ym = N.zeros((nstep,ntraj-1))
  zm = N.zeros((nstep,ntraj-1))

# Find mid-points of circ lines....

#   timer_ComputeCirc = myTimer(name = "Time to compute circulation")
# 
#   for n in N.arange(ntraj-1):
#     xm[:,n] = (x[:,n] + x[:,n+1])/2.0
#     ym[:,n] = (y[:,n] + y[:,n+1])/2.0
#     zm[:,n] = (z[:,n] + z[:,n+1])/2.0
#     
#   u1 = trilinear(xm.ravel(), ym.ravel(), zm.ravel(), state['u0'][:], xe, yc, zc, 0.0, 0.0, missing)
#   u2 = trilinear(xm.ravel(), ym.ravel(), zm.ravel(), state['u1'][:], xe, yc, zc, 0.0, 0.0, missing)
# 
#   v1 = trilinear(xm.ravel(), ym.ravel(), zm.ravel(), state['v0'][:], xc, ye, zc, 0.0, 0.0, missing)
#   v2 = trilinear(xm.ravel(), ym.ravel(), zm.ravel(), state['v1'][:], xc, ye, zc, 0.0, 0.0, missing)
# 
#   w1 = trilinear(xm.ravel(), ym.ravel(), zm.ravel(), state['w0'][:], xc, yc, ze, 0.0, 0.0, missing)
#   w2 = trilinear(xm.ravel(), ym.ravel(), zm.ravel(), state['w1'][:], xc, yc, ze, 0.0, 0.0, missing)
#     
# #   b1 = trilinear(xm.ravel(), ym.ravel(), zm.ravel(), state['b0'][:], xc, yc, zc, 0.0, 0.0, missing)
# #   b2 = trilinear(xm.ravel(), ym.ravel(), zm.ravel(), state['b1'][:], xc, yc, zc, 0.0, 0.0, missing)
# 
#   circ = N.zeros((nstep))
#   buoy = N.zeros((nstep))
#   
#   for n in N.arange(nstep):
#     um = t10[n]*u1.reshape(nstep,ntraj-1)[n,:] + t11[n]*u2.reshape(nstep,ntraj-1)[n,:]
#     vm = t10[n]*v1.reshape(nstep,ntraj-1)[n,:] + t11[n]*v2.reshape(nstep,ntraj-1)[n,:]
#     wm = t10[n]*w1.reshape(nstep,ntraj-1)[n,:] + t11[n]*w2.reshape(nstep,ntraj-1)[n,:]
#     bm = t10[n]*b1.reshape(nstep,ntraj-1)[n,:] + t11[n]*b2.reshape(nstep,ntraj-1)[n,:]
#     up = data['u'][n,:]
#     vp = data['v'][n,:]
#     wp = data['w'][n,:]
#     bp = data['b'][n,:]
#     xp = x[n,:]
#     yp = y[n,:]
#     zp = z[n,:]
# 
#     circ[n] = N.sum((up[0:ntraj-1] + 4.0*um[0:ntraj-1] + up[1:ntraj]) * (xp[1:ntraj] - xp[0:ntraj-1])/6.0 \
#                   + (vp[0:ntraj-1] + 4.0*vm[0:ntraj-1] + vp[1:ntraj]) * (yp[1:ntraj] - yp[0:ntraj-1])/6.0 \
#                   + (wp[0:ntraj-1] + 4.0*wm[0:ntraj-1] + wp[1:ntraj]) * (zp[1:ntraj] - zp[0:ntraj-1])/6.0 )
#       
#     buoy[n] = N.sum((bp[0:ntraj-1] + 4.0*bm[0:ntraj-1] + bp[1:ntraj]) * (zp[1:ntraj] - zp[0:ntraj-1])/6.0 )
#     
#   timer_ComputeCirc.printit()
    
  return data['u'], data['v'], data['w'], data['b'], 0.0, 0.0
 
#---------------------------------------------------------------------------------------------------------------------------------
#
#
def traj_compute(x, y, z, dt, state):

  ntraj = x.shape[0]
  nx = state['nx']
  ny = state['ny']
  nz = state['nz']
  
  xc = state['xc'][:]
  xe = state['xe'][:]
  yc = state['yc'][:]
  ye = state['ye'][:]
  zc = state['zc'][:]
  ze = state['ze'][:]
  
  u0 = state['u0'][:]
  v0 = state['v0'][:]
  w0 = state['w0'][:]
  u1 = state['u1'][:]
  v1 = state['v1'][:]
  w1 = state['w1'][:]
  
  t0 = state['t0']
  t1 = state['t1']
  
# Stuff to make sure the time step goes evenly into t0 & t1
  
  step = (t1 - t0) / dt
  
  chk = step.as_integer_ratio()
  
# We add 1 to nstep to cleanly step between t0 and t1
  
  if chk[0] == chk[1]:
    nstep = 1 + N.int(step)
    ndt   = dt
  else:
    ndt   = (t1 - t0) / step
    nstep = 1 + N.int((t1 - t0) / ndt)
    print "\nTraj_Integrate:  dt does not divide into interval, new dt = %f\n" % ndt

  if debug:
    print "\nTraj_Integrate\n"
    print nx, ny, nz, ntraj, nstep
    print t0, t1, ndt
    print xe[0], xc[0], xc[-1], xe[-1]
    print ye[0], yc[0], yc[-1], ye[-1]
    print ze[0], zc[0], zc[-1], ze[-1]

# Create arrays for output

  times = t0 + dt * N.arange(nstep)
  xnew  = N.zeros((nstep,ntraj))
  ynew  = N.zeros((nstep,ntraj))
  znew  = N.zeros((nstep,ntraj))

# For each trajectory, integrate nsteps
  
  if debug:  print "\n Now calling fortran subroutine to integrate trajectories \n"
  
# Call fortran program to compute RK4 time integration

  xnew, ynew, znew, unew, vnew, wnew = traj_integrate(x, y, z, u0, v0, w0, u1, v1, w1, xc, xe, yc, ye, zc, ze, ndt, t0, t1, 0.0, 0.0, nstep)
  
  return times, xnew, ynew, znew

#---------------------------------------------------------------------------------------------------------------------------------
#
#
def traj_start(time, trajfile):

# Open the trajectory file

  tfile = ncdf.Dataset(trajfile, "r")
  
# Read the x/y/z location for the input time....

  try:
    tindex = N.where(time == tfile.variables['time'][:])
    tindex = tindex[0][0]
  except:
    print "TRAJ --> Input time value is:  %d " % (time)
    print "TRAJ --> Time values in trajectory file are %f" % (tfile.variables['time'][:])
    print "TRAJ --> exiting traj"
    sys.exit(-1)
  
  x = tfile.variables['x'][tindex,:]
  y = tfile.variables['y'][tindex,:]
  z = tfile.variables['z'][tindex,:]
  
  tfile.close()

  return x, y, z  

#---------------------------------------------------------------------------------------------------------------------------------
#
#
def traj_end(times, xnew, ynew, znew, u, v, w, b, circ, buoy, trajfile):

# Open the trajectory file

  tfile = ncdf.Dataset(trajfile, "r+")
  
# Read the x/y/z location for the input time....

  try:
    tindex = N.where(times[0] == tfile.variables['time'][:])
    tindex = tindex[0][0]
  except:
    print "TRAJ --> Input time value is:  %d " % (times[0])
    print "TRAJ --> Time values in trajectory file are %f" % (tfile.variables['time'][:])
    print "TRAJ --> exiting traj"
    sys.exit(-1)
    
  tfile.variables['time'][tindex:] = times
  tfile.variables['x'][tindex:,:]  = xnew
  tfile.variables['y'][tindex:,:]  = ynew
  tfile.variables['z'][tindex:,:]  = znew
  tfile.variables['u'][tindex:,:]  = u
  tfile.variables['v'][tindex:,:]  = v
  tfile.variables['w'][tindex:,:]  = w
#   tfile.variables['b'][tindex:,:]  = b
#   tfile.variables['circ'][tindex:] = circ
#   tfile.variables['buoy'][tindex:] = buoy
  
  tfile.close()

  return 0 
#---------------------------------------------------------------------------------------------------------------------------------
#
#  
def traj_check(time, trajfile, state=None):

# Open the trajectory file

  tfile = ncdf.Dataset(trajfile, "r")
  
# Read the x/y/z location for the input time....

  try:
    tindex = N.where(time == tfile.variables['time'][:])
    tindex = tindex[0][0]
  except:
    print "TRAJ --> Input time value is:  %d " % (time)
    print "TRAJ --> Time values in trajectory file are %f" % (tfile.variables['time'][:])
    print "TRAJ --> exiting traj"
    sys.exit(-1)
  
  t = tfile.variables['time'][tindex:]
  x = tfile.variables['x'][tindex:,:]
  y = tfile.variables['y'][tindex:,:]
  z = tfile.variables['z'][tindex:,:]
  u = tfile.variables['u'][tindex:,:]
  v = tfile.variables['v'][tindex:,:]
  w = tfile.variables['w'][tindex:,:]
#   b = tfile.variables['b'][tindex:,:]
  
  ntraj = x.shape[1]
  nstep = x.shape[0]

  fig = P.figure(figsize=(10,10))

  ax = fig.add_subplot(111)
  
  if state !=None:
    clevels = N.arange(-100,105,10)
    ax.contourf(state['xe'][:]/1000., state['yc'][:]/1000., state['u0'][10,:,:], clevels, cmap=ctables.Not_PosDef_Default)
    ax.contour(state['xc'][:]/1000., state['ye'][:]/1000., state['v0'][10,:,:], colors='k', levels=clevels[::2])
  
  ax.set_title('Test Traj   XY Plot      Z = %.1f km' % (1000.),fontsize=10 )
  ax.set_xlabel('East-West Distance (km)',fontsize=10)
  ax.set_ylabel('North-South Distance (km)',fontsize=10)
  
  for n in N.arange(ntraj):
    xp = x[:nstep:10,n]/1000.
    yp = y[:nstep:10,n]/1000.
    ax.plot(xp, yp,color="k", linewidth=2.)
    
  ax.set_xlim(state['xe'][0]/1000.,state['xe'][-1]/1000.)
  ax.set_ylim(state['ye'][0]/1000.,state['ye'][-1]/1000.)

  P.show()
  
# Check for interpolation accuracy.....

  x0 = state['x0']
  y0 = state['y0']
  z0 = state['z0']
  omega = state['omega']

# Print out ending stats for each trajectory

  print "Stats for each trajectory"
  for n in N.arange(ntraj):
    b2 = N.sqrt( (x[-1,n]-x0)**2 + (y[-1,n]-y0)**2 + (z[-1,n]-z0)**2 )
    u2 = - omega*(y[-1,n]-y0)
    v2 =   omega*(x[-1,n]-x0)
    w2 =   (x[n,0] - x0) / 2500.
    print "Traj # %d  x:%.3f  y:%.3f  z:%.3f  P-B:  %.3f  A-B:  %.3f  P-U:  %.3f  A-U:  %.3f  P-V:  %.3f  A-V:  %.3f  P-W:  %.3f  A-W:  %.3f" % \
         (n, x[-1,n], y[-1,n], z[-1,n], b[-1,n],b2, u[-1,n],u2, v[-1,n],v2, w[-1,n], w2)


  print "Stats for 0th trajectory at T=0 and T=3600"

  for n in [0,nstep-1]:
    b2 = N.sqrt( (x[n,0]-x0)**2 + (y[n,0]-y0)**2 + (z[n,0]-z0)**2 )
    u2 = - omega*(y[n,0]-y0)
    v2 =   omega*(x[n,0]-x0)
    w2 =   (x[n,0] - x0) / 2500.
    print "Time %.2f  x:%.3f  y:%.3f  z:%.3f  P-B:  %.3f  A-B:  %.3f  P-U:  %.3f  A-U:  %.3f  P-V:  %.3f  A-V:  %.3f  P-W:  %.3f  A-W:  %.3f" % \
         (t[n], x[n,0]-x0, y[n,0]-y0, z[n,0], b[n,0],b2, u[n,0],u2, v[n,0],v2, w[n,0], w2)
  
  print "The X / Y /z start and finish positions differing by less than 0.1 meters \n"

#---------------------------------------------------------------------------------------------------------------------------------
#
#  
def create_fake_velo():

# Create analytical flow fields to test the trajectory 

  xdomain = 100000.
  ydomain = 100000.
  zdomain = 10000.
  
  delta    = 3600.
  omega    = 2.0*N.pi / delta
  
  nx = 100
  ny = 100
  nz = 50
  
  x0 = xdomain/2.0
  y0 = ydomain/2.0
  z0 = zdomain/2.0
  
  dx = xdomain / N.float(nx)
  dy = ydomain / N.float(ny)
  dz = zdomain / N.float(nz)
  
  xc = dx/2.0 + dx*N.arange(nx)
  yc = dy/2.0 + dy*N.arange(ny)
  zc = dz/2.0 + dz*N.arange(nz)
  xe = dx*N.arange(nx+1)
  ye = dy*N.arange(ny+1)
  ze = dz*N.arange(nz+1)
  
  print "XE-min:  %f km  XE-max:  %f km"   % (xe[0], xe[nx])
  print "YE-min:  %f km  YE-max:  %f km"   % (ye[0], ye[ny])
  print "ZE-min:  %f km  ZE-max:  %f km\n" % (ze[0], ze[nz])
  print "XC-min:  %f km  XC-max:  %f km"   % (xc[0], xc[nx-1])
  print "YC-min:  %f km  YC-max:  %f km"   % (yc[0], yc[ny-1])
  print "ZC-min:  %f km  ZC-max:  %f km\n" % (zc[0], zc[nz-1])
  
# Create a circular flow field in center of the domain

  u0 = N.ones((nz,ny,nx+1))
  v0 = N.ones((nz,ny+1,nx))
  w0 = N.ones((nz+1,ny,nx))
  b0 = N.ones((nz,ny,nx))

  for j in N.arange(ny):
    u0[:,j,:] = - omega * (yc[j] - y0)
    
  for i in N.arange(nx):
    v0[:,:,i] = omega * (xc[i] - x0)

  for i in N.arange(nx):  
    w0[:,:,i] = (xc[i] - x0) / 2500.
  
  for k in N.arange(nz):  
    for j in N.arange(ny):
      for i in N.arange(nx):
        b0[k,j,i] =  N.sqrt( (xc[i]-x0)**2 + (yc[j]-y0)**2 + (zc[k]-z0)**2 )

  print "\nCreate fake data  U-min: %f m/s  U-max:  %f m/s" % (u0.min(), u0.max())
  print "\nCreate fake data  V-min: %f m/s  V-max:  %f m/s" % (v0.min(), v0.max())
  print "\nCreate fake data  W-min: %f m/s  W-max:  %f m/s" % (w0.min(), w0.max())
  print "\nCreate fake data  B-min: %f m/s  B-max:  %f m/s" % (b0.min(), b0.max())

  u1 = u0.copy()
  v1 = v0.copy()
  w1 = w0.copy()
  b1 = b0.copy()
  t0 = 3600.
  t1 = t0 + delta
  dt = (t1 - t0) / delta
  
  print t0, t1, dt
  
  return {'t0': t0, 't1': t1, 'nx': nx, 'ny': ny, 'nz': nz, 
          'xc': xc, 'yc': yc, 'zc': zc, 'xe': xe, 'ye': ye, 'ze': ze,
          'u0': u0, 'v0': v0, 'w0': w0, 'u1': u1, 'v1': v1, 'w1': w1,
          'b0': b0, 'b1': b1,
          'omega': omega, 'x0':x0, 'y0':y0, 'z0':z0, 'dt':dt}

#####################################################################################################
if __name__ == "__main__":


  parser = OptionParser()
  parser.add_option("-f", "--file",       dest="file",      type="string",  help = "Name of netCDF trajectory file")
  parser.add_option("-t", "--time",       dest="time",      type="float",   help = "Start time of the trajectories")
  parser.add_option(      "--prefix",     dest="prefix",    type="string",  help = "Prefix of input files")
  parser.add_option("-p", "--path",       dest="path",      type="string",  help = "Path for input files")
  parser.add_option(      "--dt",         dest="dt",        type="float",   help = "Time step for the trajectory calcs")
  parser.add_option(      "--um",         dest="um",        type="float",   help = "Storm motion in X-direction (m/s)")
  parser.add_option(      "--vm",         dest="vm",        type="float",   help = "Storm motion in Y-direction (m/s)")
  parser.add_option(      "--nstep",      dest="nstep",     type="int",     help = "Number of time steps to take")
  parser.add_option(      "--ntraj",      dest="ntraj",     type="int",     help = "Number of initial trajectories, overwrites trajfile")
  parser.add_option(      "--test",       dest="test",      action="store_true", help = "Run testing of traj code")
  
  
  (options, args) = parser.parse_args()
  
  if options.test:
    print "\n ====>TESTING TRAJECTORY CODE WITH ANALYTICAL TRAJECTORIES\n"
    state = create_fake_velo()
    trajinit(3600, 1, trajfile = "fake_traj.nc", xcntr = 50000., ycntr = 50000., zcntr = 5000., radius = 5000.)
    x, y, z = traj_start(3600, "fake_traj.nc")
    times, xnew, ynew, znew = traj_compute(x, y, z, state['dt'], state)
    u, v, w, b = traj_interpolate(xnew, ynew, znew, times, state)
    ret = traj_end(times, xnew, ynew, znew, u, v, w, b, "fake_traj.nc")
    traj_check(3600, "fake_traj.nc", state=state)
    print "\n ====>COMPLETED TESTING TRAJECTORY CODE WITH ANALYTICAL TRAJECTORIES\n" 
    sys.exit(0)

  if options.prefix == None:
    print "\n ====>No prefix for run files specified, exiting....\n"
    parser.print_help()
    sys.exit(-1)
  else:
    prefix = options.prefix
    
  if options.file == None:
    print "\n ====>No trajectory file specified, exiting....\n"
    parser.print_help()
    sys.exit(-1)
  else:
    file = options.file
      
  if options.path == None:
    path = "./"
  else:
    path = options.path
    
  if options.time == None:
    print "\n ====>No initial trajectory time specified, exiting....\n"
    parser.print_help()
    sys.exit(-1)
  else:
    time = options.time
    
  if options.dt == None:
    print "\n ====>No time step set for trajectory calculations, exiting....\n"
    parser.print_help()
    sys.exit(-1)
  else:
    dt = options.dt

  if options.ntraj != None:
    print "\n ====>User has specified number of trajectories, creating a new file....\n"
    trajinit_JD(time, trajfile = file)

  if options.nstep == None:
    print "\n ====>Number of steps for trajectory calculations = 0, exiting....\n"
    parser.print_help()
    sys.exit(-1)
  else:
    nstep = options.nstep
    
  if options.um == None:
    print "\n ====>No U-box motion supplied, assuming zero\n"
    um = 0.0
  else:
    um = options.um

  if options.vm == None:
    print "\n ====>No V-box motion supplied, assuming zero\n"
    vm = 0.0
  else:
    vm = options.vm
  
# Initialize timers

  timer_Total           = myTimer(name = "Total CPU Time")
  timer_ReadFiles       = myTimer(name = "Time to read in 3D Files")
  timer_IntegrateTraj   = myTimer(name = "Time to Integrate Trajectories")
  timer_InterpolateTraj = myTimer(name = "Time to Interpolate fields to Trajectories")

  timer_Total.start()
  
# Using prefix, create list of files to do integration.....

  timer_ReadFiles.start()
  
  state   = read_raw_cm1_state(time, nstep, dt, um=um, vm=vm)
  
  timer_ReadFiles.printit()
  
# Open the trajectory file

  x, y, z = traj_start(time, file)

# Integrate trajs

  timer_IntegrateTraj.start()
  
  times, xnew, ynew, znew = traj_compute(x, y, z, dt, state)
  
  timer_IntegrateTraj.printit()
  
# Integrate trajs

  timer_InterpolateTraj.start()

  u, v, w, b, circ, buoy = traj_interpolate(xnew, ynew, znew, times, state)

  timer_InterpolateTraj.printit()
    
# Write the data back out....

  ret = traj_end(times, xnew, ynew, znew, u, v, w, b, circ, buoy, file)
  
  print "\n Trajectories written out to the trajectory file, exiting....."
  
  timer_Total.printit()


 

   

  
