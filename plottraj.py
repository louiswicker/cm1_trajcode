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
import matplotlib.pylab as P
import matplotlib as mpl
from optparse import OptionParser
import ctables
from mpl_toolkits.mplot3d import Axes3D
import trajsubs

y_limits = [96.,108.]
x_limits = [62.,76.]
z_limits = [0.,3.]

cmin, cmax = (0,3.0)
cstep = 0.2
_nskip = 1  # typical number of trajs to skip

# definitions for the plot layout with multiple panels
left, width = 0.1, 0.5
bottom, height = 0.1, 0.5
bottom_h = left_h = left+width+0.03

rectC = [left, bottom, width, height]
rectX = [left, bottom_h, width, 0.2]
rectY = [left_h, bottom, 0.2, height]

#---------------------------------------------------------------------------------------------------------------------------------
#
#
def read_raw_cm1_state(prefix, time, um=0.0, vm=0.0):

# Read in needed data from CM1 Files....

  print time

  file_time = 542 + (time - 4980)/2

  t0 = time/60
  nx = 512
  ny = 512
  nz = 86
  dx = 250.
  dy = 250.
  dz = 200.

  xc = dx*(0.5 + N.arange(nx))
  yc = dy*(0.5 + N.arange(ny))
  xe = dx*(N.arange(nx+1))
  ye = dy*(N.arange(ny+1))
  z  = N.loadtxt("z_input.data")
  zc = 1000.*z[0:nz]
  ze = 1000.*z[nz:]

  print "XE-min:  %f km  XE-max:  %f km"   % (xe[0], xe[nx])
  print "XC-min:  %f km  XC-max:  %f km"   % (xc[0], xc[nx-1])
  print "YE-min:  %f km  YE-max:  %f km"   % (ye[0], ye[ny])
  print "YC-min:  %f km  YC-max:  %f km"   % (yc[0], yc[ny-1])
  print "ZE-min:  %f km  ZE-max:  %f km\n" % (ze[0], ze[nz])
  print "ZC-min:  %f km  ZC-max:  %f km\n" % (zc[0], zc[nz-1])

  u0, v0, w0 = trajsubs.read_cm1_raw(file_time, nx, ny, nz)

  u0 = u0 - um
  v0 = v0 - vm

  print "\nRead CM1 time: %d (min) File: %s_%.6d.nc  W-min: %f m/s  w-max:  %f m/s" % (t0, file_time, t0, w0.max(), w0.min())

  return {'t0': t0*60, 'nx': nx, 'ny': ny, 'nz': nz,
          'xc': xc, 'yc': yc, 'zc': zc, 'xe': xe, 'ye': ye, 'ze': ze,
          'u': u0, 'v': v0, 'w': w0, 'b': None, 'dbz': None }
#---------------------------------------------------------------------------------------------------------------------------------
#
#  
def read_cm1_state(prefix, time,  um=0.0, vm=0.0):

# Read in needed data from CM1 Files....

  t0 = 1 + time/60
  
  file = "%s_%.6d.nc" % (prefix, t0)
  f0   = ncdf.Dataset(file, "r")
  
  xc = f0.variables['xh'][:]
  yc = f0.variables['yh'][:]
  zc = f0.variables['z'][:]
  xe = f0.variables['xf'][:]
  ye = f0.variables['yf'][:]
  ze = f0.variables['zf'][:]
  
  nx = len(f0.dimensions['ni'])
  ny = len(f0.dimensions['nj'])
  nz = len(f0.dimensions['nk'])
  
  print "XE-min:  %f km  XE-max:  %f km"   % (xe[0], xe[nx])
  print "YE-min:  %f km  YE-max:  %f km"   % (ye[0], ye[ny])
  print "ZE-min:  %f km  ZE-max:  %f km\n" % (ze[0], ze[nz])

  u0 = f0.variables['u'][:] - um
  v0 = f0.variables['v'][:] - vm
  w0 = f0.variables['w'][:]
  b0 = f0.variables['b'][:]
  d0 =  f0.variables['dbz'][:]
  
  print "\nRead CM1 time: %d (min) File: %s_%.6d.nc  W-min: %f m/s  w-max:  %f m/s" % (t0, prefix, t0, w0.max(), w0.min())
  
  return {'t0': t0*60, 'nx': nx, 'ny': ny, 'nz': nz, 
          'xc': xc, 'yc': yc, 'zc': zc, 'xe': xe, 'ye': ye, 'ze': ze,
          'u': u0, 'v': v0, 'w': w0, 'b': b0, 'dbz': d0 }
#---------------------------------------------------------------------------------------------------------------------------------
#
#
def traj_read(time, trajfile):

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
# b = tfile.variables['b'][tindex:,:]
  
  tfile.close()

  return t, x, y, z, u, v, w, w 
  
#---------------------------------------------------------------------------------------------------------------------------------
#
#
def plot_w(w, x, y, ax=None, z = 1.0, xlim=None, ylim=None):

  if ax is None:
    ax = P.gca()

  clevels = N.arange(-20,22,2)

  ax.contourf(x/1000., y/1000., w, clevels, cmap=ctables.Not_PosDef_Default)
  ax.contour(x/1000., y/1000., w, colors='k', levels=clevels[::2])
  
  ax.set_title('%s     XY Plot      Z = %.1f km   Max/Min/CINT = (%.1f / %.1f / %.1f)' % ('W m/s', z, w.max(), w.min(), 2.0 ),fontsize=10 )
  ax.set_xlabel('East-West Distance (km)',fontsize=10)
  ax.set_ylabel('North-South Distance (km)',fontsize=10)
  
  if xlim != None:  ax.set_xlim(xlim[0],xlim[1])
  if ylim != None:  ax.set_ylim(ylim[0],ylim[1])

  return ax
#---------------------------------------------------------------------------------------------------------------------------------
#
#
def plot_b(b, x, y, ax=None, z = 1.0, xlim=None, ylim=None):

  if ax is None:
    ax = P.gca()

  clevels = N.arange(-0.2,0.22,0.02)

  ax.contourf(x/1000., y/1000., b, clevels, cmap=ctables.Not_PosDef_Default)
  ax.contour(x/1000., y/1000., b, colors='k', levels=clevels[::2])
  
  ax.set_title('%s     XY Plot      Z = %.1f km   Max/Min/CINT = (%.1f / %.1f / %.1f)' % ('B m^2/s', z, b.max(), b.min(), 2.0 ),fontsize=10 )
  ax.set_xlabel('East-West Distance (km)',fontsize=10)
  ax.set_ylabel('North-South Distance (km)',fontsize=10)
  
  if xlim != None:  ax.set_xlim(xlim[0],xlim[1])
  if ylim != None:  ax.set_ylim(ylim[0],ylim[1])

  return ax

#---------------------------------------------------------------------------------------------------------------------------------
#
#
def plot_wz(u, v, x, y, ax=None, z = 1.0, xlim=None, ylim=None):

  if ax is None:
    ax = P.gca()

# Compute vertical vorticity on the staggered grid....

  nx = v.shape[1]
  ny = u.shape[0]
  wz   = N.zeros((ny+1, nx+1))
  
# assume dx and dy are constant

  dx = x[2] - x[1]
  dy = y[2] - y[1]
  
  wz[1:ny-1, 1:nx-1] = (v[1:ny-1, 2:nx  ] - v[1:ny-1, 1:nx-1]) / dx \
                     - (u[2:ny,   1:nx-1] - u[1:ny-1, 1:nx-1]) / dy
    
# fill boundaries to make the plot look nice

  wz[0, :] = wz[1,   :]
  wz[ny,:] = wz[ny-1,:]
  wz[:, 0] = wz[:,   1]
  wz[:,nx] = wz[:,nx-1]
                    
  wz = wz * 10000.
  
  clevels = N.arange(-600,620,20)

  ax.contourf(x/1000., y/1000., wz, clevels, cmap=ctables.Not_PosDef_Default)
  ax.contour(x/1000., y/1000., wz, colors='k', levels=clevels[::2])
  
  ax.set_title('%s     XY Plot      Z = %.1f km   Max/Min/CINT = (%.1f / %.1f / %.1f)' % ('Vert Vort x 10^4/s', z, wz.max(), wz.min(), 2.0 ),fontsize=12, position = (0.6, -0.125))
  ax.set_xlabel('East-West Distance (km)',fontsize=10)
  ax.set_ylabel('North-South Distance (km)',fontsize=10)
  ax.grid(True)
  
  if xlim != None:  ax.set_xlim(xlim[0],xlim[1])
  if ylim != None:  ax.set_ylim(ylim[0],ylim[1])

  return ax
  
#####################################################################################################
if __name__ == "__main__":

  parser = OptionParser()
  parser.add_option("-f", "--file",       dest="file",      type="string",  help = "Name of netCDF trajectory file")
  parser.add_option("-t", "--time",       dest="time",      type="float",  help = "Start time of the trajectories")
  parser.add_option(     "--klevel",      dest="klevel",    type="int",  help = "K-index of horizontal plane to plot as background")
  parser.add_option(      "--prefix",     dest="prefix",    type="string",  help = "Prefix of input files")
  
  (options, args) = parser.parse_args()

  if options.file == None:
    print "\n ====>No trajectory file specified, exiting....\n"
    parser.print_help()
    sys.exit(-1)
  else:
    file = options.file

  if options.time == None:
    print "\n ====>No trajectory time specified, exiting....\n"
    parser.print_help()
    sys.exit(-1)
  else:
    time = options.time

  if options.prefix == None:
    print "\n ====>No prefix for cm1 file specified, exiting....\n"
    parser.print_help()
    sys.exit(-1)
  else:
    prefix = options.prefix

  if options.klevel == None:
    print "\n ====>No k-level cm1 file specified, setting klevel = 1\n"
    klevel = 3
  else:
    klevel = options.klevel
        
  t, x, y, z, u, v, w, b = traj_read(time, file)
  
  state = read_raw_cm1_state(prefix, time, um=10, vm=3.)
  
  ntraj = x.shape[1]
  nstep = x.shape[0]
  print "NTRAJ:  ",ntraj
  print "NSTEP:  ",nstep

  fig = P.figure(figsize=(12,12))

# Using contourf to provide my colorbar info, then clearing the figure
# Setting up a colormap that's a simple transtion
  mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',['red','blue'])
  nothing = [[0,0],[0,0]]
  levels = N.arange(cmin,cmax+cstep,cstep)
  CS3 = P.contourf(nothing, levels, cmap=mymap)
  P.clf()

  axC = P.axes(rectC)
  axX = P.axes(rectX)
  axY = P.axes(rectY)

  plot = plot_wz(state['u'][klevel,:,:], state['v'][klevel,:,:], state['xe'], state['ye'], z = state['zc'][klevel], ax=axC, xlim=x_limits, ylim=y_limits)

  for n in N.arange(0,ntraj,_nskip):
    xp = x[:nstep:10,n]/1000.
    yp = y[:nstep:10,n]/1000.
    b = (z[nstep-1,n]/1000.-cmin)/(cmax-cmin)
    g = 0.0
    r = 1-b
    axC.plot(xp, yp, color=(r,g,b), linewidth=1.)

  axC.set_xlim(x_limits[0], x_limits[1])
  axC.set_ylim(y_limits[0], y_limits[1])

  cb = P.colorbar(CS3)
  cb.set_label("Height of Parcel Location (km) at Final Location")

# XZ projection

  for n in N.arange(0,ntraj,_nskip):
    xp = x[:nstep:10,n]/1000.
    zp = z[:nstep:10,n]/1000.
    b = (z[nstep-1,n]/1000.-cmin)/(cmax-cmin)
    g = 0.0
    r = 1-b
    axX.plot(xp, zp, color=(r,g,b), linewidth=1.)

  axX.set_xlim(x_limits[0], x_limits[1])
  axX.set_ylim(z_limits[0], z_limits[1])
  axX.set_title("XZ Projection")

# YZ projection

  for n in N.arange(0,ntraj,_nskip):
    yp = y[:nstep:10,n]/1000.
    zp = z[:nstep:10,n]/1000.
    b = (z[nstep-1,n]/1000.-cmin)/(cmax-cmin)
    g = 0.0
    r = 1-b
    axY.plot(zp, yp, color=(r,g,b), linewidth=1.)

  axY.set_xlim(z_limits[0], z_limits[1])
  axY.set_ylim(y_limits[0], y_limits[1])
  axY.set_title("YZ Projection")

  P.show()

  fig.savefig("plottraj.pdf")
