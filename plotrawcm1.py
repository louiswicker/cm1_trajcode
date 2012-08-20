#!/usr/bin/env python32
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
from optparse import OptionParser
import ctables
import trajsubs
from matplotlib.figure import Figure

x_limits = [32.,96]
y_limits = [64.,128.]
klevel   = 10
#---------------------------------------------------------------------------------------------------------------------------------
#
#  
def read_raw_cm1_state(time, um=0.0, vm=0.0):

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
  zc = 1000. * z[0:85]
  ze = 1000. * z[86:]
  
  print "XE-min:  %f km  XE-max:  %f km"   % (xe[0], xe[nx])
  print "YE-min:  %f km  YE-max:  %f km"   % (ye[0], ye[ny])
  print "ZE-min:  %f km  ZE-max:  %f km\n" % (ze[0], ze[nz])

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
  
  clevels = N.arange(-600,640,40)

  ax.contourf(x/1000., y/1000., wz, clevels, cmap=ctables.Not_PosDef_Default)
  ax.contour(x/1000., y/1000., wz, colors='k', levels=clevels[::2])
  
  ax.set_title('%s     XY Plot      Z = %.1f km   Max/Min/CINT = (%.1f / %.1f / %.1f)' % ('Vert Vort x 10^4/s', z, wz.max(), wz.min(), 2.0 ),fontsize=10 )
  ax.set_xlabel('East-West Distance (km)',fontsize=10)
  ax.set_ylabel('North-South Distance (km)',fontsize=10)
  ax.grid(True)
  
  if xlim != None:  ax.set_xlim(xlim[0],xlim[1])
  if ylim != None:  ax.set_ylim(ylim[0],ylim[1])

  return ax
#---------------------------------------------------------------------------------------------------------------------------------
#
#
def plot_w(w, x, y, ax=None, z = 1.0, xlim=None, ylim=None):

  if ax is None:
    ax = P.gca()

  clevels = N.arange(-30,32,2)

  ax.contourf(x/1000., y/1000., w, clevels, cmap=ctables.Not_PosDef_Default)
  ax.contour(x/1000., y/1000., w, colors='k', levels=clevels[::2])
  
  ax.set_title('%s     XY Plot      Z = %.1f km   Max/Min/CINT = (%.1f / %.1f / %.1f)' % ('W m/s', z, w.max(), w.min(), 2.0 ),fontsize=10 )
  ax.set_xlabel('East-West Distance (km)',fontsize=10)
  ax.set_ylabel('North-South Distance (km)',fontsize=10)
  ax.grid(True)
  
  if xlim != None:  ax.set_xlim(xlim[0],xlim[1])
  if ylim != None:  ax.set_ylim(ylim[0],ylim[1])

  return ax
  
#---------------------------------------------------------------------------------------------------------------------------------
#
#
def plot_u(u, x, y, ax=None, z = 1.0, xlim=None, ylim=None):

  if ax is None:
    ax = P.gca()

  clevels = N.arange(-30,32,2)

  ax.contourf(x/1000., y/1000., u, clevels, cmap=ctables.Not_PosDef_Default)
  ax.contour(x/1000., y/1000., u, colors='k', levels=clevels[::2])
  
  ax.set_title('%s     XY Plot      Z = %.1f km   Max/Min/CINT = (%.1f / %.1f / %.1f)' % ('U m/s', z, u.max(), u.min(), 2.0 ),fontsize=10 )
  ax.set_xlabel('East-West Distance (km)',fontsize=10)
  ax.set_ylabel('North-South Distance (km)',fontsize=10)
  ax.grid(True)

  if xlim != None:  ax.set_xlim(xlim[0],xlim[1])
  if ylim != None:  ax.set_ylim(ylim[0],ylim[1])


  return ax

#---------------------------------------------------------------------------------------------------------------------------------
#
#
def plot_v(v, x, y, ax=None, z = 1.0, xlim=None, ylim=None):

  if ax is None:
    ax = P.gca()

  clevels = N.arange(-30,32,2)

  ax.contourf(x/1000., y/1000., v, clevels, cmap=ctables.Not_PosDef_Default)
  ax.contour(x/1000., y/1000., v, colors='k', levels=clevels[::2])
  
  ax.set_title('%s     XY Plot      Z = %.1f km   Max/Min/CINT = (%.1f / %.1f / %.1f)' % ('V m/s', z, v.max(), v.min(), 2.0 ),fontsize=10 )
  ax.set_xlabel('East-West Distance (km)',fontsize=10)
  ax.set_ylabel('North-South Distance (km)',fontsize=10)
  ax.grid(True)
  if xlim != None:  ax.set_xlim(xlim[0],xlim[1])
  if ylim != None:  ax.set_ylim(ylim[0],ylim[1])


  return ax

#---------------------------------------------------------------------------------------------------------------------------------
#
#
def plot_dbz(dbz, x, y, ax=None, z = 1.0, xlim=None, ylim=None):

  if ax is None:
    ax = P.gca()

  clevels = N.arange(0,75,5)

  ax.contourf(x/1000., y/1000., dbz, clevels, cmap=ctables.NWSRef)
  ax.contour(x/1000., y/1000., dbz, colors='k', levels=clevels[::2])
  
  ax.set_title('%s     XY Plot      Z = %.1f km   Max/Min/CINT = (%.1f / %.1f / %.1f)' % ('Ref', z, dbz.max(), dbz.min(), 5.0 ),fontsize=10 )
  ax.set_xlabel('East-West Distance (km)',fontsize=10)
  ax.set_ylabel('North-South Distance (km)',fontsize=10)
  ax.grid(True)

  if xlim != None:  ax.set_xlim(xlim[0],xlim[1])
  if ylim != None:  ax.set_ylim(ylim[0],ylim[1])

  return ax

#####################################################################################################
if __name__ == "__main__":

  parser = OptionParser()
  parser.add_option("-t", "--time",       dest="time",      type="float",  help = "Start time of the trajectories")
  
  (options, args) = parser.parse_args()

  if options.time == None:
    print "\n ====>No trajectory time specified, exiting....\n"
    parser.print_help()
    sys.exit(-1)
  else:
    time = options.time

  state   = read_raw_cm1_state(time)
  
  z = state['zc'][klevel]/1000.
  
# Init some local variables

  fig = P.figure(figsize=(14,14))

  ax = fig.add_subplot(221)
  plot = plot_wz(state['u'][klevel,:,:], state['v'][klevel,:,:], state['xe'][:], state['ye'], z = z, ax=ax, xlim=x_limits, ylim=y_limits)
  ax = fig.add_subplot(222)
  plot = plot_w(state['w'][klevel+1,:,:], state['xc'][:], state['yc'], z = z, ax=ax, xlim=x_limits, ylim=y_limits)
  ax = fig.add_subplot(223)
  plot = plot_u(state['u'][klevel,:,:], state['xe'][:], state['yc'], z = z, ax=ax, xlim=x_limits, ylim=y_limits)
  ax = fig.add_subplot(224)
  plot = plot_v(state['v'][klevel,:,:], state['xc'][:], state['ye'], z = z, ax=ax, xlim=x_limits, ylim=y_limits)

  P.show()
  
  
