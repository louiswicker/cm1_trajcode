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

x_limits = [40.,80.]
y_limits = [30.,70.]

#---------------------------------------------------------------------------------------------------------------------------------
#
#  
def read_cm1_state(prefix, time,  um=0.0, vm=0.0):

# Read in needed data from CM1 Files....

  t0 = time/60
  
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
def traj_read(time, trajfile, circ=False, buoy=False):

# Open the trajectory file

  print trajfile

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
  b = tfile.variables['b'][tindex:,:]
  
  if circ == True:
    circ = tfile.variables['circ'][tindex:]
  if buoy == True:
    buoy = tfile.variables['buoy'][tindex:]
  
  tfile.close()

  return t, x, y, z, u, v, w, b, circ, buoy
  
#---------------------------------------------------------------------------------------------------------------------------------
#
#
def plot_w(w, x, y, ax=None, z = 1.0, xlim=None, ylim=None):

  if ax is None:
    ax = P.gca()

  clevels = N.arange(-20,22,2)

  ax.contourf(x/1000., y/1000., w, clevels, cmap=ctables.Not_PosDef_Default)
  ax.contour(x/1000., y/1000., w, colors='k', levels=clevels[::2])
  
#   ax.set_title('%s     XY Plot      Z = %.1f km   Max/Min/CINT = (%.1f / %.1f / %.1f)' % ('W m/s', z, w.max(), w.min(), 2.0 ),fontsize=10 )
  ax.set_xlabel('East-West Distance (km)',fontsize=10)
  ax.set_ylabel('North-South Distance (km)',fontsize=10)
  
  if xlim != None:  ax.set_xlim(xlim[0],xlim[1])
  if ylim != None:  ax.set_ylim(ylim[0],ylim[1])

  return ax
 
#####################################################################################################
if __name__ == "__main__":

  parser = OptionParser()
  parser.add_option("-f", "--file",       dest="file",      type="string",  help = "Name of netCDF trajectory file")
  parser.add_option("-t", "--time",       dest="time",      type="float",  help = "Start time of the trajectories")
  parser.add_option(      "--prefix",     dest="prefix",    type="string",  help = "Prefix of input files")
  parser.add_option("-p", "--path",       dest="path",      type="string",  help = "Path for input files")
  parser.add_option(      "--dt",         dest="dt",        type="float",   help = "Time step for the trajectory calcs")
  parser.add_option(      "--nplot",      dest="nplot",     type="int",     help = "Number of time steps to plot")
  
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

  if options.nplot == None:
    nplot = 4
  else:
    nplot = options.nplot
        
  t, x, y, z, u, v, w, b, circ2, buoy2 = traj_read(time, file, circ=True, buoy=True)
  
  state = read_cm1_state(prefix, time, um=11.2, vm=2.2)
  
  ntraj = x.shape[1]
  nstep = x.shape[0]
  
  np = N.arange(nplot) * nstep / (nplot)

  fig = P.figure(figsize=(16,16))

  ax = fig.add_subplot(221)

  plot = plot_w(state['w'][10,:,:], state['xc'][:], state['yc'], ax=ax, xlim=x_limits, ylim=y_limits)  
  
  xp = x[np[0],0:ntraj]/1000.
  yp = y[np[0],0:ntraj]/1000.
  ax.plot(xp, yp,color="k", linewidth=2.)
  ax.set_title('Time = %d ' % (N.int(t[np[0]])), fontsize=10 )
    
  ax.set_xlim(x_limits[0],x_limits[1])
  ax.set_ylim(y_limits[0],y_limits[1])
  
  ax = fig.add_subplot(222)

  plot = plot_w(state['w'][10,:,:], state['xc'][:], state['yc'], ax=ax, xlim=x_limits, ylim=y_limits)  
  
  xp = x[np[1],0:ntraj]/1000.
  yp = y[np[1],0:ntraj]/1000.
  ax.plot(xp, yp,color="k", linewidth=2.)
  ax.set_title('Time = %d ' % (N.int(t[np[1]])),fontsize=10 )

    
  ax.set_xlim(x_limits[0],x_limits[1])
  ax.set_ylim(y_limits[0],y_limits[1])

  ax = fig.add_subplot(223)

  plot = plot_w(state['w'][10,:,:], state['xc'][:], state['yc'], ax=ax, xlim=x_limits, ylim=y_limits)  
  
  xp = x[np[2],0:ntraj]/1000.
  yp = y[np[2],0:ntraj]/1000.
  ax.plot(xp, yp,color="k", linewidth=2.)
  ax.set_title('Time = %d ' % (N.int(t[np[2]])),fontsize=10 )
    
  ax.set_xlim(x_limits[0],x_limits[1])
  ax.set_ylim(y_limits[0],y_limits[1])
 
  ax = fig.add_subplot(224)

  plot = plot_w(state['w'][10,:,:], state['xc'][:], state['yc'], ax=ax, xlim=x_limits, ylim=y_limits)  
  
  xp = x[np[3],0:ntraj]/1000.
  yp = y[np[3],0:ntraj]/1000.
  ax.plot(xp, yp,color="k", linewidth=2.)
  ax.set_title('Time = %d ' % (N.int(t[np[3]])),fontsize=10 )
    
  ax.set_xlim(x_limits[0],x_limits[1])
  ax.set_ylim(y_limits[0],y_limits[1])

  P.show()
  
  circ  = N.zeros((nstep))
  buoy  = N.zeros((nstep))
  cbuoy = N.zeros((nstep))
  cbuoy2 = N.zeros((nstep))
  
  for n in N.arange(nstep):  
    up = u[n,:]
    vp = v[n,:]
    wp = w[n,:]
    xp = x[n,:]
    yp = y[n,:]
    zp = z[n,:]
    bp = b[n,:]

    circ[n] = N.sum((up[0:ntraj-1] + up[1:ntraj]) * (xp[1:ntraj] - xp[0:ntraj-1])/2.0 \
                  + (vp[0:ntraj-1] + vp[1:ntraj]) * (yp[1:ntraj] - yp[0:ntraj-1])/2.0 \
                  + (wp[0:ntraj-1] + wp[1:ntraj]) * (zp[1:ntraj] - zp[0:ntraj-1])/2.0 )
      
    buoy[n] = N.sum((bp[0:ntraj-1] + bp[1:ntraj]) * (zp[1:ntraj] - zp[0:ntraj-1])/2.0 )
        
  
  cbuoy[nstep-1]  = circ[nstep-1]
  cbuoy2[nstep-1] = circ2[nstep-1]
  
  for n in N.arange(nstep-2,-1,-1):
    cbuoy[n]  = cbuoy[n+1]  + (t[n] - t[n+1]) * (buoy[n+1]+buoy[n])/2.0
    cbuoy2[n] = cbuoy2[n+1] + (t[n] - t[n+1]) * (buoy2[n+1]+buoy2[n]/2.0)
    
  fig = P.figure(figsize=(16,16))
  ax  = fig.add_subplot(111)

#   lines = P.plot(t, circ, t, cbuoy, t, circ2, t, cbuoy2)
#   P.setp(lines[0], color="k",  linewidth = 2.0)
#   P.setp(lines[1], color="b",  linewidth = 2.0)
#   P.setp(lines[2], color="g",  linewidth = 2.0)
#   P.setp(lines[3], color="r",  linewidth = 2.0)
#   P.legend( ('C-Trap', 'CircB-Trap', 'C-Simpson', 'CircB-Simpson') ,loc = 'upper left')

  lines = P.plot(t, circ2, t, cbuoy2)
  P.setp(lines[0], color="k",  linewidth = 2.0)
  P.setp(lines[1], color="b",  linewidth = 2.0)
  P.legend( ('C-Raw', 'B-CircB') ,loc = 'upper left')

  P.xlabel('Trajectory Time (sec)', fontsize=10)
  P.ylabel('Circulation m**2/sec', fontsize=10)
  P.ylim(0.0, 60000.)

  P.show()
