"""
github.com/itsmoosh/fun/ylm_visual/ylm_plot.py

For visual representation of spherical harmonics,
	especially linear combinations of spherical harmonics.

This program is meant for plotting deviations from
	spherical symmetry, but pure spherical harmonics can be
	plotted by setting d_shell_km equal to R_km.

Requires:
	harmonics.py (for Ylm function definitions)
	matplotlib
	numpy

Functions:
	plot_r()
	pure()
	make_plot()

	All arguments are optional. For a list of optional arguments,
		see the README or look below.

Author: Marshall 'Moosh' Styczinski

Last updated: July 31, 2018
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d import axes3d
from matplotlib.patches import Circle
from matplotlib import cm
from matplotlib.collections import LineCollection

import harmonics as ylm


#############################################################
#	Plot r(theta,phi) for deviation from spherical symmetry	#
#############################################################

def plot_r(
  R_km = 2634.,	# Planetary radius
  d_shell_km = 100,	# Shell thickness
  a = 1., b_1 = 0., b0 = 0., b1 = 0.,
  c_2 = 0., c_1 = 0., c0 = 0., c1 = 0., c2 = 0.,
  scaling = True,	# Whether to exaggerate radial deviation
  do_3d = True,	# Whether to plot in interactive 3D
  save = False,	# Whether to save a .png of the initial view
  r_min = 0.8	# Center point radius for exaggerated scaling
  ):
	r_title = '$\,R_P$'

	# Turn scaling off when it will trash the plot
	if(d_shell_km > 0.1*R_km):
		if(scaling):
			print('Warning: scaling disabled because otherwise your plot will be garbled.')
			print('  Decrease shell thickness d or the scaling minimum radius r_min.')
		scaling = False

	if(scaling):
		rm = r_min
		units = ''	# Scaled distance is % from r_min to R_P, not % R_P
		title_start = 'Scaled s'
	else:
		rm = 0.0
		title_start = 'S'

	if(b_1 == 0 and b0 == 0 and b1 == 0 and c_2 == 0 and c_1 == 0 and c0 == 0 and c1 == 0 and c2 == 0):
		print('All harmonic coefficients zero. This will plot a uniform sphere.')
		print('  Specify a harmonic deviation to plot, e.g. ylm_plot.plot_r(c1=1.0)')

	plt_t = title_start + 'hell boundary, ' + str(rm) + r_title + ' to 1.0' + r_title
	if(do_3d): plt_t = plt_t + ',\n'	# 3D plot is narrower
	else:  plt_t = plt_t + ','
	plt_t = plt_t + '$\ d = $' + str(d_shell_km) + '$\,$km shell'

	make_plot( R_km, d_shell_km, a, b_1, b0, b1, c_2, c_1, c0, c1, c2, scaling,
	  do_3d, save, r_min=rm, plot_title=plt_t )


#################################################################
#	Plot r(theta,phi) for an arbitrary linear combination of	#
#		spherical harmonics up to degree 2						#
#################################################################

def pure(
  a = 0., b_1 = 0., b0 = 0., b1 = 0.,
  c_2 = 0., c_1 = 0., c0 = 0., c1 = 0., c2 = 0.,
  do_3d = True,	# Whether to plot in interactive 3D
  save = False	# Whether to save a .png of the initial view
  ):
	if(a == 0 and b_1 == 0 and b0 == 0 and b1 == 0 and c_2 == 0 and c_1 == 0 and c0 == 0 and c1 == 0 and c2 == 0):
		print('All coefficients zero. Specify a harmonic to plot, e.g. ylm_plot.pure(c2=1.0)')
		print('  Defaulting to Y22.')
		c2 = 1.0
	R_km = 1.
	d_shell_km = 1.
	scaling = False

	make_plot( R_km, d_shell_km, a, b_1, b0, b1, c_2, c_1, c0, c1, c2, scaling,
	  do_3d, save,
	  r_min = 0,
	  cbar_title='$r$',
	  x_append='',
	  y_append='',
	  z_append='',
	  xylabel='$xy$-plane',
	  yzlabel='$yz$-plane',
	  xzlabel='$xz$-plane',
	  plot_title = 'Real-valued spherical harmonics' )


#####################################################
#	Generate a plot of r(theta,phi) representing	#
#		a sum of spherical harmonics				#
#####################################################

def make_plot(
  R_km = 2634., d_shell_km = 100,
  a = 1., b_1 = 0., b0 = 0., b1 = 0.,
  c_2 = 0., c_1 = 0., c0 = 0., c1 = 0., c2 = 0.,
  scaling = True,	# Whether to exaggerate radial deviation
  do_3d = True,	# Whether to plot in interactive 3D
  save = False,	# Whether to save a .png of the initial view
  r_min = 0.9,	# Center point radius for exaggerated scaling
  colormap = 'seismic',	# Choice of color scheme for figures. seismic is a diverging blue/red color scale that gets darker for larger values.
  dot_color = 'black',		# Color of center indicator for 2D plots
  grid_color = 'gray',		# Color of grid lines in 2D plots
  bndry_color = 'green',	# Color of average shell inner boundary indicator for 2D plots
  grid_thk = 0.5,			# Grid line thickness in pt
  cbar_title = '$\Delta r$ (km)', # Title text for colorbar
  print_r_desc = True,		# Whether to print a qualitative description for r in terms of Ylm coefficients
  show_grid = True,		# Whether to draw polar grid lines, for illustrative purposes
  show_bndry = True,		# Whether to draw a circle indicating the average shell inner boundary depth
  r_steps = 4,				# Number of divisions for drawing grids with draw_grids
  ang_steps = 4,			# Number of angular grid lines to draw (ΔΦ = π/ang_steps)
  n_phis = 1000,			# Number of interpolated points around each polar plot
  n_phis_3d = 100,			# Square this for total # points; increase at cost to your machine (400+ looks good but runs slow)
  fig_dpi = 300,			# DPI to use for output png
  fig_size = (13.5,5.0),	# In inches
  fig3d_size = (8.0,8.25),	# Aspect ratio of axes3d objects can't be manually set due to a bug, so it's tied to window size.
  cbar_pos = [0.92,0.2,0.02,0.65],	# Bottom-left corner x,y, w,h for colorbar, in % of fig_size
  cbar3d_pos = [0.9,0.35,0.02,0.4],
  # Axis/plane labels
  x_append = r'$\,(R_P)$, corotation $\rightarrow$',
  y_append = r'$\,(R_P)$, parent body $\rightarrow$',
  z_append = r'$\,(R_P)$, spin axis $\circlearrowleft$',
  xylabel = '$xy$-plane (northern hemi)',
  yzlabel = '$yz$-plane (leading hemi)',
  xzlabel = '$xz$-plane (antiparent hemi)',
  plot_title = 'Spherical harmonics'
  ):

	d_shell_frac = d_shell_km/R_km
	xlabel = '$x$' + x_append
	ylabel = '$y$' + y_append
	zlabel = '$z$' + z_append
	rm = r_min

	#####################################
	#	Spherical harmonic coefficients	#
	#####################################
	#
	#	d/norm terms are chosen to limit shell thickness
	#	for a 2-term expansion (i.e., r = a*Y00 + c0*Y20) to be >= 0,
	#	for a pre-factor of 1.0. This is equivalent to a maximum
	#	deviation up to the planetary surface.
	#
	a = a * (1. - d_shell_frac) / ylm.l0_norm
	b_1 = b_1 * d_shell_frac / ylm.l1_norm
	b0 = b0 * d_shell_frac / ylm.l1_norm
	b1 = b1 * d_shell_frac / ylm.l1_norm
	c_2 = c_2 * d_shell_frac / ylm.l2_norm
	c_1 = c_1 * d_shell_frac / ylm.l2_norm
	c0 = c0 * d_shell_frac / (ylm.l2_norm * np.sqrt(12.))
	c1 = c1 * d_shell_frac / ylm.l2_norm
	c2 = c2 * d_shell_frac / ylm.l2_norm


	#####################
	#	Populate arrays	#
	#####################

	if(do_3d):
		# Switch to common numbers of points
		num_phis = n_phis_3d
		num_thetas = int(n_phis_3d/2)
	else:
		num_phis = n_phis
		num_thetas = int(n_phis/2)

	# Create 1D arrays for our gridded (theta,phi) values
	thetas_lin = np.linspace( 0.0, np.pi, num_thetas)
	phis_lin = np.linspace( 0.0, 2.0*np.pi, num_phis)

	if(do_3d):
		# Create 2D array of (theta,phi) values for making the surface
		thetas, phis = np.meshgrid(thetas_lin, phis_lin)
		# Define r grid using the 2D (theta,phi) grid (and subtract min value in case of scaling)
		r = a*ylm.Y00(thetas,phis) + b_1*ylm.Y1_1(thetas,phis) + b0*ylm.Y10(thetas,phis) + b1*ylm.Y11(thetas,phis) + c_2*ylm.Y2_2(thetas,phis) + c_1*ylm.Y2_1(thetas,phis) + c0*ylm.Y20(thetas,phis) + c1*ylm.Y21(thetas,phis) + c2*ylm.Y22(thetas,phis) - rm
		# Record deviations from spherical symmetry for coloration
		dr = (r+rm - (1.0-d_shell_frac))*R_km
		# Force r to be positive for benefit of plotting pure harmonics
		r = abs(r)
		# Find xyz for ease of data handling/plotting
		x = r*np.sin(thetas)*np.cos(phis) / (1.0-rm)
		y = r*np.sin(thetas)*np.sin(phis) / (1.0-rm)
		z = r*np.cos(thetas) / (1.0-rm)
	else:
		thetas, phis = thetas_lin, phis_lin

		# xy-plane: one specific theta value
		theta = np.pi/2.
		# Define r array at that theta, using array of phi values
		r = a*ylm.Y00(theta,phis) + b_1*ylm.Y1_1(theta,phis) + b0*ylm.Y10(theta,phis) + b1*ylm.Y11(theta,phis) + c_2*ylm.Y2_2(theta,phis) + c_1*ylm.Y2_1(theta,phis) + c0*ylm.Y20(theta,phis) + c1*ylm.Y21(theta,phis) + c2*ylm.Y22(theta,phis) - rm
		dr_xy = (r+rm - (1.0-d_shell_frac))*R_km
		r = abs(r)
		x_xy = r*np.sin(theta)*np.cos(phis) / (1.0-rm)
		y_xy = r*np.sin(theta)*np.sin(phis) / (1.0-rm)

		# yz-plane: 2 specific phi values
		phi = np.pi/2.
		# Define r array at this phi, using array of theta values
		r = a*ylm.Y00(thetas,phi) + b_1*ylm.Y1_1(thetas,phi) + b0*ylm.Y10(thetas,phi) + b1*ylm.Y11(thetas,phi) + c_2*ylm.Y2_2(thetas,phi) + c_1*ylm.Y2_1(thetas,phi) + c0*ylm.Y20(thetas,phi) + c1*ylm.Y21(thetas,phi) + c2*ylm.Y22(thetas,phi) - rm
		dr_yz1 = (r+rm - (1.0-d_shell_frac))*R_km
		r = abs(r)
		y_yz1 = r*np.sin(thetas)*np.sin(phi) / (1.0-rm)
		z_yz1 = r*np.cos(thetas) / (1.0-rm)

		phi = 3.*np.pi/2.
		# Reverse the order of theta array for smooth curve in r(theta) plot
		thetas = np.pi - thetas
		# Define r array at this phi, using array of theta values
		r = a*ylm.Y00(thetas,phi) + b_1*ylm.Y1_1(thetas,phi) + b0*ylm.Y10(thetas,phi) + b1*ylm.Y11(thetas,phi) + c_2*ylm.Y2_2(thetas,phi) + c_1*ylm.Y2_1(thetas,phi) + c0*ylm.Y20(thetas,phi) + c1*ylm.Y21(thetas,phi) + c2*ylm.Y22(thetas,phi) - rm
		dr_yz2 = (r+rm - (1.0-d_shell_frac))*R_km
		r = abs(r)
		y_yz2 = r*np.sin(thetas)*np.sin(phi) / (1.0-rm)
		z_yz2 = r*np.cos(thetas) / (1.0-rm)

		# Join arrays from each of the two phi values
		y_yz = np.concatenate( (y_yz1, y_yz2) )
		z_yz = np.concatenate( (z_yz1, z_yz2) )
		dr_yz = np.concatenate( (dr_yz1, dr_yz2) )

		# xz-plane: 2 specific phi values
		phi = np.pi
		# Define r array at this phi, using array of theta values
		r = a*ylm.Y00(thetas,phi) + b_1*ylm.Y1_1(thetas,phi) + b0*ylm.Y10(thetas,phi) + b1*ylm.Y11(thetas,phi) + c_2*ylm.Y2_2(thetas,phi) + c_1*ylm.Y2_1(thetas,phi) + c0*ylm.Y20(thetas,phi) + c1*ylm.Y21(thetas,phi) + c2*ylm.Y22(thetas,phi) - rm
		dr_xz1 = (r+rm - (1.0-d_shell_frac))*R_km
		r = abs(r)
		x_xz1 = r*np.sin(thetas)*np.cos(phi) / (1.0-rm)
		z_xz1 = r*np.cos(thetas) / (1.0-rm)

		phi = 0.0
		# Reverse the order of theta array again for smooth curve in r(theta) plot
		thetas = np.pi - thetas
		# Define r array at this phi, using array of theta values
		r = a*ylm.Y00(thetas,phi) + b_1*ylm.Y1_1(thetas,phi) + b0*ylm.Y10(thetas,phi) + b1*ylm.Y11(thetas,phi) + c_2*ylm.Y2_2(thetas,phi) + c_1*ylm.Y2_1(thetas,phi) + c0*ylm.Y20(thetas,phi) + c1*ylm.Y21(thetas,phi) + c2*ylm.Y22(thetas,phi) - rm
		dr_xz2 = (r+rm - (1.0-d_shell_frac))*R_km
		r = abs(r)
		x_xz2 = r*np.sin(thetas)*np.cos(phi) / (1.0-rm)
		z_xz2 = r*np.cos(thetas) / (1.0-rm)

		# Join arrays from each of the two phi values
		x_xz = np.concatenate( (x_xz1, x_xz2) )
		z_xz = np.concatenate( (z_xz1, z_xz2) )
		dr_xz = np.concatenate( (dr_xz1, dr_xz2) )


	#########################################
	#	Create graph objects -- 3D plots	#
	#########################################

	if(do_3d):
		fig = plt.figure(figsize=fig3d_size)
		ax = fig.gca(projection='3d')
		gp = 0.65	# Gray % to fill in axes backgrounds on 3D plots
		ax.w_xaxis.set_pane_color((gp+0.14,gp+0.14,gp+0.14,gp+0.14))
		ax.w_zaxis.set_pane_color((gp+0.07,gp+0.07,gp+0.07,gp+0.07))
		ax.w_yaxis.set_pane_color((gp,gp,gp,gp))
		ax.set_xlabel(xlabel)
		ax.set_ylabel(ylabel)
		ax.set_zlabel(zlabel)

		# Define colorbar to be symmetric since we are looking at diverging values
		c_lim = max( abs(dr.min()), abs(dr.max()) )
		cmap_got = plt.get_cmap(colormap)
		# Create an array from -1 to +1 for mapping later
		norm = plt.Normalize(-c_lim, c_lim)
		# Define light source direction for 3D shading
		ls = matplotlib.colors.LightSource(315, 45)
		# Map dr value for each point to specific colors in the colormap we chose
		cbar_map = cmap_got(norm(dr))

		# Generate the surface
		surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=cbar_map, linewidth=0, antialiased=False, shade=True)

		# Generate a colorbar to show deviation value scale
		m = matplotlib.cm.ScalarMappable(cmap=colormap)
		m.set_array( (cbar_map-0.5)*2*c_lim )
		cbar_ax = fig.add_axes(cbar3d_pos)
		cbar = fig.colorbar(m, ax=ax, cax=cbar_ax)
		cbar.ax.set_title(cbar_title,size=14, pad=10)


	#########################################
	#	Create graph objects -- 2D plots	#
	#########################################

	else:
		# Create a figure with 3 side-by-side subplots with defined spacing
		fig,sp_axes = plt.subplots(1,3,figsize=fig_size)
		fig.subplots_adjust(left=0.06, right=0.9, wspace=0.25, hspace=0.05, top=0.95, bottom=0.05)
		sp_axes = np.reshape(sp_axes,3)
		ax1, ax2, ax3 = sp_axes

		# Get colorbar limit as max deviation in any of the 3 planes
		c_lim = max( abs(dr_xy.min()), abs(dr_xy.max()), abs(dr_yz.min()), abs(dr_yz.max()), abs(dr_xz.min()), abs(dr_xz.max()) )
		norm = plt.Normalize(-c_lim, c_lim)

		# Configure 1D arrays into line segments we can recolor individually
		points_xy = np.array([x_xy, y_xy]).T.reshape(-1, 1, 2)
		segments_xy = np.concatenate([points_xy[:-1], points_xy[1:]], axis=1)
		lc_xy = LineCollection(segments_xy, cmap=colormap, norm=norm)
		lc_xy.set_array(dr_xy)
		xyplane = ax1.add_collection(lc_xy)

		# Same for yz
		points_yz = np.array([y_yz, z_yz]).T.reshape(-1, 1, 2)
		segments_yz = np.concatenate([points_yz[:-1], points_yz[1:]], axis=1)
		lc_yz = LineCollection(segments_yz, cmap=colormap, norm=norm)
		lc_yz.set_array(dr_yz)
		yzplane = ax2.add_collection(lc_yz)

		# Same for xz
		points_xz = np.array([x_xz, z_xz]).T.reshape(-1, 1, 2)
		segments_xz = np.concatenate([points_xz[:-1], points_xz[1:]], axis=1)
		lc_xz = LineCollection(segments_xz, cmap=colormap, norm=norm)
		lc_xz.set_array(dr_xz)
		xzplane = ax3.add_collection(lc_xz)

		# Generate colorbar to show deviation scale
		cbar_ax = fig.add_axes(cbar_pos)
		cbar = plt.colorbar(xyplane, ax=sp_axes, cax=cbar_ax)
		cbar.ax.set_title(cbar_title,size=14, pad=10)


	#############################
	#	Display and formatting	#
	#############################

	if(do_3d):
		ax.set_xlim(-1.0, 1.0)
		ax.set_ylim(-1.0, 1.0)
		ax.set_zlim(-1.0, 1.0)
	else:
		# Set aspect ratio and axes limits
		[ ax.set_aspect('equal') for ax in sp_axes ]
		[ ax.set_xlim(-1.1, 1.1) for ax in sp_axes ]
		[ ax.set_ylim(-1.1, 1.1) for ax in sp_axes ]
		# Add dot to indicate center location
		[ ax.scatter(0,0,color=dot_color,zorder=4) for ax in sp_axes ]

		# Label each plot and their axes
		ax1.set_title(xylabel,fontsize=16)
		ax2.set_title(yzlabel,fontsize=16)
		ax3.set_title(xzlabel,fontsize=16)
		ax1.set_xlabel(xlabel)
		ax1.set_ylabel(ylabel)
		ax2.set_xlabel(ylabel)
		ax2.set_ylabel(zlabel)
		ax3.set_xlabel(xlabel)
		ax3.set_ylabel(zlabel)

		if(show_grid):
			# Print a polar/radial grid
			for r_mark in np.arange(0.,1.,1./r_steps):
				# Radial grid indicators (circles)
				[ ax.add_artist(Circle((0.0,0.0), r_mark+1./r_steps, edgecolor=grid_color, fill=False, linewidth=grid_thk)) for ax in sp_axes ]
			for ang_mark in np.arange(0.,2.0*np.pi,np.pi/ang_steps):
				# Angular grid indicators (lines)
				[ ax.plot([-np.cos(ang_mark), np.cos(ang_mark)], [-np.sin(ang_mark), np.sin(ang_mark)], color=grid_color, linewidth=grid_thk) for ax in sp_axes ]
		else:
			# We're not drawing the main grid indicators, but it's still helpful
			#	to have something to guide the eye and judge deviations.
			# Draw a circle on each plot to indicate the planetary surface (unit circle)
			[ ax.add_artist(Circle((0.0,0.0), 1., edgecolor=grid_color, linewidth=grid_thk, fill=False)) for ax in sp_axes ]

		if(show_bndry):
			# Draw a circle on each plot to indicate the average shell thickness at depth d
			[ ax.add_artist(Circle((0.0,0.0), (1.-rm-d_shell_frac)/(1.-rm), edgecolor=bndry_color, linewidth=grid_thk, fill=False)) for ax in sp_axes ]		

	if(scaling):
		# Hide cartesian grid outlines/tick marks to avoid confusion
		if(do_3d):
			[s.set_visible(False) for s in ax.spines.values()]
			ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
		else:
			for axis in (ax1,ax2,ax3):
				[s.set_visible(False) for s in axis.spines.values()]
				axis.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
		
	# Add the title to the plot
	plt.suptitle( plot_title, fontsize=20, bbox=dict(facecolor='white') )

	if(print_r_desc):
		# Combine r definition, in terms of Ylms, into a human-readable expression
		pdiv = '$\quad+\quad$'
		ndiv = '$\quad-\quad$'
		pre  = '$\, d\cdot$'
		# Print l = 0 term
		if(abs(a) <= 0.001): a_desc = ''
		else: a_desc = '{:.3f}'.format( a*ylm.l0_norm ) + pre + '$Y_{00}$'
		# For l = 1 and l = 2 terms, print nothing if the term is 0.
		#	Print a + or - alone, then the coefficient and scaling, then Ylm.
		if(b_1 == 0):  b_1_desc = ''
		elif(b_1 < 0): b_1_desc = ndiv + '{:.3f}'.format( abs(b_1)*ylm.l1_norm ) + pre + '$Y_{1-1}$'
		elif(b_1 > 0): b_1_desc = pdiv + '{:.3f}'.format(     b_1 *ylm.l1_norm ) + pre + '$Y_{1-1}$'
		if(b0 == 0):   b0_desc  = ''
		elif(b0 < 0):  b0_desc  = ndiv + '{:.3f}'.format( abs(b0)*ylm.l1_norm )  + pre + '$Y_{10}$'
		elif(b0 > 0):  b0_desc  = pdiv + '{:.3f}'.format(     b0 *ylm.l1_norm )  + pre + '$Y_{10}$'
		if(b1 == 0):   b1_desc  = ''
		elif(b1 < 0):  b1_desc  = ndiv + '{:.3f}'.format( abs(b1)*ylm.l1_norm )  + pre + '$Y_{11}$'
		elif(b1 > 0):  b1_desc  = pdiv + '{:.3f}'.format(     b1 *ylm.l1_norm )  + pre + '$Y_{11}$'
		if(c_2 == 0):  c_2_desc = ''
		elif(c_2 < 0): c_2_desc = ndiv + '{:.3f}'.format( abs(c_2)*ylm.l2_norm ) + pre + '$Y_{2-2}$'
		elif(c_2 > 0): c_2_desc = pdiv + '{:.3f}'.format(     c_2 *ylm.l2_norm ) + pre + '$Y_{2-2}$'
		if(c_1 == 0):  c_1_desc = ''
		elif(c_1 < 0): c_1_desc = ndiv + '{:.3f}'.format( abs(c_1)*ylm.l2_norm ) + pre + '$Y_{2-1}$'
		elif(c_1 > 0): c_1_desc = pdiv + '{:.3f}'.format(     c_1 *ylm.l2_norm ) + pre + '$Y_{2-1}$'
		if(c0 == 0):   c0_desc  = ''
		elif(c0 < 0):  c0_desc  = ndiv + '{:.3f}'.format( abs(c0)*ylm.l2_norm*np.sqrt(12) ) + pre + '$Y_{20}$'
		elif(c0 > 0):  c0_desc  = pdiv + '{:.3f}'.format(     c0*ylm.l2_norm *np.sqrt(12) ) + pre + '$Y_{20}$'
		if(c1 == 0):   c1_desc  = ''
		elif(c1 < 0):  c1_desc  = ndiv + '{:.3f}'.format( abs(c1)*ylm.l2_norm ) + pre + '$Y_{21}$'
		elif(c1 > 0):  c1_desc  = pdiv + '{:.3f}'.format(     c1 *ylm.l2_norm ) + pre + '$Y_{21}$'
		if(c2 == 0):   c2_desc  = ''
		elif(c2 < 0):  c2_desc  = ndiv + '{:.3f}'.format( abs(c2)*ylm.l2_norm ) + pre + '$Y_{22}$'
		elif(c2 > 0):  c2_desc  = pdiv + '{:.3f}'.format(     c2 *ylm.l2_norm ) + pre + '$Y_{22}$'

		r_id = '$r\quad=\quad$' + a_desc + b_1_desc + b0_desc + b1_desc + c_2_desc + c_1_desc + c0_desc + c1_desc + c2_desc
		plt.figtext(0.5, 0.03, r_id, horizontalalignment='center', fontsize=16)

	# Optionally save, then display the results
	if(save):
		fig.savefig('ylm_output.png', format='png', dpi=fig_dpi)
	plt.show()
