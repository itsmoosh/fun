For visual representation of spherical harmonics,
	especially linear combinations of spherical harmonics.

Open the Jupyter notebook for easy access!

This program is meant for either plotting deviations from
	spherical symmetry or for pure harmonics.

Requires:
	matplotlib
	numpy

Files:
	Spherical_harmonics_visualization.ipynb
	harmonics.py
	ylm_plot.py

ylm_plot contains these functions:
	plot_r()
	pure()
	make_plot()

All function arguments are optional. See below for a list of optional
	arguments and a description for each argument.

Author: Marshall 'Moosh' Styczinski

Last updated: July 31, 2018


########################################################################
########################################################################


plot_r():
	
	Wrapper for make_plot(). Easy plotting of deviations from spherical
	symmetry along spherical harmonics.

	Optional argument (default value)	# explanation

		R_km (2634.)		# Planetary radius (what r = 1.0
							# corresponds to in plots)

		d_shell_km (100)	# Shell thickness (deviations are from
							# R_km - d_shell_km, as this is where a
							# subsurface ocean will be for icy moons)

		a (1.)				# Y00 coefficient. aka s, monopole term.
		b_1 (0.)			# Y1-1 coefficient. aka py, dipole term.
		b0 (0.)				# Y10 coefficient. aka pz, dipole term.
		b1 (0.)				# Y11 coefficient. aka px, dipole term.
		c_2 (0.)			# Y2-2 coefficient. aka d, quadrupole term.
		c_1 (0.)			# Y2-1 coefficient. aka d, quadrupole term.
		c0 (0.)				# Y20 coefficient. aka d, quadrupole term.
		c1 (0.)				# Y21 coefficient. aka d, quadrupole term.
		c2 (0.)				# Y22 coefficient. aka d, quadrupole term.

		scaling (True)		# Whether to exaggerate radial deviation.
							# r_min is subtracted from r(theta,phi) and
							# The result is what's plotted. This simply
							# stretches the deviation along r, but does
							# not correctly handle surfaces dipping
							# below r_min.

		do_3d (True)		# Whether to plot in interactive 3D if True,
							# or in the x, y, and z = 0 planes if False.

		save (False)		# Whether to save a .png of the initial view
							# to the local directory.

		r_min (0.8)			# Fraction of R_km to scale from when
							# scaling = True.


########################################################################
########################################################################


pure():
	
	Wrapper for make_plot(). Easy plotting of pure spherical harmonics
	and linear combinations of spherical harmonics.
	Note that although the default arguments are all 0, this is only
	for convenience. Some coefficients should be set to non-zero values
	manually--if none are set, a default is chosen.

	Optional argument (default value)	# explanation

		a (0.)				# Y00 coefficient. aka s, monopole term.
		b_1 (0.)			# Y1-1 coefficient. aka py, dipole term.
		b0 (0.)				# Y10 coefficient. aka pz, dipole term.
		b1 (0.)				# Y11 coefficient. aka px, dipole term.
		c_2 (0.)			# Y2-2 coefficient. aka d, quadrupole term.
		c_1 (0.)			# Y2-1 coefficient. aka d, quadrupole term.
		c0 (0.)				# Y20 coefficient. aka d, quadrupole term.
		c1 (0.)				# Y21 coefficient. aka d, quadrupole term.
		c2 (0.)				# Y22 coefficient. aka d, quadrupole term.

		do_3d (True)		# Whether to plot in interactive 3D if True,
							# or in the x, y, and z = 0 planes if False.

		save (False)		# Whether to save a .png of the initial view
							# to the local directory.


########################################################################
########################################################################


make_plot():

	Generates 2D or 3D plots based on the input harmonic coefficients.
	Many optional arguments for versatility. R_km and d_shell are
	selected to serve as analogues for Ganymede.

	Optional argument (default value)	# explanation

	##############
	##	Primary	##
	##############

		R_km (2634.)		# Planetary radius (what r = 1.0
							# corresponds to in plots)

		d_shell_km (100)	# Shell thickness (deviations are from
							# R_km - d_shell_km, as this is where a
							# subsurface ocean will be for icy moons)

		a (1.)				# Y00 coefficient. aka s, monopole term.
		b_1 (0.)			# Y1-1 coefficient. aka py, dipole term.
		b0 (0.)				# Y10 coefficient. aka pz, dipole term.
		b1 (0.)				# Y11 coefficient. aka px, dipole term.
		c_2 (0.)			# Y2-2 coefficient. aka d, quadrupole term.
		c_1 (0.)			# Y2-1 coefficient. aka d, quadrupole term.
		c0 (0.)				# Y20 coefficient. aka d, quadrupole term.
		c1 (0.)				# Y21 coefficient. aka d, quadrupole term.
		c2 (0.)				# Y22 coefficient. aka d, quadrupole term.

		scaling (True)		# Whether to exaggerate radial deviation.
							# r_min is subtracted from r(theta,phi) and
							# The result is what's plotted. This simply
							# stretches the deviation along r, but does
							# not correctly handle surfaces dipping
							# below r_min.

		do_3d (True)		# Whether to plot in interactive 3D if True,
							# or in the x, y, and z = 0 planes if False.

		save (False)		# Whether to save a .png of the initial view
							# to the local directory.

		r_min (0.9)			# Fraction of R_km to scale from when
							# scaling = True.

	##########################
	##	Plotting constants	##
	##########################

		colormap ('seismic')	# Choice of color scheme for figures.
								# seismic is a diverging blue/red color
								# scale that gets darker for
								# larger values.

		dot_color ('black')		# Color of center indicator for 2D plots

		grid_color ('gray')		# Color of grid lines in 2D plots

		bndry_color ('green')	# Color of average shell inner
								# boundary indicator for 2D plots

		grid_thk (0.5)			# Grid line thickness in pt

		cbar_title ('$\Delta r$ (km)') # Title text for colorbar

		print_r_desc (True)		# Whether to print a qualitative
								# description for r in terms of
								# Ylm coefficients

		show_grid (True)		# Whether to draw polar grid lines,
								# for illustrative purposes

		show_bndry (True)		# Whether to draw a circle indicating
								# the average shell inner boundary depth

		r_steps (4)				# Number of divisions for drawing grids
								# with draw_grids

		ang_steps (4)			# Number of angular grid lines to draw
								# with draw_grids (ΔΦ = π/ang_steps)

		n_phis (1000)			# Number of points around
								# each 2D polar plot

		n_phis_3d (100)			# Number of divisions along circles in
								# theta and phi for 3D plots. Square
								# this number for the total number of
								# points plotted. Increasing will make
								# surfaces smoother but rapidly
								# increases memory and processing needs.

		fig_dpi (300)			# DPI to use for output png

		fig_size ((13.5,5.0))	# Figure window size in inches

		fig3d_size ((8.0,8.25))	# Figure window size for 3D images; size
								# of window is tied to aspect ratio, so
								# be careful when changing size here.

		cbar_pos ( [0.92,0.2,0.02,0.65] )	# Bottom-left corner [ x,y,
											# w,h ] for colorbar,
											# as fraction of fig_size

		cbar3d_pos ( [0.9,0.35,0.02,0.4] )	# Same for 3D; in a slightly
											# different location due to
											# plot size difference.

	######################
	##	Plot label text	##
	######################

	# ax_append strings are tacked on after '$x$', '$y$', '$z$':

		x_append (r'$\,(R_P)$, corotation $\rightarrow$')

		y_append (r'$\,(R_P)$, parent body $\rightarrow$')

		z_append (r'$\,(R_P)$, spin axis $\circlearrowleft$')

		xylabel ('$xy$-plane (northern hemi)')

		yzlabel ('$yz$-plane (leading hemi)')

		xzlabel ('$xz$-plane (antiparent hemi)')

		plot_title ('Spherical harmonics')
