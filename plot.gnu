## This file is for generating plots from the resulted data of the simulation data in GNUPLOT
## Created by Jonathan Henning

## Set the output to a PNG file (comment to see directly in screen)
#set terminal pngcairo dashed size 800,800

# Remove the subtitle
unset key

set size square

# Frames of the simulation (you can get the number of lines of the resulted data file)
total_frames = 650

do for [line=0:total_frames:1] {
	## Comment the output to file if you want to see directly in the screen
	#set output 'pict/'.line.'.png'
	plot [-5e8:5e8] [-5e8:5e8] "< ./sep ".line." 3" using ($1):($2) with points pt 7 ps 0.2
	## If you want to see directly in the screen is better to pause or it will be very quick (unless the data file is huge)
	pause 0.1
}
