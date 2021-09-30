To start the ciruclar checkerboard stimulation with Landolt C task, execute
the following two steps:

Add COGENT (\Toolbox) to path 
http://www.vislab.ucl.ac.uk/cogent.php  
e.g. 
  addpath('C:\Cogent2000v1.32\Toolbox')

Then start
    circcheck_run_highc_demo_direct('SBJNAME',1,0)
where 1 is the number of the current run (each 5:30 mins, with 100 images)
, and 0 means keyboard input (1 means serial port input, in this case set 
portnum in file).

If multiple runs should be shown, call 
  circcheck_run_highc_demo_direct('SBJNAME',1,0)
  circcheck_run_highc_demo_direct('SBJNAME',2,0)
    etc. currently up to 
  circcheck_run_highc_demo_direct('SBJNAME',10,0)

Startup might take a bit. Log files will be created in \logs. Press any key
to start after loading. The files demo_direct and demo_direct_ramcache 
allow to show the images directly during loading (if
"show_images_during_startup = 1" is set in the scripts).

The script runs correctly if you a fast flickering ciruclar checkerboard 
(10Hz, currently hardcoded) with a Landolt C task in the middle. Flickering 
exchanegs white to black and vice versa. 

In the beginning, the checkboard has low contrast. Then 100 flickering
checkerboards with different highlighted sectors (medium and high contrast)
are shown, and finally the uniform low contrast checkboard again. If the 
script does not flicker, reduce the number of images (num_images) or use a 
computer with a better graphics card (with more memory).

The script has been tested with versions Cogent2000v1.32 and 2000v1.29 on 
Windows XP, Matlab2010b and Cogent2000v1.32, Windows 8, Matlab 2015a and 
R2011b.

If you experience problems with loading the images, this might be due to 
insufficient memory on the graphics card. Try the function 
  circcheck_run_highc_demo_direct('SBJNAME',1,0)
instead.

On old computers, that are too slow to show the stimulation and still 
support a 256 color mode, you can try 
  circcheck_run_highc_demo_pal('SBJNAME',1,0)
instead. NOTE: This demo lacks much of the logging capacities of the direct
demo. If the palette mode stimulation should be used, the logging needs
to be adapted.

Courtesy to Jakob Heinzle and John-Dylan Haynes.
Kai, 2015/08/03