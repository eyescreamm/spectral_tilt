# tilt~
"tilt~" is real-time spectral tilt manipulation object for Pure-data. 
 
# DEMO
 (image)
 
# Requirement
If you use "tilt~" object for Pure-data, you can compile it with "makefile".
After saving files in my repository, compile it using the "make" command.
 
# Usage
Once you have compiled it, you can use tilt~ on pure-data. open pure-data and create an object. object name is "tilt~". Set the argument to the number of poles and zeros. You can also specify the number of poles and zeros in the third inlet. The first inlet will contain the mono audio signal and the second inlet will contain the supectral tilt value from -6 to 6. The outlets will output mono audio with an arbitrary spectral tilt depending on the object.

(image)
 
# Author
 
* Syuto Yoshida
 
# License
 
"tilt~" is under [GNU General Public License v3.0]
Please see LICENSE file for details.
