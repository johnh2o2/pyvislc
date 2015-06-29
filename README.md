#Pyvislc: Visualize lighrcurves

* **PyVislc** is currently only setup to read in lightcurves that are formatted like HAT lightcurves. PyVislc will be more generalized in later versions. 

##Installation

* Run `python setup.py install`.

##Usage

* To visualize a lightcurve stored in `<file>`, run:

`python /path/to/vislc.py --file <file>`

* To visualize several lightcurves, store the list of filepaths in a separate file (`<listfile>`), and then run:

`python /path/to/vislc.py --list <listfile`

* You can optionally specify your own flags and the path of an output file for you to store the results; for example:

`python /path/to/vislc.py --list <listfile> --logfile <logfile> --flags flag1/keyboardshortcut1 flag2/keyboardshortcut2 ...`

