deconrank
======
Package to deconvolute, rank & score precursors prior to fragmentation


Installation
------------
Either from the source or from conda on the tomnl channel

Command line
------------
::

    # python -m deconrank -i [in dir] -o [out dir] -p [polarity OPTIONAL] -w [list of weights: adduct,intensity,purity,clustn ]

    $ python -m deconrank -i /path/2/camera_out.csv -o /path/2/out_dir/ -w 0.3,0.3,0.2,0.2

    $ python -m deconrank --help



Developers & Contributors
-------------------------
 - Thomas N. Lawson (tnl495@bham.ac.uk) - `University of Birmingham (UK) <http://www.birmingham.ac.uk/index.aspx>`_
 - Martin R. Jones (m.r.Jones.1@bham.ac.uk) - `University of Birmingham (UK) <http://www.birmingham.ac.uk/index.aspx>`_
 - Ralf J. M. Weber (r.j.weber@bham.ac.uk) - `University of Birmingham (UK) <http://www.birmingham.ac.uk/index.aspx>`_
 



License
-------
Released under the GNU General Public License v3.0 (see `LICENSE file <https://github.com/computational-metabolomics/deconrank/blob/master/LICENSE>`_)
