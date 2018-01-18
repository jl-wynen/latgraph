# latgraph
Utilities for converting and manipulating files containing lattices.

This is a program to read and write lattices from/to different file formats and manipulate them in the process.
A file can be converted using
```Shell
python latgraph.py infile -o outfile
```
The file formats are deduced from the extensions; see list below.
Use `python latgraph.py -h` to get a list of supported arguments.

Some example files are included under [lattices](lattices):
- ``c60_ipr.w3d``, ``c60_ipr.w2d``: C60 fullerene obeying the isolated pentagon rule. Generated using [CaGe](https://caagt.ugent.be/CaGe/).
- ``c60_ipr_yml``: Like above but relabelled using the ``anticlockwise`` method. Contains some additional metadata.

## Requirements
- Python 3
- Numpy
- Matplotlib
- PyYaml

## Supported Formats
- writegraph3d [.w3d], writegraph2d [.w2d]
  - Adjacency graph
  - 3D/2D positions
  - *No* hopping strengths (defaults to 1)
  - *No* information about time
  - *No* metadata
- YAML [.yml/.yaml]
  - Adjacency graph
  - 3D/2D positions
  - Hopping strengths
  - Number of time slices
  - Name, comment metadata

## Relabelling
Graphs can be relabelled to change the order of sites.
For this, use the arguments ``-l file -m method`` to select a file to determine the labels from and a method.
The file must be consistent with the main input file but may have a different embedding (e.g. 2D vs. 3D).
Supported methods are
- ``innermost``: Start at the site clostest to the centre.
   Then move to the neighbour of the current site which is closest the the centre and has not been visited yet; repeat.
- ``anticlockwise``: Like 'innermost' but always move around the centre in anti-clockwise order.

Depending on the label file, the program can get stuck and will abort in that case.

## Plotting
Graphs can be drawn as 3D or 2D meshes.
Use ``-p`` to plot the main graph and ``-P`` to plot the label graph.
In addition, the adjacency matrix can be plotted using the argument ``-a``.
All plots show the state after relabelling if a labelfile was specified.
