# latgraph
Utilities for creating, converting, and manipulating files containing lattices.

This program can read and write lattices from/to different file formats and manipulate them in the process.
Alternatively, new lattices can be generated.
Use `python latgraph.py -h` to get a list of supported arguments.

Some example files are included under [lattices](lattices):
- ``c60_ipr.w3d``, ``c60_ipr.w2d``: C60 fullerene obeying the isolated pentagon rule. Generated using [CaGe](https://caagt.ugent.be/CaGe/).
- ``c60_ipr_yml``: Like above but relabelled using the ``anticlockwise`` method. Contains some additional metadata.

## Requirements
- Python 3
- Numpy
- Matplotlib
- PyYaml

## Generating
Lattices can be generated using `python latgraph.py --generate <gen>`, where `<gen>` is the name of a supported generator.
Currently, only `tube` is supported for generating carbon nano tubes.

All arguments specified after `--generate <gen>` are passed directly to the generator. Any of the base argments need to be specified before generate. Use `--generate <gen> -h` to get generator specific help.

### Tube
This generator creates carbon nano tubes. Call it via
```
python latgraph.py -o tube.yml --generate tube 2,3 5
```
The first set of arguments (`2,3`, no spaces!) is the chirality of the tube.
The second argument (`5`) is the number of unit cells along the tube.
Important optional arguments are
- `--bc_sh`, `--bc_t` for boundary conditions along circumference and tube, respectively. Supported values are `[o]pen` and `[p]eriodic`, default is `periodic`
- `--emb` for the desired embedding; can be `2d` or `3d`.

The tube generator can also be used to create nano ribbons by specifying a 2D embedding (and open boundary conditions in this example):
```
python latgraph.py -o ribbon.yml --generate tube 3,3 10 --bc_ch=o --bc_t=o --emb=2d --comment="Hey a ribbon!"
```
Note however that the size of the ribbon must be set via the cirality and the number of tube(!) unit cells, not the number of hexagons.

## Supported Formats
Formats are deduced from file suffixes, see brackets below.
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
