# Distance/Occlusion Library + Tool

<a href="data/tiger.png"            ><img src="data/tiger.png"             width="192"/></a>
<a href="images/tiger-sdf.png"      ><img src="images/tiger-sdf.png"       width="192"/></a>
<a href="images/tiger-sdf-rg.png"   ><img src="images/tiger-sdf-rg.png"    width="192"/></a>
<a href="images/tiger-sdf-delta.png"><img src="images/tiger-sdf-delta.png" width="192"/></a>

<a href="data/kiwi.png"             ><img src="data/kiwi.png"              width="192"/></a>
<a href="images/kiwi-sdf.png"       ><img src="images/kiwi-sdf.png"        width="192"/></a>
<a href="images/kiwi-sdf-rg.png"    ><img src="images/kiwi-sdf-rg.png"     width="192"/></a>
<a href="images/kiwi-sdf-delta.png" ><img src="images/kiwi-sdf-delta.png"  width="192"/></a>

From left to right: original, signed distance with zero at 0.5, red/green SDF,
delta vectors to closest boundary point.

## Intro

The `DistanceFields.*` and `OcclusionFields.*` files here implement a number of algorithms for

* Distance field generation
* Signed distance field generation
* Scalar and directional occlusion generation

for both 2d (images) and 3d (volumes). Bitmask inputs are supported to keep
memory overhead down, particularly important for volumes.

The included demo tool, `distocc`, will process 2D images or 3D .obj files, and
can output the results in a number of different formats. By default the tool
uses the 'Danielsson' algorithm, and signed border initialisation, in order to
generate accurate results as quickly as possible, but other algorithms can be be
selected via -m (see below). These include:

    FastSweep        O(n), most accurate delta technique
    Danielsson       O(n), ~30% faster than FastSweep, decent accuracy
    JumpFlood        O(n log n), GPU-oriented algorithm, included to check accuracy
    BruteForce       O(n^2), exact, for checking other algorithms at low n.
    Chamfer          O(n), has issues on diagonals. Cheapest but quality suffers.
    Felzenszwalb     O(n) but much larger overhead than the others. Accurate.

The first four algorithms provide deltas from cell to closest point, which can
be used to evaluate distance, but also to answer nearest-point queries, or
generate Voronoi diagrams. The last two are strictly distance-only. Finally, as
a comparison signed distances can also be generated via the classic two-pass
interior/exterior algorithm, but this is twice as slow and uses twice as much
memory as the default.

## Components

The library comes in three parts:

- `MaskUtilities.hpp/cpp`. Contains utilities for operating on 1-bit images and
  volumes, as well as a reasonably efficient 2D/3D triangle -> bitmask
  rasteriser which may be useful for other purposes.

- `DistanceFields.hpp/cpp`. Contains distance and signed-distance generation
  techniques for mask images/volumes. There is some naive OpenMP support for the
  JumpFlood and Felzenszwalb algorithms if you compile with -fopenmp.

- `OcclusionFields.hpp/cpp`. Contains some fast occlusion generation techniques
  for masks, including (coarse) directional occlusion. Also contains a simple
  cone cast helper, which can be useful for CPU-side visiblity checks.

There is also `VLMini.hpp`, a small 2D/3D vector library, which you can replace
with your own, image read/write support in `DOFImage.*`, and basic .obj file
read support in `MeshSupport.*`.


## Building and Running

To build the library and tool, run 'make', or on Windows include the files in
your visual studio project. The code has been tested under clang, gcc, and
Visual C. C++11 is required. (Not for anything fancier than initialisers
though.)

Then to see the tool options run

    ./distocc

which should give you something like:

    Options:
      -d: generate distance field. Methods: 0, Danielsson; 1, Fast Sweep; 2, Jump Flood; 3, Brute Force; 4, Chamfer; 5, Felzenszwalb
      -s: generate signed distance field in a single pass. Methods: as for -d but Chamfer/Felzenszwalb unsupported.
      -S: generate signed distance field via orthodox two distance field passes.
      -x: generate occlusion. Methods 0: standard, 1: no self-occlusion, 2: directional components only
      -w <dx> <dy>: swept shadow

      -m <int>: select method variant for above

      -w <w:int> [<h:int>]        : set dimensions of image source
      -W <w:int> [<h:int> <d:int>]: set dimensions of volume source
      -p [<count:int> <seed:int>] : add random points to image/volume (default)
      -b <sides:int>              : add test box with the given number of sides to image/volume

      -f <path>: specify input image
      -t <int> : specify 0-255 threshold for creating mask from image
      -r       : reverse so white is occupied rather than black

      -F <path>: use given obj file to define mask.

      -o <name>: set base name for output file(s)
      -i <int> : select image output variant: 0, greyscale; 1, deltas; 2, r=exterior/g=interior
      -v       : log detailed output
      -c       : run checks on output


## Examples

Try the following:

    ./distocc -d -p 100                         # Generate a DF for 100 random points
    ./distocc -d -W 16 -b 6                     # Generate a DF for a test box within a 16^3 volume
    ./distocc -s -f data/tiger.png              # Generate an SDF for the given image
    ./distocc -s -f data/tiger.png -p 20 -i 1   # Generate an SDF for the given image, plus 20 random points, and save as deltas to surface
    ./distocc -d -W 64 -F data/plant2.obj       # Generate a 64^3 distance field for the given 3D mesh

