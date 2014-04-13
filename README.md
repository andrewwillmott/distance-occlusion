The DistanceOcclusionLib.* files here include various utilities for generating 2D and 3D distance fields and occlusion fields. There's more info here: http://www.andrewwillmott.com/tech-notes#TOC-Distance-Fields-and-Volume-AO-Generation, and there's a related talk here: http://www.andrewwillmott.com/talks/fast-generation-of-directional-occlusion-volumes

To build the full version, run

    clang++ -D USE_TGA=1 -D USE_OBJ=1 DistanceOcclusionLib.cpp DistanceOcclusionTest.cpp MeshSupport.cpp targa.c -o distocc

Then to see the tool options run

    ./distocc

which should give you something like:

    Options:
      -d: generate image distance field. Methods 0: Danielsson, 1: Sweep, 2: 2 x Danielsson, 3: 2 x Sweep
      -D: generate volume distance field. Methods 0: Danielsson, 1: Sweep, 2: 4 x Danielsson, 3: 4 x Sweep
      -o: generate image occlusion field
      -O: generate volume occlusion field. Methods 0: standard, 1: no self-occlusion
      -m <int>: select method variant

      -s w [h [d]]: set dimensions of source
      -x <int>: use random point image/volume with the given seed
      -b <int>: use test block image/volume with the given number of sides
      -f <filename>: use tga as input image
      -f <filename>: use obj to define input volume

      -l : log detailed output

E.g., 

    ./distocc -O -f plant2.obj 

should generate occlusion and directional occlusion slices for the given tree mesh in AOV-NNN.tga and AOVDir-NNN.tga respectively.
