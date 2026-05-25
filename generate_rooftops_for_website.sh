#!/bin/bash
OUTDIR="-o $HOME/self/pi3/dds-web/polyh/roof_tops/off"
EDITDIR="-o dds_to_edit"
BASENAME="-b rt_"
POSITION="-r 105 -20 0"

supersemicupola() {
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 7 2 7 1 --use_half_roof -w
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 7 3 7 2 --use_half_roof -w
}

# 4 1:
rt_4_1() {
  # Note: the following couldn't be attached to the base according to the rules
  # - {4/2}
  # - {6/3}
  # - {7/3}
  # - {8/4}
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 3 1 4 1 -w --add_extra_triangles --no_base --add_extra_base
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 3 1 4 1 -w --add_extra_triangles --no_base --add_extra_base -i 1
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 4 1 4 1 -w --add_extra_triangles --no_base                       # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 4 1 4 1 -w --add_extra_triangles --no_base -i 1                  # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 5 1 4 1 -w --add_extra_triangles --no_base --add_extra_base      # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 5 1 4 1 -w --add_extra_triangles --no_base --add_extra_base -i 1 # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 5 2 4 1 -w --add_extra_triangles --no_base --add_extra_base      # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 5 2 4 1 -w --add_extra_triangles --no_base --add_extra_base -i 1 # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 6 1 4 1 -w --add_extra_triangles --no_base --add_extra_base      # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 6 1 4 1 -w --add_extra_triangles --no_base --add_extra_base -i 1 # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 6 2 4 1 -w --add_extra_triangles --no_base --add_extra_base      # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 6 2 4 1 -w --add_extra_triangles --no_base --add_extra_base -i 1 # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 7 1 4 1 -w --add_extra_triangles --no_base --add_extra_base      # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 7 1 4 1 -w --add_extra_triangles --no_base --add_extra_base -i 1 # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 7 2 4 1 -w --add_extra_triangles --no_base --add_extra_base      # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 7 2 4 1 -w --add_extra_triangles --no_base --add_extra_base -i 1 # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 8 1 4 1 -w --add_extra_triangles --no_base --add_extra_base      # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 8 1 4 1 -w --add_extra_triangles --no_base --add_extra_base -i 1 # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 8 2 4 1 -w --add_extra_triangles --no_base                       # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 8 2 4 1 -w --add_extra_triangles --no_base -i 1                  # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 8 3 4 1 -w --add_extra_triangles --no_base --add_extra_base      # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 8 3 4 1 -w --add_extra_triangles --no_base --add_extra_base -i 1 # s
  # in reality: extra triangles and double base
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 8 3 4 1 -w --add_extra_base -i 1 -t _alt
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 9 1 4 1 -w --add_extra_triangles --no_base --add_extra_base      # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 9 1 4 1 -w --add_extra_triangles --no_base --add_extra_base -i 1 # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 9 2 4 1 -w --add_extra_triangles --no_base --add_extra_base      # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 9 2 4 1 -w --add_extra_triangles --no_base --add_extra_base -i 1 # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 9 3 4 1 -w --add_extra_triangles --no_base --add_extra_base      # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 9 3 4 1 -w --add_extra_triangles --no_base --add_extra_base -i 1 # s
}

# THESE NEED SOME MANUAL EDITING
edits_4_1() {
  # Finish by using squares instead of triangles.
  python roof_top_polyhedra.py $POSITION $BASENAME $EDITDIR 3 1 4 1 -w --no_base -i 1 -t _alt

  # Finish by using squares instead of triangles.
  # Note that four faces come together in one edge, but these need to be split:
  # different ways are possible
  python roof_top_polyhedra.py $POSITION $BASENAME $EDITDIR 4 1 4 1 -w --add_extra_triangles --no_base -i 1 # s
}

rt_5_1() {
  # Cannot be closed:
  # - (3/1} opposite
  # No solutions for {5/2}, {7/3}, {8/3}, {9/4}, {10/4} {11/5} {12/5}
  # Only one solution for {8/2}
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 3 1 5 1 -w --use_half_roof --del_no_of_opposite_triangles 2 -i 1 # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 5 1 5 1 -w --use_half_roof # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 5 1 5 1 -w --use_half_roof -t "_alt" --del_no_of_base_triangles 5 --no_base # o

  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 9 3 5 1 -w --use_half_roof --del_no_of_opposite_triangles 2 -i 1 # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 10 2 5 1 -w --use_half_roof # o
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 10 2 5 1 -w --use_half_roof -t "_alt" --del_no_of_base_triangles 10 --no_base # o

  # COMPOUND of 2:
  python roof_top_polyhedra.py $POSITION $BASENAME $EDITDIR 6 2 5 1 -w --use_half_roof --del_no_of_opposite_triangles 4 -i 1 # s
  # COMPOUND of 3:
  python roof_top_polyhedra.py $POSITION $BASENAME $EDITDIR 12 4 5 1 -w --use_half_roof --del_no_of_opposite_triangles 9 -i 1 # s
}

# THESE NEED SOME MANUAL EDITING
edits_5_1() {
  # Add a pair of concave decagons that look like pentagrams:
  python roof_top_polyhedra.py $POSITION $BASENAME $EDITDIR 10 1 5 1 -w --add_extra_triangles --add_extra_base -i 1 # o
  python roof_top_polyhedra.py $POSITION $BASENAME $EDITDIR 10 1 5 1 -w --add_extra_triangles --del_no_of_opposite_triangles 10 -i 1 -t _alt # o
}

rt_5_2() {
  # It is nicer to give the opposite triangle parallel to the base the same colour as the base
  # but for the site it is more consistent to leave it this way
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 3 1 5 2 -w --use_half_roof --del_no_of_opposite_triangles 2
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 6 2 5 2 -w --use_half_roof --del_no_of_opposite_triangles 4
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 9 3 5 2 -w --use_half_roof --del_no_of_opposite_triangles 6
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 12 4 5 2 -w --use_half_roof --del_no_of_opposite_triangles 8

  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 5 2 5 2 -w -i 1 --use_half_roof --no_base --del_no_of_base_triangles 5  # s
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 10 4 5 2 -w -i 1 --use_half_roof --no_base --del_no_of_base_triangles 10

  # double edges, but still a polyhedron according to Grünbaum
  #python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 5 2 5 2 -w --use_half_roof --add_extra_triangles --no_base --add_extra_base --del_no_of_opposite_triangles 5 -i 1 # s
}

# THESE NEED SOME MANUAL EDITING
edits_5_2() {
  # Add extra faces: add two small non-regular intersection decagrams. It needs to cover the two pentagons in the centre twice,
  # go from outer vertex to second inner vertex.
  # Note that there are edges that have 4 faces, these need to be seen as seperate edges:
  # E.g. where just {5/2} meet and where the triangles meet.
  python roof_top_polyhedra.py $POSITION $BASENAME $EDITDIR 10 3 5 2 -w --allow_holes --add_extra_triangles --no_base --add_extra_base
  # then this makes more sense (add same decagrams) though technically not a roof top anymore
  python roof_top_polyhedra.py $POSITION $BASENAME $EDITDIR 10 3 5 2 -w --allow_holes --add_extra_triangles --no_base --del_no_of_opposite_triangles 10 -t _alt
}

# {6/1}, {6/2}, {6/3} nothing found

rt_7_x() {
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 7 2 7 1 -w --use_half_roof
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 14 4 7 1 -w --use_half_roof
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 7 3 7 2 -w --use_half_roof
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 14 6 7 2 -w --use_half_roof
}

# The failed case for 993
rt_9_nok() {
  python roof_top_polyhedra.py $POSITION $BASENAME $OUTDIR 9 2 9 1 --use_half_roof --add_extra_triangles -w
}

#supersemicupola

#rt_4_1
#edits_4_1

#rt_5_1
#edits_5_1

#rt_5_2
#edits_5_2

rt_7_x

#rt_9_nok
