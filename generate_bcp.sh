#!/bin/sh
# m = 1
ns="3 5 6 7 8 9 10 11 12"
for n in $ns; do python bcp.py -wo out/bcp $n; done

# m = 2
ns="5 6 7 9 10 11 12"
for n in $ns; do python bcp.py -wo out/bcp $n -m 2; done

# m = 3
ns="7 8 9 10 11"
for n in $ns; do python bcp.py -wo out/bcp $n -m 3; done

# m = 4
ns="9 10 11 12"
for n in $ns; do python bcp.py -wo out/bcp $n -m 4; done

# m = 5
ns="11 12"
for n in $ns; do python bcp.py -wo out/bcp $n -m 5; done
