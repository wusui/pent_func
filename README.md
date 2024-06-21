# Functional Programming version of Pentomino solver.

This set of python scripts finds all rectangular tilings using one set of pentominos ([Wikipedia link](https://en.wikipedia.org/wiki/Pentomino#Tiling_puzzle_(2D)))

### grow_tree.py

Generates a list of linked nodes which represents all possible pentomino orientations.

### solve.py

Solves all rectangular tilings using the data from grow_tree.py.  Returns a dict indexed by rectangle dimensions ('3x20', '4x15', '5x12', '6x10').  Each entry of the dict will consist of lists of layouts, where each solution is a string of tile square values with rows separated by new lines.

### tester.py

Calls solve.py and compares solutions against data scraped from the Internet to insure that all possible solutions have been found.  Thanks to Stephen Chapman at Isomer Design for providing this information  (https://isomerdesign.com/Pentomino/)
