BedGraphParser
==============

Simple &amp; fast bedGraph parser

This is simple parser for bed graph format: http://genome.ucsc.edu/goldenPath/help/bedgraph.html

* Input: bed graph file with no track headers
* Output: returns dictionary of chromosomes, and values as numpy array of score for each bin (generally 20)

# Install
Download and then run:
> python3 setup.py build_ext --inplace

or under MS Windows:
> python.exe setup.py build_ext --inplace --compiler=msvc --plat-name=win-amd64

(Tested with VS 2012: to use new version of msvc compiler you may have to modify msvc9compiler.py in distutils and specifly get_build_version and PLAT_TO_VCVARS as described in http://www.xavierdupre.fr/blog/2013-07-07_nojs.html)

# Usage
See example use_bg.py how to load specific bg file.
```python
import BedGraphReader 
import numpy as np

chrom_sizes = readChromeSize()  # initialize dictionary of chromosome to chromosome size
genome = BedGraphReader.load_bedgraph(bg_file, chrom_sizes)  # actually read bed graph file

# show number
print("Number of lines with score 1 in chr1: %i" %np.sum(genome['chr1']==1))
```

# Performance
We assume good formatted input file, from trusted resource and we don't check for validty. As such this script should read bed graph very fast.
For 510Mb file it runs in ~13 secs.
Use use_bg.py to compare naive python implementation to this script. 
Example output:

>Enter bed graph file name: UW.CD14_Primary_Cells.ChromatinAccessibility.RO_01727.DS17391.bedGraph         
>Comparing bedgraph reader to naive python parsing
>Starting bed graph reader
>finished... Now starting naive implementation
>***Bg: 13.576324	 Naive: 116.989906***