#!/usr/bin/python

### automatically copy and modify required files from preseq to preseqR
### change header file names from '*hpp' to '*.h'
### change ' #include "header.hpp" ' to ' #include "header.h" ' in src files

import sys
import shutil
import os
import fileinput
import re

## set the paths for preseq and preseqR directory
if len(sys.argv) == 3:
  preseqPath = sys.argv[1]
  preseqRPath = sys.argv[2]
else:
  sys.stderr.write('Usage: ./add_preseq_src.py [preseq path] [preseqR path]\n')
  sys.exit(1)

## check any file including keyWords in its name
## matchingFile contains file names with keywords 
filenames = os.listdir(preseqPath)
keyWords = "continued_fraction"
matchingFile = [i for i in filenames if keyWords in i]

## copy required files from preseq into preseqR/src directory
preseqRSrcPath = os.path.join(preseqRPath, 'src')
for i in matchingFile:
  srcFile = os.path.join(preseqPath, i)
  desFile = os.path.join(preseqRSrcPath, i)
  shutil.copy(srcFile, desFile)

## change header suffix from 'hpp' with 'h'
os.chdir(preseqRSrcPath)
filenames = os.listdir('./')
headerFeature = re.compile('.+\.hpp')
headerFile = []
for i in filenames:

  ## check whether a file is a header file with '.hpp'
  if headerFeature.match(i):
    os.rename(i, i.replace('.hpp', '.h'))
    headerFile.append(i)

## modify #include statement since the name of the header changes
## #include "*.hpp" ----> #include "*.h"
srcFeature = re.compile('.+\.cpp')
lineFeature = [re.compile('\#include.*' + i + '.*') for i in headerFile]

for i in filenames:

  ## check whether a file is a src file with '.cpp'
  if srcFeature.match(i):
    for line in fileinput.input(i, inplace = True):
      
      ## check whether src code include a header which has changed name
      ## if does, rewrite the line to include the new header
      if any(feature.match(line) for feature in lineFeature):
        line = line.replace('.hpp', '.h')
      sys.stdout.write(line)

sys.exit(0)
