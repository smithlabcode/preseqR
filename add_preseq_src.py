#!/usr/bin/env python

##  Copyright (C) 2014
##  University of Southern California, Andrew D. Smith, Chao Deng, Timothy Daley
##
##  Authors: Chao Deng
##
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.


### automatically copy and modify required files from preseq to
### preseqR change header file names from '*hpp' to '*.h' change '
### #include "header.hpp" ' to ' #include "header.h" ' in src files

import sys
import shutil
import os
import fileinput
import re

## check that the right number of arguments are specified
if len(sys.argv) != 3:
  sys.stderr.write('Usage: ./add_preseq_src.py [preseq path] [preseqR path]\n')
  sys.exit(1)

## set the paths for preseq and preseqR directory
preseqPath = os.path.abspath(sys.argv[1])
preseqRPath = os.path.abspath(sys.argv[2])

## check any file including keyWords in its name
## matchingFile contains file names with keywords 
filenames = os.listdir(preseqPath)
keyWords = "continued_fraction"
matchingFile = [i for i in filenames if keyWords in i]

## copy required files from preseq into preseqR/src directory
preseqRSrcPath = os.path.join(preseqRPath, 'src')

## catch .o files
exe = re.compile('.*\.o')

for i in matchingFile:
  
  ## exclude .o file when copying
  if not exe.match(i):
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
