#!/usr/bin/python
import sys
import shutil
# set the path for preseq directory
if len(sys.argv) == 2:
	src_path = sys.argv[1]
else:
	sys.stderr.write('Please specify the path of preseq!\n')
	sys.exit(1)
# src_path is the path of preseq
src_path = src_path.strip()
if src_path[-1] != '/':
	src_path += '/'
# set the file names
cppfile = 'continued_fraction.cpp'
headerfile = 'continued_fraction.hpp'
# des_path is the path of preseqR/src
des_path = './src/'
# copy files
shutil.copy(src_path + cppfile, des_path + cppfile)
shutil.copy(src_path + headerfile, des_path + 'continued_fraction.h')
# read all lines
f = open(des_path + cppfile, 'r')
lines = f.readlines()
f.close()
# replace the line '#include "continued_fraction.hpp"' with "continued_fraction.h"'
f = open(des_path + cppfile, 'w')
# set the keywords to target the line
keywords = ('#include', '"continued_fraction.hpp"', '#include"continued_fraction.cpp"')
for line in lines:
	words = line.strip().split()
	# check whether the line is #include "continued_fraction.hpp"
	if keywords[0] in words and keywords[1] in words:
		f.write('#include "continued_fraction.h"\n')
	# no space between #include and "continued_fraction.hpp"
	elif keywords[2] in words:
		f.write('#include "continued_fraction.h"\n')
	else:
		f.write(line)

f.close()
sys.exit(0)
