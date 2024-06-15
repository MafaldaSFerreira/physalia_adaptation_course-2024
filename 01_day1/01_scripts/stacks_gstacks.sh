#!/bin/bash
###stacks_gstacks.sh
cd
mkdir -p stacks
cd stacks
mkdir -p gstacks
gstacks -I ~/bamfiles -M ~/scripts/popmap_all_day1.txt -O gstacks -t 3
