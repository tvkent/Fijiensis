#!/bin/bash

/cap1/tyler.kent/Software/LDhelmet_v1.7/ldhelmet table_gen --num_threads 20 -t 0.01 -r 0.0 0.1 10.0 1.0 100.0 -c $1 -o $2
