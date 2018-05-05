#!/bin/bash -v
rm -rvf build && mkdir -p build && cd build && cmake ../ && make && ./path_planning
