#!/usr/bin/env bash

for d in `find . -type d -name 'mdc*[0-9]' -maxdepth 1`; do
    python -W ignore createsummary.py $d/$d.xml $d/coinc.xml 
done
