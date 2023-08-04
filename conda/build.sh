#!/bin/bash

mkdir -p ${PREFIX}/bin

chmod +x mlst.py
cp mlst.py ${PREFIX}/bin/mlst.py

# copy script to download database
chmod +x ${RECIPE_DIR}/download-mlst-db.sh
cp ${RECIPE_DIR}/download-mlst-db.sh ${PREFIX}/bin/download-mlst-db.sh
