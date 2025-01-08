#!/bin/bash

INDEX=at429630

INCLUDE=(
  src/*
  test/*
  lib/*
  matmul.cpp
  CMakeLists.txt
)
ZIP_NAME=${INDEX}.zip

rm -f $ZIP_NAME

TEMP_ZIPPING_DIR=${INDEX}
mkdir -p $TEMP_ZIPPING_DIR
# shellcheck disable=SC2068
find ${INCLUDE[@]} -exec cp -r --parents {} $TEMP_ZIPPING_DIR \;

pushd report
mkdir -p build
pdflatex -output-directory=build report.tex 2>&1 > /dev/null
popd

cp report/build/report.pdf $TEMP_ZIPPING_DIR

zip -r $ZIP_NAME $TEMP_ZIPPING_DIR
rm -rf $TEMP_ZIPPING_DIR
