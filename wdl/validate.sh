#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l Cutesv.wdl
java -jar ${WOMTOOL_PATH} validate -l Sawfish.wdl
java -jar ${WOMTOOL_PATH} validate -l SVision.wdl
java -jar ${WOMTOOL_PATH} validate -l Pbsv.wdl
java -jar ${WOMTOOL_PATH} validate -l Sniffles.wdl
java -jar ${WOMTOOL_PATH} validate -l MapCCS.wdl
java -jar ${WOMTOOL_PATH} validate -l FilterDipcall.wdl
java -jar ${WOMTOOL_PATH} validate -l Subsample.wdl
