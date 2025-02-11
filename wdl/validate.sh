#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l pav.wdl
java -jar ${WOMTOOL_PATH} validate -l TrgtMerge.wdl
java -jar ${WOMTOOL_PATH} validate -l Hifiasm.wdl
java -jar ${WOMTOOL_PATH} validate -l BcftoolsMergeSvSnp.wdl
java -jar ${WOMTOOL_PATH} validate -l Trgt.wdl
java -jar ${WOMTOOL_PATH} validate -l Panpop.wdl
java -jar ${WOMTOOL_PATH} validate -l JasmineIntersample1.wdl
java -jar ${WOMTOOL_PATH} validate -l JasminePrime.wdl
java -jar ${WOMTOOL_PATH} validate -l Jasmine.wdl
java -jar ${WOMTOOL_PATH} validate -l Truvari.wdl
java -jar ${WOMTOOL_PATH} validate -l PbsvIntersample.wdl
java -jar ${WOMTOOL_PATH} validate -l Pbsv.wdl
java -jar ${WOMTOOL_PATH} validate -l SnifflesIntersample.wdl
java -jar ${WOMTOOL_PATH} validate -l BcftoolsMergeIntersample.wdl
java -jar ${WOMTOOL_PATH} validate -l BcftoolsMergeIntrasample.wdl
java -jar ${WOMTOOL_PATH} validate -l Resolve.wdl
java -jar ${WOMTOOL_PATH} validate -l BedtoolsMergeDipcall.wdl
java -jar ${WOMTOOL_PATH} validate -l BcftoolsMergeDipcall.wdl
java -jar ${WOMTOOL_PATH} validate -l BcftoolsMergeNoDuplicates.wdl
java -jar ${WOMTOOL_PATH} validate -l Debreak.wdl
java -jar ${WOMTOOL_PATH} validate -l Svim.wdl
java -jar ${WOMTOOL_PATH} validate -l Nanovar.wdl
java -jar ${WOMTOOL_PATH} validate -l Cutesv.wdl
java -jar ${WOMTOOL_PATH} validate -l Sawfish.wdl
java -jar ${WOMTOOL_PATH} validate -l SVision.wdl
java -jar ${WOMTOOL_PATH} validate -l Sniffles.wdl
java -jar ${WOMTOOL_PATH} validate -l MapCCS.wdl
java -jar ${WOMTOOL_PATH} validate -l FilterDipcall.wdl
java -jar ${WOMTOOL_PATH} validate -l Subsample.wdl
