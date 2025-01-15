#!/bin/bash
set -euxo

if [ $# -eq 0 ]; then
  echo "Need to provide tag for the docker image."
  exit 1
fi
TAG=$1

#cp ../scripts/*.java .
docker build --progress=plain -t fcunial/hapestry_experiments .
docker tag fcunial/hapestry_experiments fcunial/hapestry_experiments:${TAG}
docker push fcunial/hapestry_experiments:${TAG}
#rm -f *.java *.class
