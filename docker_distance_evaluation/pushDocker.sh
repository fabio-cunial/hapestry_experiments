#!/bin/bash
set -euxo

if [ $# -eq 0 ]; then
  echo "Need to provide tag for the docker image."
  exit 1
fi
TAG=$1

cp ../scripts/*.java .
cp ../scripts/*.py .
podman build --progress=plain -t fcunial/hapestry_distance_evaluation .
podman tag fcunial/hapestry_distance_evaluation fcunial/hapestry_distance_evaluation:${TAG}
podman push fcunial/hapestry_distance_evaluation:${TAG}
rm -f *.java *.class *.py
