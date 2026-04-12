#!/bin/bash
set -euxo

if [ $# -eq 0 ]; then
  echo "Need to provide tag for the docker image."
  exit 1
fi
TAG=$1

docker build --progress=plain -t fcunial/truvari_refine .
docker tag fcunial/truvari_refine fcunial/truvari_refine:${TAG}
docker push fcunial/truvari_refine:${TAG}
