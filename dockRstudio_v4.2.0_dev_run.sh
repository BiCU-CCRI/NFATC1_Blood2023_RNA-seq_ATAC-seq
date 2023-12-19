#!/bin/bash
 
# Description: bash script to start docker Rstudio server
#    Rstudio docker image: ccribioinf/dockrstudio:4.2.0-v1
#    docker image is based on tidyverse: rocker/tidyverse:4.2.0
#
# Usage: bash dockRstudio_v4.2.0_run.sh
#
# Parameters:
#    new_port (-p ${new_port}:8787) Map TCP port 8787 (default for Rstudio server) in the container to port $new_port on the Docker host.
#    -e USERID; -e GROUPID - specifying user id and group id for correct file ownership
#    -e PASSWORD - password for the Rstudio server; default user is 'rstudio', but this can be changed using, for example, -e USER=$(whoami) 
#    --volume - binding volumes for persistent data storage
#    --workdir - specifying working directory; this may not work in Rstudio server!

new_port=
RENV_PATHS_CACHE_HOST=xxxx  # the path to an renv cache on the host (global) machine
RENV_PATHS_CACHE_CONTAINER=/renv_cache  # the path to the cache within the container
# potentially bind to local cache to speed up package installation

docker run --rm \
  -p ${new_port}:8787 \
  -e USERID=$(id -u) \
  -e GROUPID=$(id -g) \
  -e PASSWORD=test0 \
  -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" \
  -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}" \
  --volume=$(pwd):"/home/rstudio/workspace" \
  --workdir="/home/rstudio/workspace" \
  ccribioinf/dockrstudio:4.2.0-v1

