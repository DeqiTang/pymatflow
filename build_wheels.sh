#!/bin/bash
# Build and Distribute Pymatflow
# References
#   https://realpython.com/python-wheels/#different-types-of-wheels
#   https://github.com/pypa/manylinux
# Multi Platform Support

# Debian 9 stretch
DOCKER_IMAGE=quay.io/pypa/manylinux_2_24_x86_64
PLAT=manylinux_2_24_x86_64
docker pull "$DOCKER_IMAGE"
#docker run -v /path/to/pymatflow:/root/pymatflow  -it quay.io/pypa/manylinux_2_24_x86_64 bash
# using -i so that we can use Ctr+C to interrupt process inside docker
docker container run -it --rm \
      -e PLAT=$PLAT \
      -v "$(pwd)":/root/pymatflow \
      "$DOCKER_IMAGE" /root/pymatflow/build_wheels/build_wheels_debian.sh

#CentOS 7
DOCKER_IMAGE=quay.io/pypa/manylinux2014_x86_64
PLAT=manylinux2014_x86_64
docker pull "$DOCKER_IMAGE"
#docker run -v /path/to/pymatflow:/root/pymatflow  -it quay.io/pypa/manylinux2014_x86_64 bash
# using -i so that we can use Ctr + C to interrupt process inside docker
docker container run -it --rm \
      -e PLAT=$PLAT \
      -v "$(pwd)":/root/pymatflow \
      "$DOCKER_IMAGE" /root/pymatflow/build_wheels/build_wheels_centos7.sh
