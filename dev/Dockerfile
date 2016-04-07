# Dockerfile to create a development envirionnement for ptools.
#
# Author: Benoist LAURENT <benoist.laurent@gmail.com>
# Date: 2016/03/30
# Version: 1.0.0
#
# Usage:
#     $ cd /path/to/ptools/
#     $
#     $ # Docker container creation (first time only).
#     $ docker build -t ptools:dev ./dev
#     $
#     $ # Work in docker environment (from ptools local repository).
#     $ docker run --rm -v $(pwd):/src/ptools ptools:dev

FROM debian:latest

# Install system dependecies.
COPY requirements_system.txt /opt/
RUN apt-get update && cat /opt/requirements_system.txt | xargs apt-get install -y

# Install python dependencies
COPY requirements_python.txt /opt/
RUN pip install --upgrade pip && pip install -r /opt/requirements_python.txt

# The entrypoint script is executed everytime the container starts.
COPY entrypoint.sh /
ENTRYPOINT ["/entrypoint.sh"]