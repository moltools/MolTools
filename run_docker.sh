#!/bin/bash

# Check if at least two arguments are provided
if [ $# -lt 5 ]; then
  echo "Usage: $0 <SERVER_PORT> <CLIENT_PORT> <BROWSER_PORT> <DATABASE_PORT> <VOLUME> [ -d ]"
  exit 1
fi

# Set the environment variables
export SERVER_PORT=$1
export CLIENT_PORT=$2
export BROWSER_PORT=$3
export DATABASE_PORT=$4
export VOLUME=$5

# Check if the optional '-d' flag is provided as the third argument
if [ $# -eq 6 ] && [ "$5" == "-d" ]; then
  # Run Docker Compose with the updated environment variables in detached mode
  SERVER_PORT=$1 CLIENT_PORT=$2 BROWSER_PORT=$3 DATABASE_PORT=$4 VOLUME=$5 docker-compose up --build --force-recreate -d
else
  # Run Docker Compose with the updated environment variables in the foreground
  SERVER_PORT=$1 CLIENT_PORT=$2 BROWSER_PORT=$3 DATABASE_PORT=$4 VOLUME=$5 docker-compose up --build --force-recreate
fi
