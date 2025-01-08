#!/bin/bash
set -e

# First argument is bootstrap server
# All other are topics

BOOTSTRAP_SERVER=$1
shift

echo "Creating topics..."
for topic in "$@"; do
    kafka/bin/kafka-topics.sh --create --bootstrap-server $BOOTSTRAP_SERVER --topic $topic
done
