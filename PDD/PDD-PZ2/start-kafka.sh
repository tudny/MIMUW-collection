#!/bin/bash

if [ -z "$KAFKA_VERSION" ]; then
  KAFKA_VERSION="3.7.0"
fi
if [ -z "$SCALA_VERSION" ]; then
  SCALA_VERSION="2.12"
fi

echo "Starting Kafka $KAFKA_VERSION with Scala $SCALA_VERSION"

KAFKA_DIR="kafka"
KAFKA_TGZ="kafka_$SCALA_VERSION-$KAFKA_VERSION.tgz"

if [ ! -d $KAFKA_DIR ] && [ ! -f $KAFKA_TGZ ]; then
  echo "Downloading Kafka $KAFKA_VERSION"
  wget "https://dlcdn.apache.org/kafka/$KAFKA_VERSION/$KAFKA_TGZ"
fi

if [ ! -d $KAFKA_DIR ] && [ -f $KAFKA_TGZ ]; then
  echo "Extracting Kafka $KAFKA_VERSION"
  mkdir -p $KAFKA_DIR
  tar -xzf $KAFKA_TGZ --strip-components=1 -C $KAFKA_DIR
fi

echo "Starting Kafka $KAFKA_VERSION"
rm -f $KAFKA_DIR/logs/*
echo "Stopping Kafka $KAFKA_VERSION"
if $KAFKA_DIR/bin/kafka-server-stop.sh; then
  echo "Kafka $KAFKA_VERSION stopped"
  sleep 5
fi

echo "Cleaning up Kafka $KAFKA_VERSION"
rm -rf /tmp/kraft-combined-logs
KAFKA_CLUSTER_ID="$($KAFKA_DIR/bin/kafka-storage.sh random-uuid)"
echo "Kafka Cluster ID: $KAFKA_CLUSTER_ID"

echo "Starting Kafka $KAFKA_VERSION"
$KAFKA_DIR/bin/kafka-storage.sh format -t $KAFKA_CLUSTER_ID -c $KAFKA_DIR/config/kraft/server.properties
echo "Starting Kafka $KAFKA_VERSION in $(pwd)/$KAFKA_DIR"
$KAFKA_DIR/bin/kafka-server-start.sh -daemon $KAFKA_DIR/config/kraft/server.properties
echo "Kafka $KAFKA_VERSION started"

for i in {1..5}; do
  echo -ne "\rWaiting for Kafka $KAFKA_VERSION to start $i/5"
  sleep 1
done

echo -e "\nKafka $KAFKA_VERSION started"
