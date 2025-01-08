#!/bin/bash

echo "Starting data generator..."
. venv/bin/activate

echo "Running data generator..."
python data_generator.py "$@"

echo "Data generator finished."
