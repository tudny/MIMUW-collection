#!/bin/bash

INDEX="at429630"
FILE_NAME="${INDEX}".tar.gz
USER_HOST=${INDEX}@students.mimuw.edu.pl
SSH_PATH='$HOME/PW/Cube/'

cd src/main/java || exit 1

tar -czvf "$FILE_NAME" concurrentcube

scp "$FILE_NAME" "${USER_HOST}:${SSH_PATH}"
scp Validate.java "${USER_HOST}:${SSH_PATH}"
scp validate.sh "${USER_HOST}:${SSH_PATH}"

ssh "$USER_HOST" "cd $SSH_PATH; rm ${INDEX} -rf; ./validate.sh $INDEX"
