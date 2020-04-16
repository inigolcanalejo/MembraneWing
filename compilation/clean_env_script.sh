#!/bin/sh
# store USER, HOME and then completely clear environment
U=$USER
H=$HOME

for i in $(env | awk -F"=" '{print $1}') ;
do
    unset $i || true   # ignore unset fails
done

# restore
export USER=$U
export HOME=$H
export PATH=/usr/local/bin:/usr/bin:/bin:/sbin:/usr/sbin

echo shell environment follows:
env