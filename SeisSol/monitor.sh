#!/bin/bash
color() { GREP_COLOR=$1 grep --line-buffered --color=always '.*'; }
(tail -qf $1.out | color 32 &
tail -qf $1.err | color 31) | cat