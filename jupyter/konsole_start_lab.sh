#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${DIR}

PATH=/home/wester/node/node-v14.15.4-linux-x64/bin:$PATH

konsole --hold --workdir ${DIR} -e /home/wester/miniforge3/bin/jupyter lab


