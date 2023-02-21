#!/usr/bin/env bash

# limit the resources of the JVM
export NXF_OPTS="-Xms500M -Xmx8G"

nextflow run main.nf -resume -profile standard -params-file input_params.yaml