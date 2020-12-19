#!/bin/bash

set -e

mkdir -p plots

# do initial correction without any offset
if [ ! -e output/offsets.txt ]; then
    python ./scripts/correct_exposure_and_align.py \
        data/20x_timelapse_MMStack_Pos18.ome.tif \
        output/timelapse.tif \
        data/empty_field/
fi

# calculate the offset using a dumb, brute force method
if [ ! -e output/offsets.txt ]; then
    python ./scripts/calc_offsets_brute_force.py \
        output/timelapse.tif \
        output/offsets.txt
fi

# use the calculated offset to make a final corrected timelapse
python ./scripts/correct_exposure_and_align.py \
    data/20x_timelapse_MMStack_Pos18.ome.tif \
    output/timelapse.tif \
    data/empty_field/ \
    output/offsets.txt

# extract cell data
python ./scripts/extract_track_data_from_timelapse.py output/timelapse.tif tracks/
python ./scripts/parse_cell_divisions.py tracks/cell_divisions/

# make the plots
python ./scripts/plot_cell_trajs.py output/cell_data.pkl output/cell_div_frame_nums.pkl
