#! /bin/bash
CORAL_DIR="/Users/glucia/Projects/CorAL"
cd "$CORAL_DIR/PHEMTO"
source "$CORAL_DIR/bin/phemto" "-1d" "input/phemto_input.dat"
python3 plot.py
cd "$CORAL_DIR"