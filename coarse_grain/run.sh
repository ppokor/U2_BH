#!/bin/bash

date 

INITTOP="init_beforeekvi_oldformat.top"
INITCRD="init_beforeekvi.dat"
EKVIMC="ekvi_mc"
EKVIMD="ekvi_md"
INPFILE="run_input"

OXDNA="/u/p/ppokorna/progs/oxDNA/build/bin/oxDNA"
REPLICAS=100

for i in $(seq -f "%02g" 50 $REPLICAS); do 

  dir="run$i"
  backup="run${i}_old"

  if [[ -d "$dir" ]]; then
    # If backup exists, remove it
    if [[ -d "$backup" ]]; then
      rm -rf "$backup"
    fi
    # Rename the current directory to backup
    echo "Backuping $dir to $backup"
    mv "$dir" "$backup"
  fi

  mkdir "$dir"
  echo "Entering directory: $dir"
  cd "$dir" || { echo "Failed to enter $dir. Skipping..."; continue; }

  cp ../${EKVIMC} ../${EKVIMD} ../${INPFILE}  .
  # set the input files for ekvilibration in two steps and for the MD
  awk '!/topology|conf_file|trajectory_file/' ${EKVIMC} > tmpfile && \
  echo -e "topology = ../${INITTOP}\nconf_file = ../${INITCRD}\ntrajectory_file = MC.dat" >> tmpfile && \
  mv tmpfile ${EKVIMC}
  awk '!/topology|conf_file|trajectory_file/' ${EKVIMD} > tmpfile && \
  echo -e "topology = ../${INITTOP}\nconf_file = MC.dat\ntrajectory_file = init.dat" >> tmpfile && \
  mv tmpfile ${EKVIMD}
  awk '!/topology|conf_file|trajectory_file/' ${INPFILE} > tmpfile && \
  echo -e "topology = ../${INITTOP}\nconf_file = init.dat\ntrajectory_file = trajectory.dat" >> tmpfile && \
  mv tmpfile ${INPFILE}
  
  ${OXDNA} ${EKVIMC}
  ${OXDNA} ${EKVIMD}

  ${OXDNA} ${INPFILE}
  
  echo "Finished in $dir. Moving to next..."
  cd ..
done

echo "All jobs completed!"
date 
