#!/bin/bash

echo "both continuous and discrete chars" > run.log

for i in {1..100}
do
  echo "** tree $i **" >> run.log
  
  # get the i-th complete birth-death tree
  sed -n "$i"p  ../../simulator/bd.trees > bd.tre
  
  # generate a data file 
  ../../msim -i bd.tre -o data.nex -n 100 -d 0 >> run.log
  
  # run mrbayes to infer the parameters
  ../../mb cmd_both.nex >> run.log
  
  # rename files to avoid overwriting
  mv data.nex        data_$i.nex
  mv data.nex.con.tre sim_$i.con.tre
done
