#!/bin/bash

echo "200 continuous chars" > run.log

for i in {1..100}
do
  echo "** tree $i **" >> run.log
  
  # get the i-th complete birth-death tree
  sed -n "$i"p  ../../simulator/bd.trees > bd.tre
  
  # generate a data file 
  ../../msim -i bd.tre -o data.nex -n 200 -d 2 >> run.log
  
  # save the true branch lengths
  ../../mb cmd_temp.nex >> run.log  
  grep "con_50_majrule" data.nex.con.tre >> true.brl.log
  
  # and infer the branch lengths
  ../../mb cmd_fixt.nex >> run.log  
  grep "con_50_majrule" data.nex.con.tre >> estm.brl.log
  
  # and finally joint inference
  ../../mb cmd_cont.nex >> run.log
  
  # rename files to avoid overwriting
  mv data.nex        data_$i.nex
  mv data.nex.con.tre sim_$i.con.tre
done
