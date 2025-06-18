#!/bin/bash

for d in nonclock*/*/
do  
    echo "working in $d"
    cd $d  
    cp ../../scripts/sum_t.r .
    
    # distance metrics
    echo "reso  sja  d_qt  d_ci  d_rf" > dist_t.txt
    for i in {1..100}
    do
    	sed -n "$i"p  ../../simulator/bd.trees > bd.tre
    	Rscript sum_t.r sim_$i.con.tre >> ../../out.log 2>&1
    done

    # tree lengths
    grep "      TL"      run.log > tl_estm.txt

    rm bd.tre
    rm sum_t.r
    cd ../..
done

for d in tipdating*/*/
do  
    echo "working in $d"
    cd $d  
    cp ../../scripts/sum_t.r .
    
    # distance metrics
    echo "reso  sja  d_qt  d_ci  d_rf" > dist_t.txt
    for i in {1..100}
    do
    	sed -n "$i"p  ../../simulator/bd.trees > bd.tre
    	Rscript sum_t.r sim_$i.con.tre >> ../../out.log 2>&1
    done

    # tree heights and clock rates
    grep "    age.*\[0\]" run.log > th_estm.txt
    grep "     clockrate" run.log > cl_base.txt
    grep "     ilnvar"    run.log > cl_var.txt

    rm bd.tre
    rm sum_t.r
    cd ../..
done
