#!/bin/bash 
kma_index -i /database/senterica/senterica.fsa -o /database/senterica/senterica &> /database/senterica/senterica.log
mlst.py -i /test/test.fq.gz -o /test/ -s senterica -mp kma -x --quiet
file=/test/results_tab.tsv
DIFF=$(diff $file /test/test_results.tsv)
if [ "$DIFF" == "" ] && [ -s $file ] ;
   then     
   echo "TEST SUCCEEDED"; 
else
   echo "TEST FAILED";
fi
