#!/bin/bash 
kma_index -i /database/senterica/senterica.fsa -o /database/senterica/senterica
mlst.py -i /test/test.fq.gz -o /test/ -s senterica -mp kma -x --quiet
DIFF=$(diff /test/results_tab.tsv /test/test_results.tsv)
if [ "$DIFF" == "" ] && [ "$?" == 0 ] ;
   then     
   echo "TEST SUCCEEDED"; 
else
   echo "TEST FAILED";
fi