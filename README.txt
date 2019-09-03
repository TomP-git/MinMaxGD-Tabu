**********************************Instructions to run the TABU SEARCH algorithm********************************

To display the meaning of parameters and help, run the executable with the flag "-h"

 For TABU SEARCH algorithm: 
  ./exec 
  argv[1]: instance.txt
  argv[2]: algorithm={tabu}
  argv[3]: tenure 
  argv[4]: alpha=[0:1] 
  argv[5]: post-processing={yes, no} 
  argv[6]: max-it={0,..,10000}
  argv[7]: instance (without .txt) 
  argv[8]: out-file-summary=<string>
  argv[9]: out-file-complete=<string>
 argv[10]: time-limit={0,..,3600}
 argv[11]: max-it without improvement before diversification
 argv[12]: max diversification moves
 argv[13]: random seed
 EXAMPLE RUN: ./exec instance.txt tabu 5 1 no 10000 instance out-f1.txt out-f2.txt 60 20 10 1


OUTFILES:
The two out files have different formats.

-- The first one, "out-file-summary" reports: 
   alpha max_cross cross_sum total_time time_to best.

-- The second file collects a complete report of the execution of the algorithm, reporting:
   alpha tenure max_cross current_time,
   
   and has a line for each time the algorithm improved the objective function value.
