## MinMaxGD-Tabu
Tabu Search algorithm for the  Min-Max Graph Drawing Problem, as originally described in "Tabu Search for Min-Max Edge Crossing in Graphs" 
by Tommaso Pastore, Anna Martínez-Gavara, Antonio Napoletano, Paola Festa, and Rafael Martí.
## How to Compile
To compile invoke a make command in the "Release" folder.
## Usage
To display the meaning of parameters and help, run the executable with the flag "-h"

For TABU SEARCH algorithm: 
  ./TabuGIT 

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
 
 EXAMPLE RUN: ./TabuGIT instance.txt tabu 5 1 no 10000 instance out-f1.txt out-f2.txt 60 20 10 1


OUTFILES:
The two out files have different formats.

-- The first one, "out-file-summary" reports: alpha max_cross cross_sum total_time time_to best.

-- The second file collects a complete report of the execution of the algorithm, reporting: alpha tenure max_cross current_time, and has a line for each time the algorithm improved the objective function value.

## License
The MIT License (MIT)

Copyright (c) 2019 Tommaso Pastore, Anna Martínez-Gavara, Antonio Napoletano, Paola Festa, and Rafael Martí

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Contact
tommaso.pastore@unina.it - feel free to contact me.
