This is a cluster algorithm for the Phi^4-Theory in 4 dimension.

It is based on fermiQCD which needs to be downloaded from http://fermiqcd.net/fermiqcd/default/index. Just follow the instuctions and change the path to the fermiQCD directory in the Makefile in the main folder.

In the main folder you also find an example input file "input.in" and an input file with explainantions of the parameters "example.in".

After compiling run_cluster.cpp you can run the program with "./run_cluster -i infile.in" .

At the moment the observables can be computed are magnetisation, propagators, and correlators, which will be writtin in its own file with a name created directly from the program. 
The MonteCarlo parameters are NOT part of the fielname, so be carful not to override it, when changing this. Check the output for ASCII or Binary.

Have fun!
