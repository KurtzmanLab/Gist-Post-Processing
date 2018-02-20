#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <vector>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include "gistpp.h"

/*
    Created by: Steven Ramsey, PhD Candidate Biochemistry CUNY Graduate Center
    Contact information: vpsramsey@gmail.com
                    Room 1409 Science Hall Lehman College, 250 Bedford Park blvd Bronx NY, 10468
    Mentor: Tom Kurtzman

*/



using namespace std;

int main (int argc, char** argv) {
    int i = 0; string infile; string infile2; string outfile = "out.dx"; string operation; vector<string > options; double val1; double val2;
    /*
        No arguments presented will print out help
    */
    if (argc <= 1) {
        cerr << "\nUSAGE:\n\n"
        << "./gistpp -i infile -op operation [-i2 infile 2][-o outfile][-opt options]\n\n"
        << "where:\n\n"
        << "bracketed commands are either optional or required dependant on the operation to be performed\n\n"
        << "infile is a file containing the main data to be read, this can either be the gist output file or a dx file\n"
        << "operation is the required specification of a desired operation to be performed on the infile\n\n"
        << "[infile 2] is a second infile, necessary when using an operation which requires two input files (ex: add or clear2)\n"
        << "[outfile] is the desired name of an outfile if the operation performed produces an outfile. Default: out.dx\n"	
        << "[options] are various option flags which apply to specific operations\n\n"
	<< "To see a full list of operations run: ./gistpp -operations\n\n"
        << "To see a full list of options and the operations they apply to run: ./gistpp -options\n\n";
	dx Fin;
	Fin.citation();
        exit(0);
    }


    else if (argc == 2) {
        /*
            Check if user has requested help, this information will also be printed if no arguments used
        */
        if (!strcmp(argv[1], "-h")) {
        	cerr << "\nUSAGE:\n\n"
        	<< "./gistpp -i infile -op operation [-i2 infile 2][-o outfile][-opt options]\n\n"
        	<< "where:\n\n"
        	<< "bracketed commands are either optional or required dependant on the operation to be performed\n\n"
        	<< "infile is a file containing the main data to be read, this can either be the gist output file or a dx file\n"
        	<< "operation is the required specification of a desired operation to be performed on the infile\n\n"
        	<< "[infile 2] is a second infile, necessary when using an operation which requires two input files (ex: add or clear2)\n"
        	<< "[outfile] is the desired name of an outfile if the operation performed produces an outfile. Default: out.dx\n"		
        	<< "[options] are various option flags which apply to specific operations\n\n"
		<< "To see a full list of operations run: ./gistpp -operations\n\n"
        	<< "To see a full list of options and the operations they apply to run: ./gistpp -options\n\n";
		dx Fin;
		Fin.citation();
        	exit(0);
        }
        /*
            Check if user has requested information regarding the available operations
        */
        else if (!strcmp(argv[1], "-operations")) {
            cout << "\nFollowing is a list of the available operations provided in this code\n\n"
            << "\tThe format is as follows:\n\n"
            << "\t\t-op [operation][# infile][type infile][outfile][synopsis]\n\n"
            << "\twhere:\n\toperation is the command flag\n"
            << "\t#infile = number of infiles required (1 or 2)\n"
            << "\ttype infile is the required infile format(dx, gist, or pdb where: gist is the gist outfile)\n"
            << "\toutfile is whether or not the code requires a specified outfile (Y/N)\n"
            << "\tsynopsis is a brief explanation of the operation\n\n"
            << "\t\t -op group 1 dx N produces files which contain spatially grouped voxels with similar qualities\n"
            //<< "\t\t -op contour 2 dx Y produces an outfile which contains the values of the first provided the same voxel meets a criteria in the second file\n"
            << "\t\t -op filter1 1 dx Y produces an outfile which contains all 1's and 0's, 0 if the voxel does not fit a desired criteria in infile (ex: high energy), 1 if the voxel does \n"
            //<< "\t\t -op filter2 2 dx Y same as clear1 only requires the voxel to fit criteria in 2 dx files rather than simply 1 (ex: high energy and high g(O))\n"
            << "\t\t -op sasa 1 dx Y produces a dx file which contains 1's in each voxel which defines the solvent accessible surface area\n"
            << "\t\t -op sum 1 dx N prints in the command line the sum of all voxel quantities and the average voxel quanitity for the values found in the dx file provided\n"
            << "\t\t -op add 2 dx Y produces a dx file in which each voxel is the sum of the 2 corresponding voxels in each input file\n"
            << "\t\t -op sub 2 dx Y produces a dx file in which each voxel is the difference of the 2 corresponding voxels in each input file (file 1 - file 2)\n"
            << "\t\t -op div 2 dx Y produces a dx file in which each voxel is the quotient of the 2 corresponding voxels in each input file (file1/file2)\n"
            << "\t\t -op mult 2 dx Y produces a dx file in which each voxel is the product of the 2 corresponding voxels in each input file\n"
            << "\t\t -op addconst 1 dx Y produces a dx file in which each voxel is changed by the addition of a specified constant\n"
            << "\t\t -op multconst 1 dx Y produces a dx file in which each voxel is changed by the multiplication of a specified constant\n"
            << "\t\t -op defbp 2 dx+pdb Y produces a dx file in which each voxel is flagged as 1 when within a set distance of any ligand.pdb heavy atoms\n"
            << "\t\t -op printpdb 1 dx Y produces a pdb file containing hydrogen atoms at every voxel with the gist voxel data stored in occupancy\n"
            //<< "\t\t -op popstat 1 gist N prints in the command line statistics taken directly from the gist outfile\n"
            //<< "\t\t -op vdw 1 dx Y prints a new dx file which contains 1's in all voxels which are within the vdw sphere of the atoms\n"
	    //<< "\t\t -op histo 1 dx prints a histogram .dat file of voxel quantities. Must specify an outfile.\n"
	    //<< "\t\t -op printcol 1 dx prints a dx file as a single column for spreadsheet purposes. Must specify an outfile.\n"
            << "\t\t -op makedx gist+dx 2 prints a new dx file which represents the data stored in a specified column of the gist text output\n\n";
	    dx Fin;
	    Fin.citation();
            exit(0);
        }
        /*
            Check to see if the user has requested information about operation options
        */
        else if (!strcmp(argv[1], "-options")) {
            /*
	    cout << "\nFollowing is a list of the available options and which operations they correspond to\n\n"
            << "\tThe format is as follows:\n\n"
            << "\t\t -opt [option][##][operation][synopsis]\n\n"
            << "\twhere:\n\toption is the command flag\n"
            << "\t## indicates a numerical value must be added after the option\n"
            << "\toperation is the list of operations this option is required for\n"
            << "\tsynopsis is a brief summary of the options effect\n\n"
            << "\t\t -opt cutoff1 ## [filter1 filter2 defbp] this option specifies the cutoff value which is desired to be applied to infile or the desired distance in defbp\n"
            << "\t\t -opt cutoff2 ## [filter2 contour] this option specifies the cutoff value which is desired to be applied to infile2\n"
            << "\t\t -opt gt1 [filter1 filter2] this option specifies we are interested in values greater than cutoff1 in infile\n"
            << "\t\t -opt lt1 [filter1 filter2] this option specifies we are interested in values less than cutoff1 in infile\n"
            << "\t\t -opt gt2 [filter2 contour] this option specifies we are interested in values greater than cutoff2 in infile2\n"
            << "\t\t -opt lt2 [filter2 contour] this option specifies we are interested in values less than cutoff2 in infile2\n"
            << "\t\t -opt const ## [addconst multconst vdw makedx] this option specifies the constant to add or multiply by in the provided function\n\n";
            exit(0);
	    */
		/*
			Rewritten to reflect the new operations help output, will not include information regarding commented out text of 
		*/
		cout << "\nFollowing is a list of the available options and which operations they correspond to\n\n"
            	<< "\tThe format is as follows:\n\n"
            	<< "\t\t -opt [option][##][operation][synopsis]\n\n"
            	<< "\twhere:\n\toption is the command flag\n"
            	<< "\t## indicates a numerical value must be added after the option\n"
            	<< "\toperation is the list of operations this option is required for\n"
            	<< "\tsynopsis is a brief summary of the options effect\n\n"
            	<< "\t\t -opt cutoff1 ## [filter1] this option specifies the threshold value which is desired to be applied to infile in filter1\n"
            	//<< "\t\t -opt cutoff2 ## [filter2 contour] this option specifies the cutoff value which is desired to be applied to infile2\n"
            	<< "\t\t -opt gt1 [filter1] this option specifies we are interested in values greater than cutoff1 in infile\n"
            	<< "\t\t -opt lt1 [filter1] this option specifies we are interested in values less than cutoff1 in infile\n"
            	//<< "\t\t -opt gt2 [filter2 contour] this option specifies we are interested in values greater than cutoff2 in infile2\n"
            	//<< "\t\t -opt lt2 [filter2 contour] this option specifies we are interested in values less than cutoff2 in infile2\n"
            	<< "\t\t -opt const ## [addconst multconst makedx defbp] this option specifies the constant utilize in the provided function\n\n";
		dx Fin;
		Fin.citation();
            	exit(0);
        }
        else {
            cerr << "\nERROR: Code requires more arguments\n"
            << "For help run ./gistpp -h\n\n";
            exit(0);
        }
    }

    /*
        The following is to read in the arguments given by the user
    */
    else {
        while (i < argc) {
            if (!strcmp(argv[i], "-i")) {
                infile = argv[++i];
            }
            if (!strcmp(argv[i], "-i2")) {
                infile2 = argv[++i];
            }
            if (!strcmp(argv[i], "-op")) {
                operation = argv[++i];
            }
            if (!strcmp(argv[i], "-opt")) {
                options.push_back(argv[++i]);
                if (!strcmp(argv[i], "cutoff1") || !strcmp(argv[i], "const")) {val1 = atof(argv[++i]);}
                if (!strcmp(argv[i], "cutoff2")) {val2 = atof(argv[++i]);}
            }
            if (!strcmp(argv[i], "-o")) {
                outfile = argv[++i];
            }
            i++;
        }
		//cout << operation << endl;
    }
    /*
        Follows is the section in which the specified operation is determined by if statements. This could possibly be converted to an enumerator and run smoother.
    */
    if (!strcmp(operation.c_str(), "group")) {
        //do group
        if (infile.empty()) {
            cerr << "\nERROR: An infile must be specified in order to group voxels\n"
            << "For help run ./gistpp -h\n\n";
            exit(0);
        }

        cout << "Group of: " << infile << " will be run\n";
        dx ONE;
        ONE.readDx(infile);
        ONE.makeGroups();
    }
    else if (!strcmp(operation.c_str(), "histo")) {
		//do histogram
		if (infile.empty()) {
			cerr <<"\nERROR: An infile must be specified in order to evaluate a histogram\n"
			<<"For help run gistpp -h\n\n";
			exit(0);
		}
		if (!strcmp(outfile.c_str(), "out.dx")) {
			cout << "\nWARNING: An outfile should be specified for histogramming, preferably .txt or .dat\n"
			<< "By default the outfile will be histo.dat\n";
			outfile = "histo.dat";
		}
		
		cout << "Histogram of: " << infile << " will be written to: " << outfile << endl;
		dx ONE;
		ONE.readDx(infile);
		ONE.histogram(outfile);
		
    }
	else if (!strcmp(operation.c_str(), "printcol")) {
		if (infile.empty()) {
			cerr <<"\nERROR: An infile must be specified with flag -i to print a column\n"
			<< "For help run gistpp -h\n\n";
			exit(0);
		}
		if (!strcmp(outfile.c_str(), "out.dx")) {
			cout << "\nWARNING: An outfile should be specified for printng to, preferably .txt or .dat\n"
			<< "By default the outfile will be print.dat\n";
			outfile = "print.dat";
		}

		cout << "Column of: " << infile << " will be written to: " << outfile << endl;
		dx ONE;
		ONE.readDx(infile);
		ONE.printcol(outfile);
	}
    else if (!strcmp(operation.c_str(), "contour")) {
        //do contour
        bool test1 = false; //cutoff2
        bool test2 = false; //gt2 or lt2
        //double c2 = 0.0; //storage for cutoff value cutoff2
        char flg; //storage for gt or lt flag
        if (infile.empty()) {
            cerr << "\nERROR: An infile must be specified with flag -i in order to run contour\n"
            << "For help run ./gistpp -h\n\n";
            exit(0);
        }
        if (infile2.empty()) {
            cerr << "\nERROR: An infile 2 must be specified with flag -i2 in order to run contour\n"
            << "For help run ./gistpp -h\n\n";
            exit(0);
        }
        if (options.size() == 0) {
            cerr << "\nERROR: Must specify cutoff2 and gt2/lt2 for contour\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        for (int j = 0; j < options.size(); j++) {
            /*
                This for loop will read in the options to ensure the user
                has entered all the necessary options to run this operation
            */
            if (!strcmp(options[j].c_str(), "cutoff2")) {test1 = true; continue;}
            if (!strcmp(options[j].c_str(), "gt2") || !strcmp(options[j].c_str(), "lt2")) {test2 = true;}
            if (!strcmp(options[j].c_str(), "gt2")) {flg = 'g';}
            if (!strcmp(options[j].c_str(), "lt2")) {flg = 'l';}
        }
        if (test1 != true && test2 != true) {
            cerr << "\nERROR: Must specify cutoff2 and gt2/lt2 for contour\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }

        cout << "Contour of: " << infile << " based on voxels in: " << infile2 << " with cutoff: " << val2 << " to be written to: " << outfile << endl;
        dx ONE;
        dx TWO;
        ONE.readDx(infile);
        TWO.readDx(infile2);
        ONE.contour(TWO, val2, flg);
        ONE.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "filter1")) {
        //do clear1
        bool test1 = false; //cutoff1
        bool test2 = false; //gt1 or lt1
        char flg; //storage for gt1 or lt1

        if (infile.empty()) {
            cerr << "\nERROR: An infile must be specified in order to run filter1\n"
            << "For help run ./gistpp -h\n\n";
            exit(0);
        }
        if (options.size() == 0) {
            cerr << "\nERROR: Must specify cutoff1 and gt1/lt1 with -opt cutoff1 and -opt gt1/lt1 for filter1\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        for (int j = 0; j < options.size(); j++) {
            /*
                This for loop will read in the options to ensure the user
                has entered all the necessary options to run this operation
            */
            if (!strcmp(options[j].c_str(), "cutoff1")) {test1 = true; continue;}
            if (!strcmp(options[j].c_str(), "gt1") || !strcmp(options[j].c_str(), "lt1")) {test2 = true;}
            if (!strcmp(options[j].c_str(), "gt1")) {flg = 'g';}
            if (!strcmp(options[j].c_str(), "lt1")) {flg = 'l';}
        }
        if (test1 != true && test2 != true) {
            cerr << "\nERROR: Must specify cutoff1 and gt1/lt1 for filter1\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }

        cout << "Filter1 to be run on: " << infile << " with cutoff " << val1 << " to be written to: " << outfile << endl;
        dx ONE;
        ONE.readDx(infile);
        dx TWO = ONE.clearByOne(val1, flg);
        TWO.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "filter2")) {
        //do clear2
        bool test1 = false; //cutoff1
        bool test2 = false; //cutoff2
        bool test3 = false; //gt1/lt1
        bool test4 = false; //gt2/lt2
        char flg; //storage for gt1/lt1
        char flg2; //storage for gt2/lt2

        if (infile.empty() || infile2.empty()) {
            cerr << "\nERROR: Two infiles need to be specified for filter2\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        if (options.size() == 0) {
            cerr << "\nERROR: Must specify cutoff1, cutoff2, gt1/lt1, gt2/lt2, each with the flag -opt to run filter2\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        for (int j = 0; j < options.size(); j++) {
            /*
                This for loop will read in the options to ensure the user
                has entered all the necessary options to run this operation
            */
            if (!strcmp(options[j].c_str(), "cutoff1")) {test1 = true; continue;}
            if (!strcmp(options[j].c_str(), "cutoff2")) {test2 = true; continue;}
            if (!strcmp(options[j].c_str(), "gt1") || !strcmp(options[j].c_str(), "lt1")) {test3 = true;}
            if (!strcmp(options[j].c_str(), "gt1")) {flg = 'g';}
            if (!strcmp(options[j].c_str(), "lt1")) {flg = 'l';}
            if (!strcmp(options[j].c_str(), "gt2") || !strcmp(options[j].c_str(), "lt2")) {test4 = true;}
            if (!strcmp(options[j].c_str(), "gt2")) {flg2 = 'g';}
            if (!strcmp(options[j].c_str(), "lt2")) {flg2 = 'l';}
        }
        if (test1 != true && test2 != true && test3 != true && test4 != true) {
            cerr << "\nERROR: Must specify cutoff1, cutoff2, gt1/lt1, gt2/lt2 to run filter2\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }

        cout << "Filter2 to be run on: " << infile << " with cutoff " << val1 << " and cutoff2" << val2 << " applied to: " << infile2 << " to be written to: " << outfile << endl;
        dx ONE; dx TWO;
        ONE.readDx(infile);
        TWO.readDx(infile2);
        dx THREE = ONE.clearByTwo(TWO, val1, val2, flg, flg2);
        THREE.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "sasa")) {
        //do sasa
        if (infile.empty()) {
            cerr << "\nERROR: An infile must be specified (gist-gO.dx) in order to run sasa\n"
            << "For help run ./gistpp -h\n\n";
            exit(0);
        }

        cout << "Solvent accessible surface area to be produced based on: " << infile << " to be written to: " << outfile << endl;
        cout << "Warning!\t SASA requires a g(O) file as input, cannot be run on any other dx file\n\n";
        dx ONE;
        ONE.readDx(infile);
        dx TWO = ONE.solventAccessible();
        TWO.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "sum")) {
        //do sum
        if (infile.empty()) {
            cerr << "\nERROR: An infile must be specified with flag -i in order to run sum\n"
            << "For help run ./gistpp -h\n\n";
            exit(0);
        }

        //cout << "The sum of all voxels in: " << infile << " to be calculated\n";
        dx ONE;
        ONE.readDx(infile);
        double sum = ONE.sum();
        double avg = ONE.avg();
        //cout << "sum of dx file: " << sum << endl;
	cout << "sum of: " << infile << " is: " << sum << endl; 
        cout << "avg of: " << infile << " is: " << avg << endl;
    }
    else if (!strcmp(operation.c_str(), "add")) {
        //do add
        if (infile.empty() || infile2.empty()) {
            cerr << "\nERROR: Two infiles need to be specified for add\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }

        cout << infile << " to be added to " << infile2 << " and written to: " << outfile << endl;
        dx ONE;
        dx TWO;
        ONE.readDx(infile);
        TWO.readDx(infile2);
        ONE.add(TWO);
        ONE.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "sub")) {
        //do sub
        if (infile.empty() || infile2.empty()) {
            cerr << "\nERROR: Two infiles need to be specified with flag -i for sub\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }

        cout << infile2 << " to be subtracted from " << infile << " and written to: " << outfile << endl;
        dx ONE;
        dx TWO;
        ONE.readDx(infile);
        TWO.readDx(infile2);
        ONE.sub(TWO);
        ONE.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "div")) {
        //do div
        if (infile.empty() || infile2.empty()) {
            cerr << "\nERROR: Two infiles need to be specified for div\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }

        cout << infile << " to be divided by " << infile2 << " and written to: " << outfile << endl;
        dx ONE;
        dx TWO;
        ONE.readDx(infile);
        TWO.readDx(infile2);
        ONE.div(TWO);
        ONE.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "mult")) {
        //do mult
        if (infile.empty() || infile2.empty()) {
            cerr << "\nERROR: Two infiles need to be specified with flags -i and -i2 for mult\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }

        cout << infile << " to be multiplied by " << infile2 << " and written to: " << outfile << endl;
        dx ONE;
        dx TWO;
        ONE.readDx(infile);
        TWO.readDx(infile2);
        ONE.mult(TWO);
        ONE.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "addconst")) {
        //do addbyconst
        bool test1 = false; //const

        if (infile.empty()) {
            cerr << "\nERROR: An infile needs to be specified with flag -i for addconst\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        if (options.size() == 0) {
            cerr << "\nERROR: Need to specify constant to add by with flag -opt const:\n\n"
            << "-opt const ##\n\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        for (int j = 0; j < options.size(); j++) {
            if (!strcmp(options[j].c_str(), "const" )) {test1 = true; break;}
        }
        if (test1 == false) {
            cerr << "\nERROR: Need to specify constant to add by with flag -opt const\n\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }

        cout << "Voxels in: " << infile << " to have constant " << val1 << " added to to them, then new values written to: " << outfile << endl;
        dx ONE;
        ONE.readDx(infile);
        ONE.addbyConst(val1);
        ONE.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "multconst")) {
        //do multbyconst
        bool test1 = false; //const

        if (infile.empty()) {
            cerr << "\nERROR: An infile needs to be specified with flag -i for multconst\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        if (options.size() == 0) {
            cerr << "\nERROR: Need to specify constant to mult by with -opt const\n\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        for (int j = 0; j < options.size(); j++) {
            if (!strcmp(options[j].c_str(), "const" )) {test1 = true; break;}
        }
        if (test1 == false) {
            cerr << "\nERROR: Need to specify constant to mult by\n\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }

        cout << "Voxels in: " << infile << " to have constant " << val1 << " multiplied to to them, then new values written to: " << outfile << endl;
        dx ONE;
        ONE.readDx(infile);
        ONE.multbyConst(val1);
        ONE.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "defbp")) {
        //do defbp
        //distance to be specified by const
        bool test1 = false; //cutoff1
        if (infile.empty() || infile2.empty()) {
            cerr << "\nERROR: Two infiles need to be specified for defbp, a dx and a ligand file\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
	if (infile.substr( infile.length() - 4 ) == ".pdb" ) {
		cerr << "\nERROR: Dx file expected to be defined by flag -i, .pdb file given instead\n"
		<< "For help run ./gistpp -h\n\n";
		exit(0);
	}
	if (infile2.substr( infile.length() - 3 ) == ".dx" ) {
		cerr << "\nERROR: Ligand file expected to be defined by flag -i2, .dx file given instead\n"
		<< "For help run ./gistpp -h\n\n";
		exit(0);
	}
        if (options.size() == 0) {
            cerr << "\nERROR: Need to specify desired distance around heavy atoms with const option with the flag -opt const\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        for (int j = 0; j < options.size(); j++) {
            if (!strcmp(options[j].c_str(), "const")) {test1 = true; break;}
        }
        if (test1 == false) {
            cerr << "\nERROR: Need to specify desired distance around heavy atoms with const option with the flag -otp const\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }


        cout << "Binding pocket to be defined using heavy atoms in: " << infile2 << " will be applied to: " << infile << " with a distance of: " << val1 << " and written to: " << outfile << endl;
        lig L;
        L.readLF(infile2);
        dx ONE;
        ONE.readDx(infile);
        ONE.setBP(L.hax, L.hay, L.haz, val1);
        ONE.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "heavi")) {
        //do defbp
        //distance to be specified by const
        bool test1 = false; //cutoff1
        if (infile.empty() || infile2.empty()) {
            cerr << "\nERROR: Two infiles need to be specified for heavi, a dx map specified with -i, and a ligand file with -i2\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        if (options.size() == 0) {
            cerr << "\nERROR: Need to specify desired distance around heavy atoms with const option using flag -opt const\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        for (int j = 0; j < options.size(); j++) {
            if (!strcmp(options[j].c_str(), "const")) {test1 = true; break;}
        }
        if (test1 == false) {
            cerr << "\nERROR: Need to specify desired distance around heavy atoms with const option\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }


        cout << "Heavi pocket to be defined using heavy atoms in: " << infile2 << " will be applied to: " << infile << " with a distance of: " << val1 << " and written to: " << outfile << endl;
        cout << "The ligand structure file is expected in i2 NOT in i!!\n";
        lig L;
        L.readLF(infile2);
        dx ONE;
        ONE.readDx(infile);
        ONE.setHeavi(L.hax, L.hay, L.haz, val1);
        ONE.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "popstat")) {
        //do pop
        if (infile.empty()) {
            cerr << "\nERROR: Two infiles need to be specified for popstat\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }

        cout << "popstat to be run on: " << infile << endl;

        pop ONE;
        ONE.calcpop(infile);
    }
    else if (!strcmp(operation.c_str(), "printpdb")) {
        //do print pdb
        if (infile.empty()) {
            cerr << "\nERROR: An infile must be specified for printpdb\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }

        cout << "A pdb file based on: " << infile << " will be written to: " << outfile << endl;
        dx ONE;
        ONE.readDx(infile);
        ONE.printPdb(outfile);
    }

    else if (!strcmp(operation.c_str(), "makedx")) {
        bool test1 = false; //cutoff1
        if (infile.empty()) {
            cerr << "\nERROR: The gist text file must be specified as -i for makedx\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        if (infile2.empty()) {
            cerr << "\nERROR: A dx file from the same system must be entered in -i2 for makedx\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        if (options.size() == 0) {
            cerr << "\nERROR: Need to specify desired column to be printed as dx file, with -opt const ##\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        for (int j = 0; j < options.size(); j++) {
            if (!strcmp(options[j].c_str(), "const")) {test1 = true; break;}
        }
        if (test1 == false) {
            cerr << "\nERROR: Need to specify desired column to be printed as dx file, with -opt const ##\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }


        cout << "A dx file will be written with column: " << val1 << " from the text file: " << infile << " and written to: " << outfile << endl;
        cout << "Keep in mind columns start at 0\n";

        dx ONE;
        ONE.readDx(infile2);
        ONE.write_out_dx(infile, val1);
        ONE.writeDx(outfile);
    }

    else if (!strcmp(operation.c_str(), "vdw")) {
        bool test1 = false; //cutoff1
        if (infile.empty()) {
            cerr << "\nERROR: The group dx file for a desired vdw volume analysis must be specified with -i for vdw\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        if (options.size() == 0) {
            cerr << "\nERROR: Need to specify the vdw volume desired (1.2 for water oxygen), with -opt const ##\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }
        for (int j = 0; j < options.size(); j++) {
            if (!strcmp(options[j].c_str(), "const")) {test1 = true; break;}
        }
        if (test1 == false) {
            cerr << "\nERROR: Need to specify the vdw volume desired (1.2 for water oxygen), with -opt const ##\n"
            << "For help run ./gistpp -h\n\n";
            exit (0);
        }

        dx ONE;
        ONE.readDx(infile);
        ONE.calcVdw(val1);
        ONE.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "cat")) {
        if (infile.empty() || infile2.empty()) {
            cerr << "\nERROR: The function cat attempts to combine two .dx outputs, therefore input files must be declared with -i and -i2\n"
                    << "For help run ./gistpp -h\n\n";
            exit(0);
        }
        dx ONE;
        ONE.readDx(infile);
        dx TWO;
        TWO.readDx(infile2);
        ONE.cat(TWO);
        ONE.writeDx(outfile);
    }
    else if (!strcmp(operation.c_str(), "zero")) {
        if (infile.empty()) {
        	cerr << "\nERROR: The function zero write a .dx map containing only 0.0's and requires and infile specified with the flag -i\n"
			<< "For help run ./gistpp -h\n\n";
		exit(0);
	}
	dx ONE;
	ONE.readDx(infile);
	dx TWO = ONE.zeros();
	TWO.writeDx(outfile);
    }
    else {
        cerr << "\nERROR: No operation specified to be run\n"
        << "For help run ./gistpp -h\n\n";
        exit (0);
    }
}
