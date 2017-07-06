#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "speaker.h"

void printUsage( void )
{
	cout << "gmmtrain: help" << endl;
	cout << endl;
	cout << "-h,  --help\t\tThis message" << endl;
	cout << "-o,  --output\t\tOutput model file" << endl;
	cout << "-i,  --input\t\tInput model file or VQ codebook" << endl;
	cout << "-l,  --list\t\tFile containing data file list" << endl;
	cout << "-t,  --inittype\t\tInput init file type 1 = model, 2 = VQ codebook" << endl;
	cout << "-e,  --traintype\tTraining type type 1 = EM, 2 = MAP " << endl;
	cout << "-m,  --mixture\t\tMixture number" << endl;
	cout << "-d,  --dimension\tFeature vector dimension" << endl;
	cout << "-v,  --vfloor\t\tVariance flooring constant" << endl;
	cout << "-n,  --number\t\tNumber of feature vectors to load" << endl;
	cout << "-a,  --adapt\t\tAdapt model parameters 1 (weights), 2 (means), 4 (variances)" << endl;
	cout << "            \t\tAdd together for adaption combinations: 1+2+4 = 7 for adapting all parameters" << endl;
	cout << "-p,  --percent\t\tTermination percent (NO 100% multiplier)" << endl;
	cout << "-r,  --results\t\tOutput results to this file" << endl;
	cout << "-c   --cycle\t\tIteration number" << endl;

	exit( -1 );
}

int main( int argc, char *argv[] )
{
	int nextOption;

	const char * shortOptions = "ho:i:l:t:e:m:d:v:n:a:p:r:c:";

	const struct option longOptions[] = {
	{ "help", 0, NULL, 'h' },
	{ "output", 1, NULL, 'o' },
	{ "input", 1, NULL, 'i' },
	{ "list", 1, NULL, 'l' },
	{ "inittype", 1, NULL, 't' },
	{ "traintype", 1, NULL, 'e' },
	{ "mixture", 1, NULL, 'm' },
	{ "dimension", 1, NULL, 'd' },
	{ "vfloor", 1, NULL, 'v' },
	{ "number", 1, NULL, 'n' },
	{ "adapt", 1, NULL, 'a' },
	{ "percent", 1, NULL, 'p' },
	{ "results", 1, NULL, 'r' },
	{ "cycle", 1, NULL, 'c' },
	{ NULL, 0, NULL, 0 }
	};

	string outModelFile, inFile, listFile, resultsFile;
	unsigned int inittype, traintype, mixture, dimension, vectorNum = 1000, adaptOpt = 0, iteration = 20;
	double vfloor = 0.1, percent = 0.005;

	unsigned int check = 0;

	do {
		nextOption = getopt_long( argc, argv, shortOptions, longOptions, NULL );

		switch( nextOption )
		{
			case 'o':
				outModelFile = optarg;
				check += 1;
				break;

			case 'i':
				inFile = optarg;
				check += 2;
				break;

			case 'l':
				listFile = optarg;
				check += 4;
				break;

			case 't':
				inittype = atoi( optarg );
				check += 8;
				break;

			case 'e':
				traintype = atoi( optarg );
				check += 16;
				break;

			case 'm':
				mixture = atoi( optarg );
				check += 32;
				break;

			case 'd':
				dimension = atoi( optarg );
				check += 64;
				break;

			case 'n':
				vectorNum = atoi( optarg );
				break;

			case 'a':
				adaptOpt = atoi( optarg );
				check += 128;
				break;

			case 'v':
				vfloor = atof( optarg );
				break;

			case 'p':
				percent = atof( optarg );
				break;

			case 'r':
				resultsFile = optarg;
				check += 256;
				break;

			case 'c':
				iteration = atoi( optarg );
				break;

			case 'h':
				printUsage();

			case '?':
				printUsage();

			case -1:
				break;

			default:
				abort();
		}

	} while( nextOption != -1 );

	bool testTransaction = false;

	{
		if( ! (check & 1) )
		{
			cout << "-o, --output not set" << endl;
			testTransaction = true;
		}

		if( ! (check & 2) )
		{
			cout << "-i, --input not set" << endl;
			testTransaction = true;
		}

		if( ! (check & 4) )
		{
			cout << "-l, --list not set" << endl;
			testTransaction = true;
		}

		if( ! (check & 8) )
		{
			cout << "-t, --inittype not set" << endl;
			testTransaction = true;
		}

		if( ! (check & 16) )
		{
			cout << "-e, --traintype not set" << endl;
			testTransaction = true;
		}

		if( ! (check & 32) )
		{
			cout << "-m, --mixture not set" << endl;
			testTransaction = true;
		}

		if( ! (check & 64) )
		{
			cout << "-d, --dimension not set" << endl;
			testTransaction = true;
		}

		if( ! (check & 128) && traintype == 2 )
		{
			cout << "-a, --adapt not set" << endl;
			testTransaction = true;
		}

		if( ! (check & 256) )
		{
			cout << "-r, --results not set" << endl;
			testTransaction = true;
		}

		if( testTransaction == true )
		{
			printUsage();
		}
	}

	Speaker person( outModelFile, inFile, inittype, mixture, dimension, vfloor, vectorNum, resultsFile );

	double oldLL = 0.0, newLL = 0.0;
	unsigned int total = 0;

	do
	{
		oldLL = newLL;
		person.modifyModel( listFile, traintype, adaptOpt );
		newLL = person.LogL( listFile );

		cout << "LL\t" << newLL << endl;
		total++;

		if( total == iteration )
		{
			cout << "Max iteration reached" << endl;
			break;
		}
	} while( fabs((newLL-oldLL)/(newLL)) > percent );

	person.saveModel();

	return 0;
}
