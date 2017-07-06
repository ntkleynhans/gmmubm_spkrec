#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gmm.h"

void printUsage( void )
{
	cout << "gmmscore: help" << endl;
	cout << endl;
	cout << "-h,  --help\t\tThis message" << endl;
	cout << "-i,  --input\t\tInput model file or VQ codebook" << endl;
	cout << "-w,  --world\t\tNormalize score with this model" << endl;
	cout << "-l,  --list\t\tFile containing data file list" << endl;
	cout << "-t,  --modeltype\tInput init file type 1 = model, 2 = VQ codebook" << endl;
	cout << "-b,  --worldtype\tInput init file type 1 = model, 2 = VQ codebook" << endl;
	cout << "-m,  --mixture\t\tMixture number" << endl;
	cout << "-d,  --dimension\tFeature vector dimension" << endl;
	cout << "-v,  --vfloor\t\tVariance flooring constant" << endl;
	cout << "-n,  --number\t\tNumber of feature vectors to load" << endl;
	cout << "-r,  --results\t\tOutput results file" << endl;
	cout << "-g,  --tag\t\tAdd 'tag' before score" << endl;

	exit( -1 );
}

int main( int argc, char *argv[] )
{
	int nextOption;

	const char * shortOptions = "hi:w:l:t:b:m:d:v:n:r:g:";

	const struct option longOptions[] = {
	{ "help", 0, NULL, 'h' },
	{ "input", 1, NULL, 'i' },
	{ "world", 1, NULL, 'w' },
	{ "list", 1, NULL, 'l' },
	{ "modeltype", 1, NULL, 't' },
	{ "worldtype", 1, NULL, 'b' },
	{ "mixture", 1, NULL, 'm' },
	{ "dimension", 1, NULL, 'd' },
	{ "vfloor", 1, NULL, 'v' },
	{ "number", 1, NULL, 'n' },
	{ "results", 1, NULL, 'r' },
	{ "tag", 1, NULL, 'g' },
	{ NULL, 0, NULL, 0 }
	};

	string modelFile, listFile, worldFile, resFile, tag;
	unsigned int modeltype, worldtype, mixture, dimension, vectorNum = 1000;
	double vfloor = 0.1;
	unsigned int check = 0;

	do {
		nextOption = getopt_long( argc, argv, shortOptions, longOptions, NULL );

		switch( nextOption )
		{
			case 'i':
				modelFile = optarg;
				check += 1;
				break;

			case 'l':
				listFile = optarg;
				check += 2;
				break;

			case 't':
				modeltype = atoi( optarg );
				check += 4;
				break;

			case 'm':
				mixture = atoi( optarg );
				check += 8;
				break;

			case 'd':
				dimension = atoi( optarg );
				check += 16;
				break;

			case 'n':
				vectorNum = atoi( optarg );
				break;

			case 'v':
				vfloor = atof( optarg );
				break;

			case 'w':
				worldFile = optarg;
				check += 32;
				break;

			case 'b':
				worldtype = atoi( optarg );
				check += 64;
				break;

			case 'r':
				resFile = optarg;
				check += 128;
				break;

			case 'g':
				tag = optarg;
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
			cout << "-i, --input not set" << endl;
			testTransaction = true;
		}

		if( ! (check & 2) )
		{
			cout << "-l, --list not set" << endl;
			testTransaction = true;
		}

		if( ! (check & 4) )
		{
			cout << "-t, --inittype not set" << endl;
			testTransaction = true;
		}

		if( ! (check & 8) )
		{
			cout << "-m, --mixture not set" << endl;
			testTransaction = true;
		}

		if( ! (check & 16) )
		{
			cout << "-d, --dimension not set" << endl;
			testTransaction = true;
		}

		if( (check & 32) )
		{
			if( ! (check & 64) )
			{
				cout << "-b, --worldtype not set" << endl;
				testTransaction = true;
			}
		}

		if( ! (check & 128) )
		{
			cout << "-r, --results not set" << endl;
			testTransaction = true;
		}

		if( testTransaction == true )
		{
			printUsage();
		}
	}

	double LL = 0.0;

	ofstream Fresult( resFile.c_str(), ios_base::app );

	if( ! tag.empty() )
	{
		Fresult << tag << "\t";
	}

	if( worldFile.empty() )
	{
		GMM model( modelFile, modeltype, mixture, dimension, vfloor, vectorNum );
		LL = model.LogL( listFile );

		cout << "Model Score: " << LL << endl;
		cout << "Final Score: " << LL << endl;
		Fresult << LL << endl;
	}
	else
	{
		double WL = 0.0;

		GMM model( modelFile, modeltype, mixture, dimension, vfloor, vectorNum );
		LL = model.LogL( listFile );

		GMM world( worldFile, worldtype, mixture, dimension, vfloor, vectorNum );
		WL = world.LogL( listFile );

		cout << "Model Score: " << LL << endl;
		cout << "World Score: " << WL << endl;
		cout << "Final Score: " << LL-WL << endl;
		Fresult << LL-WL << endl;
	}

	Fresult.close();
	return 0;
}
