#include "gmm.h"

GMM::GMM( string modelInitFile, unsigned int initType, unsigned int mixtures, unsigned int length, double floor, unsigned int dataSize )
{
	MixtureNumber = mixtures;
	Dimension = length;
	vFloor = floor;
	MaxDataNumber = dataSize;

	means = new vector< valarray<double> * > (MixtureNumber);
	variances = new vector< valarray<double> * > (MixtureNumber);
	dataParm = new vector< valarray<double> * >(MaxDataNumber);

	weights = new valarray<double>( 0.0, MixtureNumber );
	globalvars = new valarray<double>( 0.0, Dimension );
	distance = new valarray<double>( 0.0, Dimension );
	PR = new valarray<double>( 0.0, MixtureNumber );

	int i = 0;

	while( i < MixtureNumber )
	{
		(*means)[i] = new valarray<double>( 0.0, Dimension );
		(*variances)[i] = new valarray<double>( 0.0, Dimension );
		i++;
	}

	i = 0;
	while( i < MaxDataNumber )
	{
		(*dataParm)[i++] = new valarray<double>( 0.0, Dimension );
	}

	if( initType == 1 )
	{
		loadModel( modelInitFile );
	}
	else if( initType == 2 )
	{
		loadVQ( modelInitFile );
	}
	else
	{
		InClassError( this, "GMM(): InitType error value specified not known.", -100 );
	}
}

void GMM::loadModel( string modelFile )
{
	Finit.open( modelFile.c_str(), ios_base::binary );

	Finit.read( reinterpret_cast <char *>( &Hmodel ), sizeof( ModelHeader ) );

	if( !Finit )
	{
		InClassError( this, "LoadModel(): Cannot open model file " + modelFile + " .", -200 );
	}

	if( Hmodel.MixtureNumber != MixtureNumber )
	{
		InClassError( this, "LoadModel(): Model file reports non-equal mixture number", -201 );
	}

	if( Hmodel.Dimension != Dimension )
	{
		InClassError( this, "LoadModel(): Model file reports non-equal feature vector dimension", -202 );
	}

	int i = 0;
	double value;

	while( i < Dimension )
	{
		Finit.read( reinterpret_cast <char *>( &value ), sizeof(double) );
		(*globalvars)[i++] = value;
	}

	int j = 0;
	i = 0;

	while( i < MixtureNumber )
	{
		Finit.read( reinterpret_cast <char *>( &value ), sizeof(double) );
		(*weights)[i] = value;

		while( j < Dimension )
		{
			Finit.read( reinterpret_cast <char *>( &value ), sizeof(double) );
			(*(*means)[i])[j] = value;

			Finit.read( reinterpret_cast <char *>( &value ), sizeof(double) );
			(*(*variances)[i])[j++] = value;
		}
		j = 0;
		i++;
	}
	
	Finit.close();
}

void GMM::loadVQ( string dataFile )
{
	Finit.open( dataFile.c_str() );

	if( !Finit )
	{
		InClassError( this, "LoadVQ(): Cannot open VQcodebook file " + dataFile + " .", -300 );
	}

	int SkipNumber = 6*Dimension*Dimension;

	// Ignore headers and mean data
	Finit.ignore( SkipNumber, '\n' );
	Finit.ignore( SkipNumber, '\n' );
	Finit.ignore( SkipNumber, '\n' );
	Finit.ignore( SkipNumber, '\n' );

	int i = 0;
	double value;

	while( i < Dimension )
	{
		Finit >> value;
		(*globalvars)[i++] = value;
	}

//	(*globalvars) = pow( (*globalvars), -1.0 );
//	(*globalvars) *= vFloor;

	int type = 0;
	Finit >> type;	// magic
	Finit >> type; // type

	if( type != 0 )
	{
		InClassError( this, "LoadVQ(): Cannot read in tree-structed VQ codebook.", -301 );
	}

	int covkind = 0;
	Finit >> covkind;

	if( covkind != 1 )
	{
		InClassError( this, "LoadVQ(): File contains non-diagonal covariance values.", -302 );
	}

	int clusterNumber = 0;
	Finit >> clusterNumber;

	if( clusterNumber != MixtureNumber )
	{
		InClassError( this, "LoadVQ(): The cluster number does not match mixture number.", -303 );
	}

	int streamNumber = 0;
	Finit >> streamNumber;

	if( streamNumber != 1 )
	{
		InClassError( this, "LoadVQ(): More than one stream in file.", -304 );
	}

	int streamWidth = 0;
	Finit >> streamWidth;

	if( streamWidth != Dimension )
	{
		InClassError( this, "LoadVQ(): Stream width not equal to feature vector dimension.", -305 );
	}

	// Readout cluster header
	Finit.ignore( 40, '\n' );
	Finit.ignore( 40, '\n' );

	// Read cluster data
	int j;
	i = 0;

	while( i < MixtureNumber )
	{
		j = 0;

		while( j < Dimension )
		{
			Finit >> value;
			(*(*means)[i])[j++] = value;
		}

		j = 0;

		while( j < Dimension )
		{
			Finit >> value;
			(*(*variances)[i])[j++] = value;
		}
		i++;

		// Readout next header
		Finit.ignore( 40, '\n' );
		Finit.ignore( 40, '\n' );
		Finit.ignore( 40, '\n' );
	}

	i = 0;

	while( i < MixtureNumber )
	{
		(*weights)[i++] = 1.0/(double)MixtureNumber;
	}

	Finit.close();
}

double GMM::LogL( string dataList )
{
	Flist.open( dataList.c_str() );

	if( !Flist )
	{
		InClassError( this, "LogL(): Cannot open data list file " + dataList,  -400 );
	}

	string dataFile;
	float htkData = 0.0f;
	unsigned int i, j;
	VectorProcessNumber = 0;
	VectorsIgnored = 0;
	LL = 0.0;

	Flist >> dataFile;
	cout << "LogL()" << endl;

	while( !Flist.eof() )
	{
		Fdata.clear();
		Fdata.open( dataFile.c_str(), ios_base::binary );

		if( !Fdata )
		{
			InClassError( this, "SetupData(): Cannot open data file " + dataFile,  -502 );
		}

		Fdata.read( reinterpret_cast<char *> ( &Htk ), sizeof( Htk ) );

		unsigned int tmpSamples = Htk.nSamples;

		if( Htk.nSamples > MaxDataNumber )
		{
			while( tmpSamples >= MaxDataNumber )
			{

				i = 0;

				while( i < MaxDataNumber )
				{
					j = 0;
					while( j < Dimension )
					{
						Fdata.read(  reinterpret_cast<char *> ( &htkData ), sizeof( float ));
						(*(*dataParm)[i])[j++] = (double)htkData;
					}
					i++;
				}

				Score( MaxDataNumber );

				tmpSamples -= MaxDataNumber;
			}
		}

// process rest of samples
		i = 0;

		while( i < tmpSamples )
		{
			j = 0;
			while( j < Dimension )
			{
				Fdata.read(  reinterpret_cast<char *> ( &htkData ), sizeof( float ));
				(*(*dataParm)[i])[j++] = (double)htkData;
			}
			i++;
		}

		Score( tmpSamples );

		Fdata.close();
		Flist >> dataFile;
	}

	cout << "VectorProcessNumber\t" << VectorProcessNumber << endl;
	cout << "VectorsIgnored\t\t" << VectorsIgnored << endl;

	Flist.close();
	Flist.clear();

	return LL/(double)VectorProcessNumber;
}

inline void GMM::Score( unsigned int VectorNumber )
{
	double value = 0.0, enorm = 0.0;
	int T = 0, i, j;
	bool overflow;

	while( T < VectorNumber )
	{
		i = 0;
		overflow = false;

		while( i < MixtureNumber )
		{
			*distance = *((*dataParm)[T]) - *((*means)[i]);
			*distance = pow( *distance, 2.0 );
			*distance /= *((*variances)[i]);
			value = (*distance).sum();
			value *= -0.5;

			if( value < -700.0 )
			{
				overflow = true;
				VectorsIgnored++;
				break;
			}

			enorm = 1.0;
			j = 0;

			while( j < Dimension )
			{
				enorm *= (*(*variances)[i])[j++];
			}

			enorm = sqrt( enorm );
			enorm *= pow( 2.0*M_PI, (double)Dimension / 2.0 );
			enorm = pow( enorm, -1.0 );
			(*PR)[i] = exp(value)*enorm*(*weights)[i];

			i++;
		}

		if( overflow != true )
		{
			VectorProcessNumber++;
			LL += log( (*PR).sum() );
		}

		T++;
	}
}

void GMM::printModel()
{
	int i = 0, j = 0;

	cout << "WEIGHTS" << endl;
	while( i < MixtureNumber )
	{
		cout << (*weights)[i++] << " ";
	}
	cout << endl;

	i = 0;
	cout << "MEANS" << endl;
	while( i < MixtureNumber )
	{
		j = 0;
		cout << i << ": ";
		while( j < Dimension )
		{
			cout << (*(*means)[i])[j++] << " ";
		}
		cout << endl;
		i++;
	}
	cout << endl;

	i = 0;
	cout << "VARS" << endl;
	while( i < MixtureNumber )
	{
		j = 0;
		cout << i << ": ";
		while( j < Dimension )
		{
			cout << (*(*variances)[i])[j++] << " ";
		}
		cout << endl;
		i++;
	}
	cout << endl;

	i = 0;
	cout << "GLOBAL VARS" << endl;
	while( i < Dimension )
	{
		cout << (*globalvars)[i++] << " ";
	}
	cout << endl;


}

GMM::~GMM()
{
	int i = 0;

	while( i < MixtureNumber )
	{
		delete (*means)[i];
		delete (*variances)[i++];
	}

	i = 0;

	while( i < MaxDataNumber )
	{
		delete (*dataParm)[i++];
	}

	delete weights;
	delete means;
	delete variances;
	delete globalvars;
	delete dataParm;
	delete distance;
	delete PR;
}
