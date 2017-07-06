#include "speaker.h"

Speaker::Speaker( string modelName, string modelInitFile, unsigned int initType, unsigned int mixtures, unsigned int length, double floor, unsigned int dataSize, string resFile )
{
	MixtureNumber = mixtures;
	Dimension = length;
	ModelName = modelName;
	vFloor = floor;
	MaxDataNumber = dataSize;

	means = new vector< valarray<double> * > (MixtureNumber);
	variances = new vector< valarray<double> * > (MixtureNumber);
	CPmeans = new vector< valarray<double> * > (MixtureNumber);
	CPvariances = new vector< valarray<double> * > (MixtureNumber);
	dataParm = new vector< valarray<double> * >(MaxDataNumber);
	EX = new vector< valarray<double> * >( MixtureNumber );
	EX2 = new vector< valarray<double> * >( MixtureNumber );
	
	weights = new valarray<double>( 0.0, MixtureNumber );
	CPweights = new valarray<double>( 0.0, MixtureNumber );
	globalvars = new valarray<double>( 0.0, Dimension );
	distance = new valarray<double>( 0.0, Dimension );
	PR = new valarray<double>( 0.0, MixtureNumber );
	N = new valarray<double>( 0.0, MixtureNumber );
	DDA = new valarray<double>( 0.0, MixtureNumber );

	unsigned int i = 0;

	while( i < MixtureNumber )
	{
		(*means)[i] = new valarray<double>( 0.0, Dimension );
		(*variances)[i] = new valarray<double>( 0.0, Dimension );
		(*CPmeans)[i] = new valarray<double>( 0.0, Dimension );
		(*CPvariances)[i] = new valarray<double>( 0.0, Dimension );
		(*EX)[i] = new valarray<double>( 0.0, Dimension );
		(*EX2)[i++] = new valarray<double>( 0.0, Dimension );
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
		InClassError( this, "Speaker(): InitType error value specified not known.", -100 );
	}

	Fresult.open( resFile.c_str() );

	if( !Fresult )
	{
		InClassError( this, "LoadModel(): Cannot open results file " + resFile + " .", -101 );
	}
}

void Speaker::loadModel( string modelFile )
{
	Finit.open( modelFile.c_str(), ios_base::binary );

	Finit.read( reinterpret_cast <char *>( &Hmodel ), sizeof( ModelHeader ) );

	if( !Finit )
	{
		InClassError( this, "LoadModel(): Cannot open model file " + modelFile + " .", -300 );
	}

	if( Hmodel.MixtureNumber != MixtureNumber )
	{
		InClassError( this, "LoadModel(): Model file reports non-equal mixture number", -301 );
	}

	if( Hmodel.Dimension != Dimension )
	{
		InClassError( this, "LoadModel(): Model file reports non-equal feature vector dimension", -302 );
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

void Speaker::loadVQ( string dataFile )
{
	Finit.open( dataFile.c_str() );

	if( !Finit )
	{
		InClassError( this, "LoadVQ(): Cannot open VQcodebook file " + dataFile + " .", -200 );
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
	(*globalvars) *= vFloor;

	int type = 0;
	Finit >> type;	// magic
	Finit >> type; // type

	if( type != 0 )
	{
		InClassError( this, "LoadVQ(): Cannot read in tree-structed VQ codebook.", -201 );
	}

	int covkind = 0;
	Finit >> covkind;

	if( covkind != 1 )
	{
		InClassError( this, "LoadVQ(): File contains non-diagonal covariance values.", -202 );
	}

	int clusterNumber = 0;
	Finit >> clusterNumber;

	if( clusterNumber != MixtureNumber )
	{
		InClassError( this, "LoadVQ(): The cluster number does not match mixture number.", -203 );
	}

	int streamNumber = 0;
	Finit >> streamNumber;

	if( streamNumber != 1 )
	{
		InClassError( this, "LoadVQ(): More than one stream in file.", -204 );
	}

	int streamWidth = 0;
	Finit >> streamWidth;

	if( streamWidth != Dimension )
	{
		InClassError( this, "LoadVQ(): Stream width not equal to feature vector dimension.", -205 );
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
			//value = pow( value, -1.0 );

			if( value < (*globalvars)[j] )
			{
				(*(*variances)[i])[j] = (*globalvars)[j];
			}
			else 
			{
				(*(*variances)[i])[j] = value;
			}
			j++;
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

void Speaker::saveModel( void )
{
	Fmodel.open( ModelName.c_str(), ios_base::binary );

	if( !Fmodel )
	{
		InClassError( this, "SaveModel(): Cannot open output model file " + ModelName,  -400 );
	}

	Hmodel.MixtureNumber = MixtureNumber;
	Hmodel.Dimension = Dimension;
	Hmodel.vFloor = vFloor;

	Fmodel.write( reinterpret_cast< char * > ( &Hmodel ), sizeof( ModelHeader ) );

	int i = 0;

	while( i < Dimension )
	{
		Fmodel.write( reinterpret_cast <char *>( &(*globalvars)[i++] ), sizeof(double) );
	}

	int j = 0;
	i = 0;

	while( i < MixtureNumber )
	{
		Fmodel.write( reinterpret_cast <char *>( &(*weights)[i] ), sizeof(double) );

		while( j < Dimension )
		{
			Fmodel.write( reinterpret_cast <char *>( &(*(*means)[i])[j] ), sizeof(double) );
			Fmodel.write( reinterpret_cast <char *>( &(*(*variances)[i])[j++] ), sizeof(double) );
		}

		j = 0;
		i++;
	}
	
	Fmodel.close();
}

void Speaker::modifyModel( string dataList, int task, unsigned int flags )
{
	Flist.open( dataList.c_str() );

	if( !Flist )
	{
		saveModel();
		InClassError( this, "SetupData(): Cannot open data list file " + dataList,  -500 );
	}

	if( task != 1 && task != 2 )
	{
		saveModel();
		InClassError( this, "SetupData(): Unknown task given",  -501 );
	}

	string dataFile;
	float htkData = 0.0f;
	unsigned int i, j;
	VectorProcessNumber = 0;
	VectorsIgnored = 0;

	Flist >> dataFile;

	cout << "ModifyModel()" << endl;
	while( !Flist.eof() )
	{
		Fdata.clear();
		Fdata.open( dataFile.c_str(), ios_base::binary );

		if( !Fdata )
		{
			InClassError( this, "SetupData(): Cannot open data file " + dataFile,  -502 );
		}

		Fdata.read( reinterpret_cast<char *> ( &Htk ), sizeof( Htk ) );

		Fresult << dataFile << endl;
		unsigned int tmpSamples = Htk.nSamples;
		SpeakerIgnored = 0;

		if( Htk.nSamples > MaxDataNumber )
		{
			while( tmpSamples > MaxDataNumber )
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

				ExpectStep( MaxDataNumber );

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

		ExpectStep( tmpSamples );

		Fdata.close();
		Flist >> dataFile;
		Fresult << SpeakerIgnored << " " <<  Htk.nSamples << endl;
	}

	cout << "VectorProcessNumber\t" << VectorProcessNumber << endl;
	cout << "VectorsIgnored\t\t" << VectorsIgnored << endl;

	Fresult << "Train" << endl;
	Fresult << "VectorProcessNumber\t" << VectorProcessNumber << endl;
	Fresult << "VectorsIgnored\t\t" << VectorsIgnored << endl;
	Fresult << "Percent \t\t" << (float) ((float)VectorsIgnored)/ ((float)VectorProcessNumber) * 100.0f << endl;
	Fresult << endl;

	i = 0;

	while( i < MixtureNumber )
	{
		(*(*EX)[i]) /= (*N)[i];
		(*(*EX2)[i]) /= (*N)[i];
		i++;
	}

	if( task == 1 )
	{
		Train();
	}
	else if( task == 2 )
	{
		Adapt( flags );
	}

	Flist.close();
	Flist.clear();

	(*N) = 0.0;
	i = 0;
	while( i < MixtureNumber )
	{
		(*(*EX)[i]) = 0.0;
		(*(*EX2)[i]) = 0.0;
		i++;
	}
}

inline void Speaker::ExpectStep( unsigned int VectorNumber )
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
				cout << "Warning: value to small " << value << endl;
				SpeakerIgnored++;
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
			(*PR) /= (*PR).sum();
			(*N) += (*PR);

			i = 0;

			while( i < MixtureNumber )
			{
				(*(*EX)[i]) += (*PR)[i] * (*(*dataParm)[T]);
				(*(*EX2)[i]) += (*PR)[i] * pow( (*(*dataParm)[T]), 2.0 );
				i++;
			}
		}

		T++;
	}
}

void Speaker::Train()
{
	int i = 0, j;

	(*N) /= VectorProcessNumber;

	while( i < MixtureNumber )
	{
		//(*(*EX2)[i]) -= pow( (*(*means)[i]), 2.0 );
		(*(*EX2)[i]) -= pow( (*(*EX)[i]), 2.0 );
		i++;
	}

	i = 0;

	while( i < MixtureNumber )
	{
		j = 0;
		while( j < Dimension )
		{
			if( (*(*EX2)[i])[j] < (*globalvars)[j] )
			{
				//cout << (*(*EX2)[i])[j] << " " << i << " " << j << endl;
				(*(*EX2)[i])[j] = (*globalvars)[j];
			}
			j++;
		}
		i++;
	}

	(*weights) = (*N);
	i = 0;

	while( i < MixtureNumber )
	{
		(*(*means)[i]) = (*(*EX)[i]);
		(*(*variances)[i]) = (*(*EX2)[i]);
		i++;
	}
}

void Speaker::Adapt( unsigned int flag )
{
	int i;
	(*DDA) = (*N) / ( (*N) + 16.0 );
	(*N) /= VectorProcessNumber;

	if( flag & 2 || flag & 4 )
	{
		i = 0;
		while( i < MixtureNumber )
		{
			(*(*CPmeans)[i]) = (*(*means)[i]);
			(*(*CPvariances)[i]) = (*(*variances)[i]);
			i++;
		}
	}

	if( flag & 1 )
	{
		(*CPweights) = (*weights);
		(*weights) = (*DDA)*(*N) + ( 1.0 - (*DDA) )*(*CPweights);
		(*weights) /= (*weights).sum();
	}

	if( flag & 2 )
	{
		i = 0;
		while( i < MixtureNumber )
		{
			(*(*means)[i]) = (*DDA)[i]*(*(*EX)[i]) + ( 1.0 - (*DDA)[i] ) * (*(*CPmeans)[i]);
			i++;
		}
	}

	if( flag & 4 )
	{
		i = 0;
		while( i < MixtureNumber )
		{
			(*(*variances)[i]) = (*DDA)[i]*(*(*EX2)[i]) + ( 1.0 - (*DDA)[i] ) * ( pow( (*(*CPmeans)[i]), 2.0 ) + (*(*CPvariances)[i]) ) - pow( (*(*means)[i]), 2.0 );
			i++;
		}

		int j;
		i = 0;
		while( i < MixtureNumber )
		{
			j = 0;
			while( j < Dimension )
			{
				if( (*(*variances)[i])[j] < (*globalvars)[j] )
				{
					//cout << (*(*EX2)[i])[j] << " " << i << " " << j << endl;
					(*(*variances)[i])[j] = (*globalvars)[j];
				}
				j++;
			}
			i++;
		}
	}
}

double Speaker::LogL( string dataList )
{
	Flist.open( dataList.c_str() );

	if( !Flist )
	{
		InClassError( this, "LogL(): Cannot open data list file " + dataList,  -600 );
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
			InClassError( this, "LogL(): Cannot open data file " + dataFile,  -601 );
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

	Fresult << "LL\t\t" << LL/(double)VectorProcessNumber << endl;
	Fresult << endl;

	Flist.close();
	Flist.clear();

	return LL/(double)VectorProcessNumber;
}

inline void Speaker::Score( unsigned int VectorNumber )
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

void Speaker::printModel()
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

Speaker::~Speaker()
{
	int i = 0;

	Fresult.close();

	while( i < MixtureNumber )
	{
		delete (*means)[i];
		delete (*variances)[i];
		delete (*CPmeans)[i];
		delete (*CPvariances)[i];
		delete (*EX)[i];
		delete (*EX2)[i++];
	}

	i = 0;

	while( i < MaxDataNumber )
	{
		delete (*dataParm)[i++];
	}

	delete weights;
	delete CPweights;
	delete means;
	delete variances;
	delete CPmeans;
	delete CPvariances;
	delete globalvars;
	delete dataParm;
	delete distance;
	delete PR;
	delete N;
	delete EX;
	delete EX2;
	delete DDA;
}
