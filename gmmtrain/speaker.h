#ifndef SPEAKER_H
#define SPEAKER_H

/**
	@author neil taylor kleynhans <ntkleynhans@csir.co.za>
*/

#include <iostream>
#include <sstream>
#include <valarray>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <getopt.h>

using std::ios_base;
using std::cout;
using std::endl;
using std::vector;
using std::valarray;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;

//! HTK Binary File header.
//! This is defined in the HTKBook.

typedef struct {
	unsigned int	nSamples,	//!< number of samples in ﬁle (4-byte integer)
			sampPeriod;	//!< sample period in 100ns units (4-byte integer)
	unsigned short	sampSize,	//!< number of bytes per sample (2-byte integer)
			parmKind;	//!< a code indicating the sample kind (2-byte integer)
}HTKHeader;

//! Model header format

typedef struct {
	unsigned int	MixtureNumber,	//!< mixture number.
			Dimension;	//!< feature vector dimension
	double		vFloor;		//!< global variance flooring value.
}ModelHeader;

//! Speaker object.
//! Handles model initialization, training or adapting and saving.

class Speaker {
	public:
		//! Constructor.
		/*!	\param Model name.
			\param Initialization file name.
			\param Initialization file type.
			\param Model mixture number.
			\param Feature vector dimension.
			\param Variance minimizing factor.
		*/
		Speaker( string =0, string =0, unsigned int =0, unsigned int =0, unsigned int =0, double =0.0, unsigned int =0, string =0 );

		void saveModel();
		void modifyModel( string, int, unsigned int );
		double LogL( string );
		void printModel();

		~Speaker();

	private:

		void loadModel( string );	// HTK binary format file : type = 1
		void loadVQ( string );		// VQ Text file format : type = 2

		void ExpectStep( unsigned int );
		void Train();
		void Adapt( unsigned int );
		void Score( unsigned int );

		ofstream Fmodel;	//!< Model file stream handle.
		ifstream Finit;		//!< Initial model file stream handle.
		ifstream Flist;		//!< Data list file stream handle.
		ifstream Fdata;		//!< Data file stream handle.
		ofstream Fresult;

		vector< valarray<double> * > *means;		//!< Model means container.
		vector< valarray<double> * > *variances;	//!< Model variances container.
		vector< valarray<double> * > *CPmeans;		//!< Copy of model means container.
		vector< valarray<double> * > *CPvariances;	//!< Copy of odel variances container.
		vector< valarray<double> * > *dataParm;		//!< Data vector container
		vector< valarray<double> * > *EX;		//!< first moment
		vector< valarray<double> * > *EX2;		//!< second moment

		valarray<double> *weights;	//!< Model weights container.
		valarray<double> *CPweights;	//!< Copy of model weights container.
		valarray<double> *globalvars;	//!< Global variances container
		valarray<double> *distance;
		valarray<double> *PR;
		valarray<double> *N;
		valarray<double> *DDA;

		unsigned int MixtureNumber;	//!< The mixture number of the model.
		unsigned int Dimension;		//!< The dimension of the feature vector.
		unsigned int VectorProcessNumber;	//!< The number of feature vectors loaded.
		unsigned int MaxDataNumber;	//!< Only load this amount of vectors at a time.
		unsigned int VectorsIgnored;	//!< The number of feature vectors during training/ adapting.
		unsigned int SpeakerIgnored;

		string ModelName;		//!< Store model name.
		double vFloor;		//!< Value multiplied by global variance values to give minimum variances values.
		double LL;

		HTKHeader Htk;
		ModelHeader Hmodel;
};

//! This function is called if a error occurs within the speaker class.
/*!	\param speaker model reference.
	\param error message string.
	\param exit code.
*/

void InClassError( Speaker *, string, int );

#endif
