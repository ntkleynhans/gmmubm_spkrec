#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct {
	unsigned int	nSamples,	//!< number of samples in ﬁle (4-byte integer)
			sampPeriod;	//!< sample period in 100ns units (4-byte integer)
	unsigned short	sampSize,	//!< number of bytes per sample (2-byte integer)
			parmKind;	//!< a code indicating the sample kind (2-byte integer)
}HTKHeader;

static float **data;
static double **old_mean, **new_mean, *sum_num, *score, *global_mean;
static double **var, *global_var, old_error, new_error;
static int dims, cluster_size, vector_num;
static char *data_list, *out_file, *res_file;
static unsigned int total;

void init( void );
void assign_mean( int );
void assign_var( int );
int findmin( double * );
void cluster( void );
void alloc_mem( void );
void free_mem( void );
void get_data( FILE *, int );
void print( int );
void output_cluster( void );
void calculate_var( void );
void load_parms( char *);

int main( int argc, char *argv[] )
{
	double diff;
	FILE *fres;
	int i, j;

	if( argc < 2 )
	{
		printf( "Input parameter file missing\n" );
		printf("Kmeans (parameter file)\n");
		exit(-1);
	}

	load_parms( argv[1] );

	old_error = 0.0;
	alloc_mem();
	init();

	do {
		cluster();

		diff = fabs( old_error-new_error );
		diff /= new_error;
		old_error = new_error;
		printf( "Error %f\n", old_error );

	} while( diff > 0.001 );

	output_cluster();

	fres = fopen( res_file, "w" );

	if( fres == NULL )
	{
		printf( "Main(): Cannot open results file %s\n", res_file );
		exit(-1);
	}

	fprintf( fres, "Cluster\tTotal\n" );
	j = 0;
	for( i = 0; i < cluster_size; i++ )
	{
		fprintf( fres, "%d\t%.2f\n", i+1, sum_num[i] );
		if( sum_num[i] == 0.0 )
		{
			j++;
		}
	}

	fclose(fres);

	if( j != 0 )
	{
		printf( "Warning: Clusters with zero elements (see RESULT file): %.1f %%\n", 100.0f*(float)j/(float)cluster_size );
	}

	free_mem();

	return 0;
}

void load_parms( char *file_name )
{
	FILE *fin;
	char string[256];
	int i;

	fin = fopen( file_name, "r" );

	if( fin == NULL )
	{
		printf( "Load_Parms(): Cannot open file %s\n", file_name );
		exit( -1 );
	}

	data_list = malloc( sizeof(char) * 256 );
	out_file = malloc( sizeof(char) * 256 );
	res_file = malloc( sizeof(char) * 256 );
	i = 0;

	fscanf( fin, "%s", string );
	while( !feof( fin ) )
	{
		if( strcmp( "VECTOR", string ) == 0 )
		{
			i++;
			fscanf( fin, "%s", string );
			vector_num = atoi( string );
		}
		else if( strcmp( "DIMS", string ) == 0 )
		{
			i++;
			fscanf( fin, "%s", string );
			dims = atoi( string );
		}
		else if( strcmp( "CLUSTER", string ) == 0 )
		{
			i++;
			fscanf( fin, "%s", string );
			cluster_size = atoi( string );
		}
		else if( strcmp( "LIST", string ) == 0 )
		{
			i++;
			fscanf( fin, "%s", string );
			strcpy( data_list, string );
		}
		else if( strcmp( "VQOUT", string ) == 0 )
		{
			i++;
			fscanf( fin, "%s", string );
			strcpy( out_file, string );
		}
		else if( strcmp( "RESULT", string ) == 0 )
		{
			i++;
			fscanf( fin, "%s", string );
			strcpy( res_file, string );
		}

		fscanf( fin, "%s", string );
	}

	if( i != 6 )
	{
		printf( "Load_Parms(): Some parameters not set - Check script file %s\n", file_name );
		exit( -1 );
	}

	fclose( fin );
}

void alloc_mem( void )
{
	int i;

	data = malloc( sizeof( double * )*vector_num );

	i = 0;
	while( i < vector_num )
		data[i++] = malloc( sizeof( double )*dims );

	old_mean = malloc( sizeof( double * )*cluster_size );
	new_mean = malloc( sizeof( double * )*cluster_size );
	var = malloc( sizeof( double * )*cluster_size );

	i = 0;
	while( i < cluster_size )
	{
		old_mean[i] = malloc( sizeof( double )*dims );
		new_mean[i] = malloc( sizeof( double )*dims );
		var[i] = malloc( sizeof( double )*dims );
		i++;
	}

	global_mean = malloc( sizeof( double )*dims );
	global_var = malloc( sizeof( double )*dims );
	score = malloc( sizeof( double )*cluster_size );
	sum_num = malloc( sizeof( double )*cluster_size );
}

void free_mem( void )
{
	int i;

	i = 0;
	while( i < vector_num )
		free(data[i++]);

	i = 0;
	while( i < cluster_size )
	{
		free( old_mean[i] );
		free( new_mean[i] );
		free( var[i] );
		i++;
	}

	free( data );
	free( old_mean );
	free( new_mean );
	free( var );
	free( global_mean );
	free( global_var );
	free( score );
	free( sum_num );
	free( data_list );
	free( out_file );
	free( res_file );
}

void init( void )
{
	FILE *flist, *fdata;
	int i, j, k, sum;
	char string[256];
	HTKHeader htk;
	double dev;

//clean out global variables
	for( k = 0; k < dims; k++ )
	{
		global_mean[k] = 0.0;
		global_var[k] = 0.0;
	}

	flist = fopen( data_list, "r" );
	if( flist == NULL )
	{
		printf( "init(): Cannot open datalist file %s\n", string );
		exit(-1);
	}

	fscanf( flist, "%s", string );

//calculate global mean
	sum = 0;
	while( !feof( flist ) )
	{
		fdata = fopen( string, "rb" );
		if( fdata == NULL )
		{
			printf( "init(): Cannot open data file %s\n", string );
			exit(-1);
		}

		fread( &htk, sizeof( HTKHeader ), 1, fdata );

		i = htk.nSamples;
		if( htk.nSamples > vector_num )
		{
			while( i >= vector_num )
			{
				get_data( fdata, vector_num );
	
				for( j = 0; j < vector_num; j++ )
				{
					for( k = 0; k < dims; k++ )
					{
						global_mean[k] += data[j][k];
					}
				}
				sum += vector_num;
				i -= vector_num;
			}
		}

		get_data( fdata, i );

		for( j = 0; j < i; j++ )
		{
			for( k = 0; k < dims; k++ )
				global_mean[k] += data[j][k];
		}
		sum += i;

		fclose( fdata );
		fscanf( flist, "%s", string );
	}

	for( k = 0; k < dims; k++ )
		global_mean[k] /= sum;

//calculate global covariances
	rewind( flist );
	fscanf( flist, "%s", string );
	sum = 0;
	while( !feof( flist ) )
	{
		fdata = fopen( string, "rb" );
		if( fdata == NULL )
		{
			printf( "init(): Cannot open data file %s\n", string );
			exit(-1);
		}
	
		fread( &htk, sizeof( HTKHeader ), 1, fdata );

		i = htk.nSamples;
		if( htk.nSamples > vector_num )
		{
			while( i >= vector_num )
			{
				get_data( fdata, vector_num );
	
				for( j = 0; j < vector_num; j++ )
				{
					for( k = 0; k < dims; k++ )
						global_var[k] += pow( global_mean[k]-data[j][k], 2.0 );
				}
				sum += vector_num;
				i -= vector_num;
			}
		}

		get_data( fdata, i );

		for( j = 0; j < i; j++ )
		{
			for( k = 0; k < dims; k++ )
				global_var[k] += pow( global_mean[k]-data[j][k], 2.0 );
		}
		sum += i;

		fclose( fdata );
		fscanf( flist, "%s", string );
	}

	for( k = 0; k < dims; k++ )
		global_var[k] /= sum;

	fclose( flist );

// calculate mean
	dev = 2.0/(double)cluster_size;

	for( i = 0; i < cluster_size; i++ )
	{
		for( j = 0; j < dims; j++ )
		{
			
			old_mean[i][j] = global_mean[j] + dev*sqrt( global_var[j] )*( (double)rand() / (double)RAND_MAX + 1.0 );
		}
	}
}

void cluster( void )
{
	FILE *flist, *fdata;
	int i, j;
	char string[256];
	HTKHeader htk;

	new_error = 0.0;
	total = 0;

	for( i = 0; i < cluster_size; i++ )
	{
		for( j = 0; j < dims; j++ )
		{
			new_mean[i][j] = 0.0;
		}
		sum_num[i] = 0.0;
	}

	flist = fopen( data_list, "r" );
	if( flist == NULL )
	{
		printf( "cluster(): Cannot open datalist file %s\n", string );
		exit(-1);
	}

	fscanf( flist, "%s", string );

//calculate mean
	while( !feof( flist ) )
	{
		fdata = fopen( string, "rb" );
		if( fdata == NULL )
		{
			printf( "cluster(): Cannot open data file %s\n", string );
			exit(-1);
		}
	
		fread( &htk, sizeof( HTKHeader ), 1, fdata );

		i = htk.nSamples;
		if( htk.nSamples > vector_num )
		{
			while( i >= vector_num )
			{
				get_data( fdata, vector_num );
				assign_mean( vector_num );	
				i -= vector_num;
				total += vector_num;
			}
		}

		get_data( fdata, i );
		assign_mean( i );
		total += i;

		fclose( fdata );
		fscanf( flist, "%s", string );
	}

	fclose( flist );

	for( i = 0; i < cluster_size; i++ )
	{
		if( sum_num[i] != 0.0 )
		{
			for( j = 0; j < dims; j++ )
			{
				old_mean[i][j] = new_mean[i][j] / sum_num[i];
			}
		}
	}
	new_error /= (double)total;
}

void assign_mean( int size )
{
	int i, j, k;

	for( i = 0; i < size; i++ )
	{
		for( j = 0; j < cluster_size; j++ )
		{
			score[j] = 0.0;
			for( k = 0; k < dims; k++ )
			{
				score[j] += pow( ( data[i][k] - old_mean[j][k] ), 2.0 );
			}
		}

		j = findmin( score );
		new_error += score[j]/(double)dims;

		for( k = 0; k < dims; k++ )
		{
			new_mean[j][k] += data[i][k];
		}
		sum_num[j]++;
	}
}

void assign_var( int size )
{
	int i, j, k;

	for( i = 0; i < size; i++ )
	{
		for( j = 0; j < cluster_size; j++ )
		{
			score[j] = 0.0;
			for( k = 0; k < dims; k++ )
			{
				score[j] += pow( ( data[i][k] - old_mean[j][k] ), 2.0 );
			}
		}

		j = findmin( score );

		for( k = 0; k < dims; k++ )
		{
			var[j][k] += pow( old_mean[j][k] - data[i][k], 2.0 );
		}
		sum_num[j]++;
	}
}

void get_data( FILE *fin, int number )
{
	int i;

	for( i = 0; i < number; i++ )
	{
		fread( data[i], sizeof( float )*dims, 1, fin );
	}
}

int findmin( double *data )
{
	int i, j;
	double var = data[0];
	j = 0;

	for( i = 1; i < cluster_size; i++ )
	{
		if( data[i] < var )
		{
			var = data[i];
			j = i;
		}
	}

	return j;
}

void print( int flag )
{
	int i, j;

	if( flag & 1 ) // means
	{
		for( i = 0; i < cluster_size; i++ )
		{
			printf( "Cluster %d\n", i );
			for( j = 0; j < dims; j++ )
			{
				printf( "%f ", old_mean[i][j] );
			}
			printf( "\n" );
		}
	}

	if( flag & 2 ) // vars
	{
		for( i = 0; i < cluster_size; i++ )
		{
			printf( "Cluster %d\n", i );
			for( j = 0; j < dims; j++ )
			{
				printf( "%f ", var[i][j] );
			}
			printf( "\n" );
		}
	}

	if( flag & 4 ) // global_mean
	{
		for( j = 0; j < dims; j++ )
		{
			printf( "%f ", global_mean[j] );
		}
		printf( "\n" );
	}

	if( flag & 8 ) // global_var
	{
		for( j = 0; j < dims; j++ )
		{
			printf( "%f ", global_var[j] );
		}
		printf( "\n" );
	}

	if( flag & 16 ) // vars
	{
		for( i = 0; i < vector_num; i++ )
		{
			printf( "Data node %d\n", i );
			for( j = 0; j < dims; j++ )
			{
				printf( "%f ", data[i][j] );
			}
			printf( "\n" );
		}
	}
}

void output_cluster( void )
{
	FILE *fout;
	int i, j;

	fout = fopen( out_file, "w" );

	if( fout == NULL )
	{
		printf( "Output_cluster(): Cannot open file %s\n", out_file );
		exit(-1);
	}

	calculate_var();

	fprintf( fout, "Global mean:\n\n" );

	for( j = 0; j < dims; j++ )
	{
		fprintf( fout, "%f\t", global_mean[j] );
	}

	fprintf( fout, "\nGlobal variance:\n" );

	for( j = 0; j < dims; j++ )
	{
		fprintf( fout, "%f\t", global_var[j] );
	}

	fprintf( fout, "\n262 0 1 %d 1 %d\n", cluster_size, dims );

	for( i = 0; i < cluster_size; i++ )
	{
		fprintf( fout, "1 %d %d 0 %d\n", cluster_size - i, i+1, i+2 );

		for( j = 0; j < dims; j++ )
		{
			fprintf( fout, "%f ", old_mean[i][j] );
		}
		fprintf( fout, "\n" );

		for( j = 0; j < dims; j++ )
		{
			fprintf( fout, "%f ", var[i][j] );
		}
		fprintf( fout, "\n\n" );
	}

	fclose( fout );
}

void calculate_var( void )
{
	FILE *flist, *fdata;
	int i, j;
	char string[256];
	HTKHeader htk;

	for( i = 0; i < cluster_size; i++ )
	{
		for( j = 0; j < dims; j++ )
		{
			var[i][j] = global_var[j];
		}
		sum_num[i] = 0.0;
	}

	flist = fopen( data_list, "r" );
	if( flist == NULL )
	{
		printf( "calculate_var(): Cannot open datalist file %s\n", string );
		exit(-1);
	}

	fscanf( flist, "%s", string );

//calculate var
	while( !feof( flist ) )
	{
		fdata = fopen( string, "rb" );
		if( fdata == NULL )
		{
			printf( "calculate_var(): Cannot open data file %s\n", string );
			exit(-1);
		}
	
		fread( &htk, sizeof( HTKHeader ), 1, fdata );

		i = htk.nSamples;
		if( htk.nSamples > vector_num )
		{
			while( i >= vector_num )
			{
				get_data( fdata, vector_num );
				assign_var( vector_num );	
				i -= vector_num;
			}
		}

		get_data( fdata, i );
		assign_var( i );

		fclose( fdata );
		fscanf( flist, "%s", string );
	}

	fclose( flist );

	for( i = 0; i < cluster_size; i++ )
	{
		if( sum_num[i] != 0.0 )
		{
			for( j = 0; j < dims; j++ )
			{
				var[i][j] /=  sum_num[i];
			}
		}
	}
}
