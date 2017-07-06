#include "gmm.h"

void InClassError( GMM *model, string Message, int ErrorCode )
{
	model->~GMM();
	cout << "Error Encountered!" << endl;
	cout << Message << endl;
	exit(ErrorCode);
}
