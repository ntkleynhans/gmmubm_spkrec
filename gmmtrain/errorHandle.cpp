#include "speaker.h"

void InClassError( Speaker *model, string Message, int ErrorCode )
{
	model->~Speaker();
	cout << "Error Encountered!" << endl;
	cout << Message << endl;
	exit(ErrorCode);
}
