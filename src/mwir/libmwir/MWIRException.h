/*
 * Filename : MWIRException.h
 * the definitions of exceptions
 */
 
#ifndef _MWIREXCEPTION_H
#define _MWIREXCEPTION_H
#include <stdlib.h>
using namespace std;

class MWIRException
{
	protected:
		string reason;
	public:
		MWIRException(string reasonInfo=""):reason(reasonInfo) {}
		string getReason() { return reason; }
};

class DataException : public MWIRException
{
	public:
		DataException(string reasonInfo=""):MWIRException(reasonInfo) {}
};

template <class T>
class DataIDException : public DataException
{
	protected:
		T ID;
	public:
		DataIDException(T dataID, string reasonInfo=""):DataException(reasonInfo),ID(dataID) {}
		T getID() { return ID; }
};

class DataSizeException : public DataException
{
	protected:
		int size;
	public:
		DataSizeException (int dataSize, string reasonInfo=""):DataException(reasonInfo) { size = dataSize;	}
		int getSize() { return size; }
};
#endif
