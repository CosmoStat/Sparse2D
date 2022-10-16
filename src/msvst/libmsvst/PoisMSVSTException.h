/*
 * Filename : PoisMSVSTException.h
 * the definitions of exceptions
 */
 
#ifndef _PoisMSVSTEXCEPTION_H
#define _PoisMSVSTEXCEPTION_H
#include <stdlib.h>
using namespace std;

class PoisMSVSTException
{
	protected:
		string reason;
	public:
		PoisMSVSTException(string reasonInfo=""):reason(reasonInfo) {}
		string getReason() { return reason; }
};

class DataException : public PoisMSVSTException
{
	public:
		DataException(string reasonInfo=""):PoisMSVSTException(reasonInfo) {}
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
