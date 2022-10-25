/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  5/03/2000 
**    
**    File:  PrimeNumber.h
**
**    Modification history :
**
**
******************************************************************************/

#ifndef _PrimeNumber_H_
#define _PrimeNumber_H_

#include "GlobalInc.h"

#define NBR_PRIME_NUMBER 10000
#define LAST_PRIME_NUMBER 104729
#define T_UL unsigned long

extern T_UL TabPrimeNumber10000[NBR_PRIME_NUMBER];

/***********************************************************************/

class CPRIME_NUMBER {
        T_UL  NbrPrime;
        T_UL  LastPrime;
        T_UL  *TabPN;
        void test_n(T_UL  N)
        {
           if (N < 0)
           {
               cout << "Error in CPRIME_NUMBER: N is < 0 : " << N << endl;
               exit(-1);
           }
           if (N > LastPrime)
           {
               cout << "Error in CPRIME_NUMBER: N is too high : " << N << endl;
               exit(-1);
           }
        }
        public:
           CPRIME_NUMBER() {NbrPrime=NBR_PRIME_NUMBER;
                            TabPN= TabPrimeNumber10000;
                            LastPrime=LAST_PRIME_NUMBER;}
           ~CPRIME_NUMBER() {};
	   Bool is_prime_number(T_UL  N) // is prime number ?
           {
               Bool ValRet = False;
               int i = 0;
               test_n(N);
               while ((i < NBR_PRIME_NUMBER) && (TabPN[i] < N)) i++;
               if ((i < NBR_PRIME_NUMBER) && (TabPN[i] == N)) ValRet = True;
               return ValRet;
           }
           T_UL kth_prime_number(int Num) // kth prime number
           {
               T_UL ValRet = 0;
               if ((Num >= 0) && (Num < NBR_PRIME_NUMBER))
                    ValRet = TabPN[Num];
               return ValRet;
           }
           T_UL next_prime_number(T_UL  N) // higher or equal prime number
           {
               T_UL ValRet = 0;
               int i = 0;
               test_n(N);
               while ((i < NBR_PRIME_NUMBER) && (TabPN[i] < N)) i++;
               if (i < NBR_PRIME_NUMBER) ValRet = TabPN[i];
               return ValRet;
           }
           T_UL previous_prime_number(T_UL  N) // lower prime number
           {
               T_UL ValRet = 0;
               int i = 0;
               test_n(N);
               while ((i < NBR_PRIME_NUMBER) && (TabPN[i] < N)) i++;
               if ((i > 0) && (i < NBR_PRIME_NUMBER)) ValRet = TabPN[i-1];
               return ValRet;
           }
};


/***********************************************************************/

#endif
