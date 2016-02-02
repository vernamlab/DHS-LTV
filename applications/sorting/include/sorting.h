/*
 * sorting.h
 *
 *  Created on: Apr 29, 2014
 *      Author: dante
 */

#ifndef SORTING_H_
#define SORTING_H_

#include "instruct.h"
#include <iostream>
#include <algorithm>

//////////////////////////////////////////////////////////////////
// Functions For Sorting Tests
void Print(string s, int *l, int size){

	cout << s << endl;
	for(int i=0; i<size; i++)
		cout << l[i] << "\t";
	cout << endl;
}

void MySort(int *out, int *in, int size){

	for(int i=0; i<size; i++)
		out[i] = in[i];

	sort(out, out + size);
}
////////////////////////////////////////////////////////////////////

void SortingList_SetEncrypt(Instruct &t, InstByte *x, int *list, int size){

	for(int i=0; i<size; i++)
		x[i] = t.InttoInstByte(list[i]);

	for(int i=0; i<size; i++)
		x[i] = t.Encrypt(x[i]);
}

void SortingList_SetDecrypt(Instruct &t, int *x, InstByte *y, int size){
	InstByte e;
	for(int i=0; i<size; i++){
		e = t.Decrypt(y[i]);
		x[i] = t.InstBytetoInt(e);
	}
}

void SortingList_SetEncrypt(Instruct &t, InstInt *x, int *list, int size){

	for(int i=0; i<size; i++)
		x[i] = t.InttoInstInt(list[i]);

	for(int i=0; i<size; i++)
		x[i] = t.Encrypt(x[i]);
}

void SortingList_SetDecrypt(Instruct &t, int *x, InstInt *y, int size){
	InstInt e;
	for(int i=0; i<size; i++){
		e = t.Decrypt(y[i]);
		x[i] = t.InstInttoInt(e);
	}
}

void TestSort(){

    SortingAlgorithm algo = (SortingAlgorithm)Sort_Alg;
    int elt_size = Sort_Size;
    int bit_size = Bit_Size;

    myTimer mt;

	Instruct t(Modulus_M, Num_Primes, Max_Prime);
	t.GenerateParameters();
	t.GenerateKeys();

    // Create Random List
	int *randList = new int[elt_size];
	srand (time(NULL));
	for(int j=0; j<elt_size; j++)
	{
        if(bit_size == 8)
            randList[j] = rand()%256;
		else
            randList[j] = rand();
    }

    // Sort the list for testing
	int *sortedList = new int [elt_size];
	MySort(sortedList, randList, elt_size);
	Print("X:", randList, elt_size);

    // Convert the input and output to FHE data types
    if(bit_size == 8)
    {
        InstByte *x = new InstByte[elt_size];
        InstByte *y = new InstByte[elt_size];
        int	*d = new int[elt_size];

        // Encrypt the input list
        SortingList_SetEncrypt(t, x, randList, elt_size);

        if(algo == 1)
        {
            mt.Start();
            t.BubbleSort(y, x, elt_size);
            mt.Stop();
        }
        else if(algo == 2)
        {
            mt.Start();
            t.BatcherSort(y, x, elt_size);
            mt.Stop();
        }
        else if(algo == 3)
        {
            mt.Start();
            t.DirectSort(y, x, elt_size);
            mt.Stop();
        }
        else if(algo == 4)
        {
            mt.Start();
            t.GreedySort(y, x, elt_size);
            mt.Stop();
        }

        mt.ShowTime("Sort Time:\t");

        // Decrypt the output list
        SortingList_SetDecrypt(t, d, y, elt_size);
        Print("Y:", d, elt_size);
        Print("Check:", sortedList, elt_size);

        delete [] x;
        delete [] y;
        delete [] d;
    }
    else if(bit_size == 32)
    {
        InstInt *x = new InstInt[elt_size];
        InstInt *y = new InstInt[elt_size];
        int	*d = new int[elt_size];

        // Encrypt the input list
        SortingList_SetEncrypt(t, x, randList, elt_size);

        if(algo == 1)
        {
            mt.Start();
            t.BubbleSort(y, x, elt_size);
            mt.Stop();
        }
        else if(algo == 2)
        {
            mt.Start();
            t.BatcherSort(y, x, elt_size);
            mt.Stop();
        }
        else if(algo == 3)
        {
            mt.Start();
            t.DirectSort(y, x, elt_size);
            mt.Stop();
        }
        else if(algo == 4)
        {
            mt.Start();
            t.GreedySort(y, x, elt_size);
            mt.Stop();
        }
        mt.ShowTime("Sort Time:\t");

        // Decrypt the output list
        SortingList_SetDecrypt(t, d, y, elt_size);
        Print("Y:", d, elt_size);
        Print("Check:", sortedList, elt_size);

        delete [] x;
        delete [] y;
        delete [] d;
    }
    else
    {
        cout << "Inappropriate bit size." << endl;
        return;
    }

	delete [] randList;
	delete [] sortedList;
}

#endif /* SORTING_H_ */




