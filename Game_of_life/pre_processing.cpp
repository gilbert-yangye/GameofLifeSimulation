//
//  main_0_processor.cpp
//  Game_of_life
//
//  Created by Darren Shan on 2020/2/13.
//  Copyright © 2020 Darren Shan. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;


 //using grid division;
void find_dimensions(int p, int &p_rows, int &p_columns)
{
    int top = sqrt(p) + 1;
    for (int i = top; i <= top; i--)
    {
        if (p%i == 0)
        {
            p_rows = i;
            p_columns = p / i;
            break;
        }
    }
}

//sort datasets in the order of processors
void assign_datasets(int p_rows, int p_cols, int* rows, int* cols, int total_rows, int total_cols, bool * values, bool * sorted_values)
{
    int temp_rows = total_rows, temp_cols= total_cols;
    for (int i = 0; i < p_rows; i++)
    {
        rows[i] = int (temp_rows/(p_rows-i));
        temp_rows -= rows[i];
    }
    for (int i = 0; i < p_cols; i++)
    {
        cols[i] = int (temp_cols/(p_cols-i));
        temp_cols -= cols[i];
    }
    
    
    int jump = 0;
    int org_row = 0;
    int org_col = 0;
    int row(0), col(0);
    for (int proc =0; proc<p_cols*p_rows; proc++)
    {
        for (int i = 0; i < rows[row]; i++)
        {
            for (int j = 0; j < cols[col]; j++)
            {
                sorted_values[i*cols[col]+j+jump] = values[(i+org_row)*total_cols+j+org_col];
            }
        }
        jump+=cols[col]*rows[row];
        if ((proc+1)%p_cols==0)
        {
            org_row += rows[row];
            org_col = 0;
            
            row += 1;
            col = 0;
        }
        else
        {
            org_col += cols[col];
            col+=1;
        }
    }
    
}

void pre_read_file(string FileName, int m, int n, bool *values)
{
    
    ifstream infile;
    infile.open(FileName);
    if (!infile)
        throw invalid_argument(FileName + " not open! Change the directory of the file!");
    
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            infile >> values[i * n + j];
        }
    }
    infile.close();
}

void pre_write_file (string outFileName, int m, int n, int jump, int p_rows, int p_cols,bool *values)
{
    ofstream outfile;
    outfile.open(outFileName);
    if (!outfile)
        throw invalid_argument(" " + outFileName + "not open! Change the directory of the file!");
    
    outfile << p_rows << " ";
    outfile << p_cols << endl;
    outfile << m << " ";
    outfile << n << endl;
    for (int i = 0; i < m*n; i++)
    {
        outfile << values[i+jump] << " ";
        if (i%n == n-1) outfile << endl;
        
    }
    outfile.close();
}


int main()
{
    int p = 4;
    int rows = 36;
    int columns = 36;
    int p_rows, p_cols;
    find_dimensions(p, p_rows, p_cols);
    
    string directory = "../Game_of_life/testdata/";
    string file_in = "testdata.txt";
    file_in = directory + file_in;
    bool* values = new bool[rows*columns];
    bool* sorted_values = new bool[rows*columns];
    
    pre_read_file(file_in, rows, columns, values);
    
//    or generate from random
//    for (int i = 0; i<rows*columns; i++) values[i] = rand()%2;
    
    int *rows_per_p, *cols_per_p;
    rows_per_p = new int[p_rows];
    cols_per_p = new int[p_cols];
    
    assign_datasets(p_rows, p_cols, rows_per_p, cols_per_p, rows, columns, values, sorted_values);
    
    string file_out = "splitdata";
    string temp;
    int jump = 0;

    for (int i = 0; i<p; i++)
    {
        temp = directory + file_out + to_string(i) +".txt";
        
        pre_write_file(temp, rows_per_p[int (i/p_cols)], cols_per_p[i%p_cols], jump, p_rows, p_cols,sorted_values);
        jump+=rows_per_p[int (i/p_cols)]*cols_per_p[i%p_cols];
        
        cout<<"finished writing "<<i<<endl;
    }
    
    
    delete[] values;
    delete[] sorted_values;
    delete[] rows_per_p;
    delete[] cols_per_p;
    return 0;
}
