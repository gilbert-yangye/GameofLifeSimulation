//
//  MPI_game.cpp
//  Game_of_life
//
//  Created by Darren Shan on 2020/2/15.
//  Copyright © 2020 Darren Shan. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cstdlib>
#include <time.h>
using namespace std;

void read_size(string FileName, int &row, int &col, int &p_rows, int &p_cols)
{

    ifstream infile;
    infile.open(FileName);
    if (!infile)
        throw invalid_argument(FileName + " not open! Change the directory of the file!");

    infile >> p_rows;
    infile >> p_cols;
    infile >> row;
    infile >> col;
    infile.close();
}

void read_file(string FileName, int row, int col, int p_rows, int p_cols, bool *values)
{

    ifstream infile;
    infile.open(FileName);
    if (!infile)
        throw invalid_argument(FileName + " not open! Change the directory of the file!");

    infile >> p_cols;
    infile >> p_cols;
    infile >> p_cols;
    infile >> p_cols;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            infile >> values[i * col + j];
        }
    }
    infile.close();
}

void write_file (string outFileName, int row, int col, int p_rows, int p_cols, bool *values)
{
    ofstream outfile;
    outfile.open(outFileName);
    if (!outfile)
        throw invalid_argument(" " + outFileName + "not open! Change the directory of the file!");

    outfile << p_rows << " ";
    outfile << p_cols << endl;
    outfile << row << " ";
    outfile << col <<endl;
    for (int i = 0; i < row*col; i++)
    {
        outfile << values[i] << " ";
        if (i%col == col-1) outfile << endl;
    }
    outfile.close();
}

void find_neighbors(int columns,bool* data, int &neighbors,int rank)
{
    neighbors += data[rank-columns-1];
    neighbors += data[rank-columns];
    neighbors += data[rank-columns+1];
    neighbors += data[rank-1];
    neighbors += data[rank+1];
    neighbors += data[rank+columns-1];
    neighbors += data[rank+columns];
    neighbors += data[rank+columns+1];
}

bool alive_cell_prop(int columns,bool* data,int rank)
{
    int alive_nbs = 0;
    find_neighbors(columns, data, alive_nbs, rank);
    if (alive_nbs==2||alive_nbs==3)
    {
        return true;
    }
    else return false;
}

bool dead_cell_prop(int columns,bool* data,int rank)
{
    int alive_nbs = 0;
    find_neighbors(columns, data, alive_nbs, rank);
    if (alive_nbs==3)
    {
        return true;
    }
    else return false;
}

void update_life(int rows, int columns,bool* data, bool* new_data)
{
    for (int i = columns+1; i<rows*columns-columns; i++) //ignore the top and bottom rows
    {
        if ((i+1)%columns==0) // ignore the left and right columns
        {
            i = i+1;
        }
        else
        {
            if (data[i]) //if alive
            {
                new_data[i] = alive_cell_prop(columns, data, i);
            }
            else //else dead
            {
                new_data[i] = dead_cell_prop(columns, data, i);
            }
        }
    }
}

// for outer circle cells, we need to specify how to update with this function
bool update_custom(bool status, bool tl,bool t, bool tr, bool l, bool r, bool bl, bool b, bool br)
{
    int neighbors = 0;
    neighbors += tl;
    neighbors += t;
    neighbors += tr;
    neighbors += l;
    neighbors += r;
    neighbors += bl;
    neighbors += b;
    neighbors += br;
    if (status)
    {
        if (neighbors == 2|| neighbors == 3) return true;
        else return false;
    }
    else
    {
        if (neighbors == 3) return true;
        else return false;
    }

}

// for outer circle cells, we need to specify how to update with this function
void update_outer_life(int rows, int columns,bool* data, bool* new_data,bool* top_row, bool* bot_row, bool *left_col, bool *right_col, bool *top_left, bool *top_right,bool *bot_left,bool *bot_right)
{
    int jump = (rows-1)*columns; //how many elements before the last row: jump

    //update corners
    new_data[0] = update_custom(data[0], top_left[0], top_row[0], top_row[1], left_col[0], data[1], left_col[1], data[columns], data[columns+1]);
    new_data[rows*columns-1] = update_custom(data[rows*columns-1],data[(rows-1)*columns-2], data[(rows-1)*columns-1], right_col[rows-2], data[rows*columns-2], right_col[rows-1], bot_row[columns-2], bot_row[columns-1], bot_right[0]);
    new_data[columns-1] = update_custom(data[columns-1],top_row[columns-2], top_row[columns-1], top_right[0], data[columns-2], right_col[0], data[2*columns-2], data[2*columns-1], right_col[1]);
    new_data[jump] = update_custom(data[jump],left_col[rows-2], data[jump-columns], data[jump-columns+1], left_col[rows-1], data[jump+1], bot_left[0], bot_row[0], bot_row[1]);

    //update rows
    for (int i = 1; i<columns-1; i++)
    {
        new_data[i] = update_custom(data[i], top_row[i-1], top_row[i], top_row[i+1], data[i-1], data[i+1], data[columns+i-1], data[columns+i], data[columns+i+1]);
        new_data[i+jump] = update_custom(data[i+jump], bot_row[i-1], bot_row[i], bot_row[i+1], data[i+jump-1], data[i+jump+1], data[i+jump-columns-1], data[i+jump-columns], data[i+jump-columns+1]);
    }

    //updata columns
    for (int i = 1; i<rows-1; i++)
    {
        new_data[i*columns] = update_custom(data[i*columns], left_col[i-1], left_col[i], left_col[i+1], data[(i-1)*columns], data[(i+1)*columns], data[(i-1)*columns+1], data[i*columns+1], data[(i+1)*columns+1]);
        new_data[i*columns+columns-1] = update_custom(data[i*columns+columns-1], right_col[i-1], right_col[i], right_col[i+1], data[i*columns-1], data[i*columns-2], data[i*columns+2*columns-1], data[i*columns+2*columns-2], data[i*columns+columns-2]);
    }
}

//if non_periodic condition, some cells should regard some data as 0
void update_outer_life_non_periodic(int id, int p_rows, int p_cols, int rows, int columns,bool* data, bool* new_data,bool* top_row, bool* bot_row, bool *left_col, bool *right_col, bool *top_left, bool *top_right,bool *bot_left,bool *bot_right)
{
    int p_i = id/p_cols;
    int p_j = id%p_cols;
    
    int jump = (rows-1)*columns; //how many elements before the last row: jump

    //update corners!!!!
    if (p_i == 0 && p_j ==0)// top left processor
    {
        new_data[0] = update_custom(data[0], 0, 0, 0, 0, data[1], 0, data[columns], data[columns+1]);
        new_data[columns-1] = update_custom(data[columns-1],0, 0, 0, data[columns-2], right_col[0], data[2*columns-2], data[2*columns-1], right_col[1]);
        new_data[jump] = update_custom(data[jump],0, data[jump-columns], data[jump-columns+1], 0, data[jump+1], 0, bot_row[0], bot_row[1]);
        new_data[rows*columns-1] = update_custom(data[rows*columns-1],data[(rows-1)*columns-2], data[(rows-1)*columns-1], right_col[rows-2], data[rows*columns-2], right_col[rows-1], bot_row[columns-2], bot_row[columns-1], bot_right[0]);
    }
    else if (p_i == 0 && p_j ==p_cols-1) // top right processor
    {
        new_data[0] = update_custom(data[0], 0, 0, 0, left_col[0], data[1], left_col[1], data[columns], data[columns+1]);
        new_data[columns-1] = update_custom(data[columns-1], 0, 0, 0, data[columns-2], 0, data[2*columns-2], data[2*columns-1], 0);
        new_data[jump] = update_custom(data[jump],left_col[rows-2], data[jump-columns], data[jump-columns+1], left_col[rows-1], data[jump+1], bot_left[0], bot_row[0], bot_row[1]);
        new_data[rows*columns-1] = update_custom(data[rows*columns-1],data[(rows-1)*columns-2], data[(rows-1)*columns-1], 0, data[rows*columns-2], 0, bot_row[columns-2], bot_row[columns-1], 0);
    }
    else if (p_i == p_rows-1 && p_j ==0) // bot left processor
    {
        new_data[0] = update_custom(data[0], 0, top_row[0], top_row[1], 0, data[1], 0, data[columns], data[columns+1]);
        new_data[columns-1] = update_custom(data[columns-1],top_row[columns-2], top_row[columns-1], top_right[0], data[columns-2], right_col[0], data[2*columns-2], data[2*columns-1], right_col[1]);
        new_data[jump] = update_custom(data[jump],0, data[jump-columns], data[jump-columns+1], 0, data[jump+1], 0, 0, 0);
        new_data[rows*columns-1] = update_custom(data[rows*columns-1],data[(rows-1)*columns-2], data[(rows-1)*columns-1], right_col[rows-2], data[rows*columns-2], right_col[rows-1], 0, 0, 0);
    }
    else if (p_i == p_rows-1 && p_j ==p_cols-1) // bot right processor
    {
        new_data[0] = update_custom(data[0], top_left[0], top_row[0], top_row[1], left_col[0], data[1], left_col[1], data[columns], data[columns+1]);
        new_data[columns-1] = update_custom(data[columns-1],top_row[columns-2], top_row[columns-1], 0, data[columns-2], 0, data[2*columns-2], data[2*columns-1], 0);
        new_data[jump] = update_custom(data[jump],left_col[rows-2], data[jump-columns], data[jump-columns+1], left_col[rows-1], data[jump+1], 0, 0, 0);
        new_data[rows*columns-1] = update_custom(data[rows*columns-1],data[(rows-1)*columns-2], data[(rows-1)*columns-1], 0, data[rows*columns-2], 0, 0, 0, 0);
    }
    else //with normal iteration methods
    {
        new_data[0] = update_custom(data[0], top_left[0], top_row[0], top_row[1], left_col[0], data[1], left_col[1], data[columns], data[columns+1]);
        new_data[columns-1] = update_custom(data[columns-1],top_row[columns-2], top_row[columns-1], top_right[0], data[columns-2], right_col[0], data[2*columns-2], data[2*columns-1], right_col[1]);
        new_data[jump] = update_custom(data[jump],left_col[rows-2], data[jump-columns], data[jump-columns+1], left_col[rows-1], data[jump+1], bot_left[0], bot_row[0], bot_row[1]);
        new_data[rows*columns-1] = update_custom(data[rows*columns-1],data[(rows-1)*columns-2], data[(rows-1)*columns-1], right_col[rows-2], data[rows*columns-2], right_col[rows-1], bot_row[columns-2], bot_row[columns-1], bot_right[0]);
    }
    
    
    //update rows!!!!!
    for (int i = 1; i<columns-1; i++)
    {
        if (p_i == 0)
        {
            new_data[i] = update_custom(data[i], 0, 0, 0, data[i-1], data[i+1], data[columns+i-1], data[columns+i], data[columns+i+1]);
        }
        else
        {
            new_data[i] = update_custom(data[i], top_row[i-1], top_row[i], top_row[i+1], data[i-1], data[i+1], data[columns+i-1], data[columns+i], data[columns+i+1]);
        }
        
        if (p_i == p_rows-1)
        {
            new_data[i+jump] = update_custom(data[i+jump], 0, 0, 0, data[i+jump-1], data[i+jump+1], data[i+jump-columns-1], data[i+jump-columns], data[i+jump-columns+1]);
        }
        else
        {
            new_data[i+jump] = update_custom(data[i+jump], bot_row[i-1], bot_row[i], bot_row[i+1], data[i+jump-1], data[i+jump+1], data[i+jump-columns-1], data[i+jump-columns], data[i+jump-columns+1]);
        }
        
    }

    //updata columns !!!!!
    for (int i = 1; i<rows-1; i++)
    {
        if (p_j == 0)
        {
            new_data[i*columns] = update_custom(data[i*columns], 0, 0, 0, data[(i-1)*columns], data[(i+1)*columns], data[(i-1)*columns+1], data[i*columns+1], data[(i+1)*columns+1]);
        }
        else
        {
            new_data[i*columns] = update_custom(data[i*columns], left_col[i-1], left_col[i], left_col[i+1], data[(i-1)*columns], data[(i+1)*columns], data[(i-1)*columns+1], data[i*columns+1], data[(i+1)*columns+1]);
        }
        
        if (p_j == p_cols-1)
        {
             new_data[i*columns+columns-1] = update_custom(data[i*columns+columns-1], 0, 0, 0, data[i*columns-1], data[i*columns-2], data[i*columns+2*columns-1], data[i*columns+2*columns-2], data[i*columns+columns-2]);
        }
        else
        {
             new_data[i*columns+columns-1] = update_custom(data[i*columns+columns-1], right_col[i-1], right_col[i], right_col[i+1], data[i*columns-1], data[i*columns-2], data[i*columns+2*columns-1], data[i*columns+2*columns-2], data[i*columns+columns-2]);
        }
    }
}

// define to which processor the information should go, and to which it should receive
void set_connection (int id, int p_rows, int p_cols, int *connect)
{
    //0 top left
    if (id == 0) connect[0] = p_rows*p_cols-1;
    else if (id%p_cols == 0) connect[0] = id  - 1;
    else if (id < p_cols) connect[0] = (p_rows-1)*p_cols+id-1;
    else connect[0] = id - p_cols - 1;

    //1 top row
    if (id - p_cols>=0) connect[1] = id - p_cols;
    else connect[1] = (p_rows-1)*p_cols+id;

    //2 top right
    if (id == p_cols-1) connect[2] = (p_rows-1)*p_cols;
    else if ((id+1)%p_cols==0) connect[2] = id + 1 - 2*p_cols;
    else if (id < p_cols) connect[2] = (p_rows-1)*p_cols+id+1;
    else connect[2] = id - p_cols + 1;

    //3 left col
    if (id%p_cols!=0) connect[3] = id -1;
    else connect[3] = id + p_cols - 1;

    //4 right col
    if ((id+1)%p_cols!=0) connect[4] = id +1;
    else connect[4] = id - p_cols + 1;

    //5 bot left
    if (id == (p_rows-1)*p_cols) connect[5] = p_cols - 1;
    else if (id%p_cols==0) connect[5] =  id + 2*p_cols - 1;
    else if (id >= (p_rows-1)*p_cols) connect[5] =  id-(p_rows-1)*p_cols-1;
    else connect[5] =  id + p_cols - 1;

    //6 bot row
    if (id + p_cols < p_rows*p_cols) connect[6] = id + p_cols;
    else connect[6] = id - (p_rows-1)*p_cols;

    //7 bot right
    if (id == p_rows*p_cols-1) connect[7] = 0;
    else if ((id+1)%p_cols==0) connect[7] = id + 1;
    else if (id >= (p_rows-1)*p_cols) connect[7] = id-(p_rows-1)*p_cols+1;
    else connect[7] = id + p_cols + 1;
}

// use mpi type to send the discontinuous column data
void buildMPIType(int rows, int cols, bool *value, MPI_Datatype &new_type)
{
    int block_lengths[rows];
    MPI_Aint displacements[rows];
    MPI_Aint add_start;
    MPI_Datatype typelist[rows];

    for (int i = 0; i<rows; i++)
    {
        typelist[i] = MPI_C_BOOL;
        block_lengths[i] = 1;
    }

    MPI_Get_address(value, &add_start);

    displacements[0] = 0;
    for (int i = 1; i < rows; i++) displacements[i] = displacements[i-1] + cols; // displacement = number of columns

    MPI_Type_create_struct(rows, block_lengths, displacements, typelist, &new_type);
    MPI_Type_commit(&new_type);
}

int id, p;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int itr = atoi(argv[1]); // get number of iterations and boundary condition from execution command line
    bool periodic = 1;
    periodic = atoi(argv[2]);
    
    int p_rows, p_cols, rows, columns;
    string Filename_in = "splitdata"+to_string(id)+".txt";
    string Filename_out = "testdata_result"+to_string(id);
    string directory = "../Game_of_life/testdata/";
    Filename_in = directory+Filename_in;
    Filename_out = directory + Filename_out;
    


    
    // read the divided data, served by pre_processing.cpp
    //firstly read the size, to initialise the dataset
    read_size(Filename_in, rows, columns,p_rows,p_cols);
    bool *data_1 = new bool[rows*columns];
    bool *data_2 = new bool[rows*columns];


    bool* top_row_rcv = new bool[columns];
    bool* bot_row_rcv = new bool[columns];
    bool *left_col_rcv= new bool[rows];
    bool *right_col_rcv= new bool[rows];
    bool *top_left_rcv = new bool[1];
    bool *top_right_rcv = new bool[1];
    bool *bot_left_rcv = new bool[1];
    bool *bot_right_rcv = new bool[1];

    //read original data ---> data_1
    read_file(Filename_in, rows, columns,p_rows,p_cols, data_1);

    bool *data_this, *data_nxt; //data_this always points at this generation; data_nxt always points at next generation
    
    MPI_Datatype MPI_type1l;
    MPI_Datatype MPI_type1r;
    buildMPIType(rows, columns, data_this, MPI_type1l);
    buildMPIType(rows, columns, &data_this[columns-1], MPI_type1r);

    //set up connectivity
    int* connect = new int[8];
    set_connection(id, p_rows, p_cols, connect);

    int tag_num = 1;

    MPI_Request* request = new MPI_Request[16*itr];

    for (int i = 0; i<itr; i++)
    {
        if (i%2==0)//when even, data_1 is this generation
        {
            data_this = &data_1[0];
            data_nxt = &data_2[0];
        }
        else//when odd, data_2 is this generation
        {
            data_this = &data_2[0];
            data_nxt = &data_1[0];
        }
        
        if (p_cols == 1 && p_rows == 1) // deal with single processor
        {
            //copy all the information if only 1 processor involved;
            for (int n = 0; n<rows; n++) left_col_rcv[n] = data_this[(n+1)*columns-1];
            for (int n = 0; n<rows; n++) right_col_rcv[n] = data_this[n*columns];
            for (int n = 0; n<columns; n++) top_row_rcv[n] = data_this[n];
            for (int n = 0; n<columns; n++) bot_row_rcv[n] = data_this[(rows-1)*columns+n];
            top_left_rcv[0] = data_this[rows*columns-1];
            top_right_rcv[0] = data_this[(rows-1)*columns];
            bot_left_rcv[0] = data_this[columns-1];
            bot_right_rcv[0] = data_this[0];
            
            update_life(rows, columns, data_this, data_nxt);
        }
        else if (p_cols == 1) // deal with small number of processors that has only 1 column (2*1)
        {
            //copy column information if only one cols of processors
            MPI_Irecv(&top_left_rcv[0], 1, MPI_C_BOOL, connect[0], tag_num+7, MPI_COMM_WORLD, &request[i*16]);
            MPI_Irecv(&top_row_rcv[0], columns, MPI_C_BOOL, connect[1], tag_num+6, MPI_COMM_WORLD, &request[i*16+1]);
            MPI_Irecv(&top_right_rcv[0], 1, MPI_C_BOOL, connect[2], tag_num+5, MPI_COMM_WORLD, &request[i*16+2]);
            MPI_Irecv(&bot_left_rcv[0], 1, MPI_C_BOOL, connect[5], tag_num+2, MPI_COMM_WORLD, &request[i*16+3]);
            MPI_Irecv(&bot_row_rcv[0], columns, MPI_C_BOOL, connect[6], tag_num+1, MPI_COMM_WORLD, &request[i*16+4]);
            MPI_Irecv(&bot_right_rcv[0], 1, MPI_C_BOOL, connect[7], tag_num, MPI_COMM_WORLD, &request[i*16+5]);


            MPI_Isend(&data_this[0], 1, MPI_C_BOOL, connect[0], tag_num, MPI_COMM_WORLD, &request[i*16+6]);
            MPI_Isend(&data_this[0], columns, MPI_C_BOOL, connect[1], tag_num+1, MPI_COMM_WORLD, &request[i*16+7]);
            MPI_Isend(&data_this[columns-1], 1, MPI_C_BOOL, connect[2], tag_num+2, MPI_COMM_WORLD, &request[i*16+8]);
            MPI_Isend(&data_this[(rows-1)*columns], 1, MPI_C_BOOL, connect[5], tag_num+5, MPI_COMM_WORLD, &request[i*16+9]);
            MPI_Isend(&data_this[(rows-1)*columns], columns, MPI_C_BOOL, connect[6], tag_num+6, MPI_COMM_WORLD, &request[i*16+10]);
            MPI_Isend(&data_this[rows*columns-1], 1, MPI_C_BOOL, connect[7], tag_num+7, MPI_COMM_WORLD, &request[i*16+11]);
            tag_num+=8;

            //do operations when waiting for data transmission;
            //update the life status onto nxt generation;
            //aft 1st iteration, data_nxt (data_2) contains new info. so in 2nd itr, data_this should point at data_2;
            update_life(rows, columns, data_this, data_nxt);
            
            for (int n = 0; n<rows; n++) left_col_rcv[n] = data_this[(n+1)*columns-1];
            for (int n = 0; n<rows; n++) right_col_rcv[n] = data_this[n*columns];
            
            //wait for everything done
            MPI_Waitall(12, &request[i*16], MPI_STATUS_IGNORE);
        }
        else if (p_rows == 1) // deal with small number of processors that has only 1 row (1*2)
        {
            // if only one row of processors, copy the rows information
            MPI_Irecv(&top_left_rcv[0], 1, MPI_C_BOOL, connect[0], tag_num+7, MPI_COMM_WORLD, &request[i*16]);
            MPI_Irecv(&top_right_rcv[0], 1, MPI_C_BOOL, connect[2], tag_num+5, MPI_COMM_WORLD, &request[i*16+1]);
            MPI_Irecv(&left_col_rcv[0], rows, MPI_C_BOOL, connect[3], tag_num+4, MPI_COMM_WORLD, &request[i*16+2]);
            MPI_Irecv(&right_col_rcv[0], rows, MPI_C_BOOL, connect[4], tag_num+3, MPI_COMM_WORLD, &request[i*16+3]);
            MPI_Irecv(&bot_left_rcv[0], 1, MPI_C_BOOL, connect[5], tag_num+2, MPI_COMM_WORLD, &request[i*16+4]);
            MPI_Irecv(&bot_right_rcv[0], 1, MPI_C_BOOL, connect[7], tag_num, MPI_COMM_WORLD, &request[i*16+5]);


            MPI_Isend(&data_this[0], 1, MPI_C_BOOL, connect[0], tag_num, MPI_COMM_WORLD, &request[i*16+6]);
            MPI_Isend(&data_this[columns-1], 1, MPI_C_BOOL, connect[2], tag_num+2, MPI_COMM_WORLD, &request[i*16+7]);
            MPI_Isend(&data_this[0], 1, MPI_type1l, connect[3], tag_num+3, MPI_COMM_WORLD, &request[i*16+8]);
            MPI_Isend(&data_this[columns-1], 1, MPI_type1r, connect[4], tag_num+4, MPI_COMM_WORLD, &request[i*16+9]);
            MPI_Isend(&data_this[(rows-1)*columns], 1, MPI_C_BOOL, connect[5], tag_num+5, MPI_COMM_WORLD, &request[i*16+10]);
            MPI_Isend(&data_this[rows*columns-1], 1, MPI_C_BOOL, connect[7], tag_num+7, MPI_COMM_WORLD, &request[i*16+11]);
            tag_num+=8;

            //do operations when waiting for data transmission;
            //update the life status onto nxt generation;
            //aft 1st iteration, data_nxt (data_2) contains new info. so in 2nd itr, data_this should point at data_2;
            update_life(rows, columns, data_this, data_nxt);

            for (int n = 0; n<columns; n++) top_row_rcv[n] = data_this[n];
            for (int n = 0; n<columns; n++) bot_row_rcv[n] = data_this[(rows-1)*columns+n];
            
            //wait for everything done
            MPI_Waitall(12, &request[i*16], MPI_STATUS_IGNORE);
        }
        else // for normal condition, where both p_rows and p_cols larger than or equals to 2;
        {
            MPI_Irecv(&top_left_rcv[0], 1, MPI_C_BOOL, connect[0], tag_num+7, MPI_COMM_WORLD, &request[i*16]);
            MPI_Irecv(&top_row_rcv[0], columns, MPI_C_BOOL, connect[1], tag_num+6, MPI_COMM_WORLD, &request[i*16+1]);
            MPI_Irecv(&top_right_rcv[0], 1, MPI_C_BOOL, connect[2], tag_num+5, MPI_COMM_WORLD, &request[i*16+2]);
            MPI_Irecv(&left_col_rcv[0], rows, MPI_C_BOOL, connect[3], tag_num+4, MPI_COMM_WORLD, &request[i*16+3]);
            MPI_Irecv(&right_col_rcv[0], rows, MPI_C_BOOL, connect[4], tag_num+3, MPI_COMM_WORLD, &request[i*16+4]);
            MPI_Irecv(&bot_left_rcv[0], 1, MPI_C_BOOL, connect[5], tag_num+2, MPI_COMM_WORLD, &request[i*16+5]);
            MPI_Irecv(&bot_row_rcv[0], columns, MPI_C_BOOL, connect[6], tag_num+1, MPI_COMM_WORLD, &request[i*16+6]);
            MPI_Irecv(&bot_right_rcv[0], 1, MPI_C_BOOL, connect[7], tag_num, MPI_COMM_WORLD, &request[i*16+7]);


            MPI_Isend(&data_this[0], 1, MPI_C_BOOL, connect[0], tag_num, MPI_COMM_WORLD, &request[i*16+8]);
            MPI_Isend(&data_this[0], columns, MPI_C_BOOL, connect[1], tag_num+1, MPI_COMM_WORLD, &request[i*16+9]);
            MPI_Isend(&data_this[columns-1], 1, MPI_C_BOOL, connect[2], tag_num+2, MPI_COMM_WORLD, &request[i*16+10]);
            //use mpitype to send column data;
            MPI_Isend(&data_this[0], 1, MPI_type1l, connect[3], tag_num+3, MPI_COMM_WORLD, &request[i*16+11]);
            MPI_Isend(&data_this[columns-1], 1, MPI_type1r, connect[4], tag_num+4, MPI_COMM_WORLD, &request[i*16+12]);
            
            MPI_Isend(&data_this[(rows-1)*columns], 1, MPI_C_BOOL, connect[5], tag_num+5, MPI_COMM_WORLD, &request[i*16+13]);
            MPI_Isend(&data_this[(rows-1)*columns], columns, MPI_C_BOOL, connect[6], tag_num+6, MPI_COMM_WORLD, &request[i*16+14]);
            MPI_Isend(&data_this[rows*columns-1], 1, MPI_C_BOOL, connect[7], tag_num+7, MPI_COMM_WORLD, &request[i*16+15]);
            tag_num+=8;

            //do operations when waiting for data transmission;
            //update the life status onto nxt generation;
            //aft 1st iteration, data_nxt (data_2) contains new info. so in 2nd itr, data_this should point at data_2;
            update_life(rows, columns, data_this, data_nxt);

            //wait for everything done
            MPI_Waitall(16, &request[i*16], MPI_STATUS_IGNORE);
        }


        //do the iteration for outside;
        // use different function for periodic and non periodic functions
        if (periodic) update_outer_life(rows, columns, data_this, data_nxt, top_row_rcv, bot_row_rcv, left_col_rcv, right_col_rcv, top_left_rcv, top_right_rcv, bot_left_rcv, bot_right_rcv);
        else update_outer_life_non_periodic(id, p_rows, p_cols, rows, columns, data_this, data_nxt, top_row_rcv, bot_row_rcv, left_col_rcv, right_col_rcv, top_left_rcv, top_right_rcv, bot_left_rcv, bot_right_rcv);


        string temp;
        temp= Filename_out+"_"+to_string(i)+".txt";
        write_file(temp, rows, columns, p_rows, p_cols, data_nxt);
        cout<<"Processor "<<id<<" finished writing data in iteration "<<i<<endl;
    }


    delete[] top_row_rcv;
    delete[] bot_row_rcv;
    delete[] left_col_rcv;
    delete[] right_col_rcv;
    delete[] top_left_rcv;
    delete[] top_right_rcv;
    delete[] bot_left_rcv;
    delete[] bot_right_rcv;

    delete[] data_1;
    delete[] data_2;
    delete[] request;
    delete[] connect;

    MPI_Finalize();
    return 0;
}
