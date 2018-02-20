#include <vector>

/*
    Created by: Steven Ramsey, PhD Candidate Biochemistry CUNY Graduate Center
    Contact information: vpsramsey@gmail.com
                    Room 1409 Science Hall Lehman College, 250 Bedford Park blvd Bronx NY, 10468
    Mentor: Tom Kurtzman

*/

struct dx {
    /*
        Struct contains the main functions and values associated with a particular dx file
    */
    int count[3];
    double origin[3];
    double delta[3];
    int totalpoints;
    int numgrpatt;
    std::vector<double > data;


    dx(); //initializer which will put it to 0
    void readDx(std::string infile); //properly reads in paramaters from dx file
    void writeDx(std::string outfile); //properly writes dx file
    bool checkSame(dx test); //Necessary so that when two dx files are used in a function their voxel positions in an array are equivalent
    dx clearByOne(double cutoff, char flag);
    dx clearByTwo(dx E, double cutoff1, double  cutoff2, char flag1, char flag2);
    void contour(dx G, double cutoff, char flag); //works the same as clearbyone, but does not copy
    void makeGroups();
    dx solventAccessible(); //input must be a g(r) dx file
    int checkNeighbor(int index, std::vector<int > &grpNum, int crnum);
    void mult(dx G);
    double sum();
    void printPdb(std::string outfile);
    void add(dx G);
    void addbyConst(double C);
    void multbyConst(double C);
    void sub(dx G);
    void div(dx G);
    double avg();
    void setBP(std::vector<double > &ha1, std::vector<double > &ha2, std::vector<double > &ha3, double D);
    void setHeavi(std::vector<double > &ha1, std::vector<double > &ha2, std::vector<double > &ha3, double D);
    void calcVdw(double vdw);
    void write_out_dx(std::string infile, int column);
    void histogram(std::string outfile);
    void printcol(std::string outfile);
    dx cat(dx N); //concatenate dx files self and N to create a new dx file.
    void citation();
    dx zeros();

};


struct lig {
    /*
        Structure merely contains array of x, y, z coordinates (separated) and a singular function to set the values
    */
    std::vector<double > hax;
    std::vector<double > hay;
    std::vector<double > haz;
    void readLF(std::string ligfile);
};

struct pop {
    /*
        struct purpose is to read in the gist out file and perform the specific calculations of interest
    */
    std::vector<int > N;
    std::vector<double > transnorm;
    std::vector<double > orientnorm;
    std::vector<double > Eswnorm;
    std::vector<double > Ewwnorm;
    std::string temp;
    double tempx;
    double average;
    double avgtnorm;
    double avgonorm;
    double avgEwwnorm;
    double avgEswnorm;
    pop();
    void calcpop(std::string gistfile); //where gist file is the column deliminited data file produced by gist output
};



