#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <string.h>
#include <vector>
#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include "gistpp.h"

/*
    Created by: Steven Ramsey, PhD Candidate Biochemistry CUNY Graduate Center
    Contact information: vpsramsey@gmail.com
                    Room 1409 Science Hall Lehman College, 250 Bedford Park blvd Bronx NY, 10468
    Mentor: Tom Kurtzman

*/

using namespace std;

void lig::readLF(string ligfile) {
    /*
        Function reads in the ligand.pdb file line by line taking x y z coordinates
         from the correct column, ignoring rows which contain hydrogen data
    */
    string temp; //storage
    double tempx; //storage
    ifstream input(ligfile.c_str());
    if (!input.is_open()) {
	cerr << "Could not open file " << ligfile << " please check it is available.\n";
        exit(0);
    }
    //ofstream output("ligandcoord.txt"); //test file to be sure ligand is read properly
    while (!input.eof()) {
        getline(input, temp);
        if (!temp.empty()){
            if (!strcmp(temp.substr(0,1).c_str(), "T")) {continue;} //skip TER and TITLE lines
            else if (!strcmp(temp.substr(0,1).c_str(), "R")) {continue;} //skip REMARK lines
            else if (!strcmp(temp.substr(0,1).c_str(), "E")) {continue;} //skip END lines
            else if (!strcmp(temp.substr(0,1).c_str(), "C")) {continue;} //skip CRYSTAL lines
	    else if (!strcmp(temp.substr(0,1).c_str(), "M")) {continue;} //skip MODEL lines
            else {
                //checked for TER, END, CRYSTAL, and REMARK flags
                if (!strcmp(temp.substr(13, 1).c_str(), "H") || !strcmp(temp.substr(13,1).c_str(), "h")) {continue;} //do not both with hydrogen data
                else {
                    //got rid of the hydrogens
                    tempx = atof(temp.substr(31, 7).c_str());
                    hax.push_back(tempx);
                    tempx = atof(temp.substr(39, 7).c_str());
                    hay.push_back(tempx);
                    tempx = atof(temp.substr(47, 7).c_str());
                    haz.push_back(tempx);
                }
            }
        }
        else {continue;}
    }
    //test prints
    //for (int i = 0; i < hax.size(); i++) {
    //    output << hax[i] << "\t" << hay[i] << "\t" << haz[i] << endl;
    //}
}


void dx::setBP(vector<double > &ha1, vector<double > &ha2, vector<double > &ha3, double D) {
    /*
        Function goes through dx file voxel by voxel and checks the voxels distance compared to each heavy atom in the ligand
    */
    double vals[3]; //storage for voxel x y z
    double dd = 0; //storage for calculated distance
    double DD = D*D; //squared distances since why do sqrt
    int pos = 0;
    for (int i = 0; i < count[0]; i++) {
        for (int j = 0; j < count[1]; j++) {
            for (int k = 0; k < count[2]; k++) {
                vals[0] = i*delta[0] + origin[0]; //find x value of voxel
                vals[1] = j*delta[1] + origin[1]; //find y value of voxel
                vals[2] = k*delta[2] + origin[2]; //find z value of voxel
                for (int a = 0; a < ha1.size(); a++) { //check each heavy atom distance
                    dd = pow(vals[0]-ha1[a], 2) + pow(vals[1] - ha2[a], 2) + pow(vals[2] - ha3[a], 2); //calc distance
                    if (dd <= DD){ //&& data[pos] > 0) {
                        //cout << "voxel in range:\n" << vals[0] << "\t" << vals[1] << "\t" << vals[2] << "\t" << ha1[a] << "\t" << ha2[a] << "\t" << ha3[a] << endl;
                        data[pos] = 1; break; //if within range stop checking and move on
                    }
                    else {data[pos] = 0;}
                }
                pos++;
            }
        }
    }



}

void dx::setHeavi(vector<double > &ha1, vector<double > &ha2, vector<double > &ha3, double D) {
    /*
        Function goes through dx file voxel by voxel and checks the voxels distance compared to each heavy atom in the ligand
    */
    double vals[3]; //storage for voxel x y z
    double dd = 0; //storage for calculated distance
    double DD = D*D; //squared distances since why do sqrt
    int pos = 0;
    for (int i = 0; i < count[0]; i++) {
        for (int j = 0; j < count[1]; j++) {
            for (int k = 0; k < count[2]; k++) {
                data[pos] = 0;
                pos++;
            }
        }
    }
    pos = 0;
    double heavi = 0;
    for (int i = 0; i < count[0]; i++) {
        for (int j = 0; j < count[1]; j++) {
            for (int k = 0; k < count[2]; k++) {
                vals[0] = i*delta[0] + origin[0]; //find x value of voxel
                vals[1] = j*delta[1] + origin[1]; //find y value of voxel
                vals[2] = k*delta[2] + origin[2]; //find z value of voxel
                for (int a = 0; a < ha1.size(); a++) { //check each heavy atom distance
                    dd = pow(vals[0]-ha1[a], 2) + pow(vals[1] - ha2[a], 2) + pow(vals[2] - ha3[a], 2); //calc distance
                    if (dd <= DD){ //&& data[pos] > 0) {
                        //cout << "voxel in range:\n" << vals[0] << "\t" << vals[1] << "\t" << vals[2] << "\t" << ha1[a] << "\t" << ha2[a] << "\t" << ha3[a] << endl;
                        heavi = D-sqrt(dd);
                        heavi = heavi/D;
                        if (heavi > data[pos]) {data[pos] = heavi;}
                        //data[pos] = 1; break; //if within range stop checking and move on
                        //can't move on once resetting, if an atom is closer it takes priority
                        
                    }
                }
                pos++;
            }
        }
    }



}


dx::dx() {
    //default constructor, sets to 0's, except data which will remain un-initialized
    for (int i = 0; i < 3; i++) {
        count[i] = 0;
        origin[i] = 0.0;
        delta[i] = 0.0;
        totalpoints = 0;
    }
    numgrpatt = 0;
}


void dx::readDx(string infile) {
    string temp1;
    double tempx;
    string temp[44];
    int head_pos = 0;
    int data_count = 0;
    bool check = true;
    ifstream input(infile.c_str());
    if (!input.is_open()) {
    	cerr << "Could not find file: " << infile << " please ensure it is available.\n";
	exit(0);
    }
    while (!input.eof()) {
        if (head_pos < 44) {
            input >> temp1;
            if (temp1[0] != '#') { //checks for default comment type in dx files
                temp[head_pos] = temp1;
                head_pos++;
            }
            /*
                At this point we've read in the header in word order so that it can be deciphered.
            */
        }
        else {
            if (check) {
            	totalpoints = atoi(temp[41].c_str());
                check = false;
            }
            if (data_count == totalpoints) {
   		break;
            }
            input >> tempx;
	    //if (temp1 == "object") {
	    //break;
	    //}
            data.push_back(tempx);
            data_count++;
	    
        }
    }
    //Set the variables based on their position in the header chunk of a dx file
    count[0] = atoi(temp[5].c_str()); count[1] = atoi(temp[6].c_str()); count[2] = atoi(temp[7].c_str());
    origin[0] = atof(temp[9].c_str()); origin[1] = atof(temp[10].c_str()); origin[2] = atof(temp[11].c_str());
    delta[0] = atof(temp[13].c_str()); delta[1] = atof(temp[18].c_str()); delta[2] = atof(temp[23].c_str());
    //totalpoints = atoi(temp[41].c_str());

    //data.pop_back();

    if (totalpoints != data.size()) {
        cerr << "\nERROR!:\n\n"
        << "header: " << totalpoints << endl
        << "file: " << data.size() << endl
        << "last file data: " << data[data.size()+1] << endl
        << "first file data: " << data[0] << endl
        << "total points from header does not match number of data points acquired!\n\n";
        exit(0);
    }
}

void dx::writeDx(string outfile) {
    //make header then print in a 3 row matrix the data
    ofstream output(outfile.c_str());
    output << "object 1 class gridpositions counts " << count[0] << " " << count[1] << " " << count[2] << "\n";
    output << "origin " << origin[0] << " " << origin[1] << " " << origin[2] << "\n";
    output << "delta " << delta[0] << " 0 0\n";
    output << "delta " << "0 " << delta[1] << " 0\n";
    output << "delta " << "0 0 " << delta[2] << "\n";
    output << "object 2 class gridconnections counts " << count[0] << " " << count[1] << " " << count[2] << "\n";
    output << "object 3 class array type float rank 0 items " << totalpoints << " data follows\n";
    int a = 0;
    output.precision(16);
    for (int i = 0; i < totalpoints; i++) {
        if (a < 2) {
            output << data[i] << " ";
            a++;
        }
        else {
            output << data[i] << "\n";
            a = 0;
        }
    }
}


bool dx::checkSame(dx test) {
    //function will check to see if the two dx files have the same header values
    bool ans = true; //we want them to be the same, so this code will attempt to disprove the sameness
    if (count[0] != test.count[0]) {
        ans = false;
        cout << "count 0 doesn't match\n";
        return ans;
    }
    if (count[1] != test.count[1]) {
        ans = false;
        cout << "count 1 doesn't match\n";
        return ans;
    }
    if (count[2] != test.count[2]) {
        ans = false;
        cout << "count 2 doesn't match\n";
        return ans;
    }
    if (origin[0] != test.origin[0]) {
        ans = false;
        cout << "origin 0 doesn't match\n";
        return ans;
    }
    if (origin[1] != test.origin[1]) {
        ans = false;
        cout << "origin 1 doesn't match\n";
        return ans;
    }
    if (origin[2] != test.origin[2]) {
        ans = false;
        cout << "origin 2 doesn't match\n";
        return ans;
    }
    if (delta[0] != test.delta[0]) {
        ans = false;
        cout << "delta 0 doesn't match\n";
        return ans;
    }
    if (delta[1] != test.delta[1]) {
        ans = false;
        cout << "delta 1 doesn't match\n";
        return ans;
    }
    if (delta[2] != test.delta[2]) {
        ans = false;
        cout << "delta 2 doesn't match\n";
        return ans;
    }
    if (totalpoints != test.totalpoints) {
        ans = false;
        cout << "totalpoints doesn't match\n";
        return ans;
    }
    return ans;
}

dx dx::solventAccessible() {
    //No longer used definitions of neighbors
    int addz = 1;
    int subz = -1*1;
    int addy = count[2];
    int suby = -1*count[2];
    int addx = count[2]*count[1];
    int subx = -1*count[2]*count[1];
    int exc = 0;


    dx SASA;
    for (int i = 0; i < 3; i++) {
        SASA.count[i] = count[i];
        SASA.origin[i] = origin[i];
        SASA.delta[i] = delta[i];
    }
    SASA.totalpoints = totalpoints;
    cout << "\nSolvent Accessible surface representation to be calculated, input dx file must be representative of g(O)\n\n";
    double high = 0.3;
    double low = 0.1;
    int numplane = 0;


    for(int i = 0; i < totalpoints; i++) {
        numplane = i/(count[2]*count[1]);
        if (data[i] > high) {
            try {
                //if (i + addz > totalpoints || i+addz < 0 /*|| i == count[2]*/) {throw exc;}
                if (i%count[2] == count[2]-1) {throw exc;}
                else {
                    if (data[i+addz] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + addy > totalpoints || i+addy < 0 /*|| (i > count[1]*(count[2]-1) && i < count[1]*count[2])*/) {throw exc;}
                if (i%(count[2]*(count[1]-1)+(numplane*count[2]*count[1]))< count[2]) {throw exc;}
                else {
                    if (data[i+addy] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + addx > totalpoints || i+addx < 0) {throw exc;}
                if (i > count[2]*count[1]*(count[0]-1) && i < count[2]*count[1]*count[0]) {throw exc;}
                else {
                    if (data[i+addx] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + subz > totalpoints || i+subz < 0) {throw exc;}
                if (i%count[2] == 0) {throw exc;}
                else {
                    if (data[i+subz] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + suby > totalpoints || i+suby < 0) {throw exc;}
                if (i%(count[2]*count[1]) < count[2]) {throw exc;}
                else {
                    if (data[i+suby] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + subx > totalpoints || i+subx < 0) {throw exc;}
                if (i >= 0 && i < count[2]*count[1]) {throw exc;}
                else {
                    if (data[i+subx] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + addz + addy > totalpoints || i+addz+addy < 0) {throw exc;}
                if ((i%count[2] == count[2]-1) || (i%(count[2]*(count[1]-1)+(numplane*count[2]*count[1])) < count[2])) {throw exc;}
                else {
                    if (data[i+addz+addy] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + addz + suby > totalpoints || i+addz+suby < 0) {throw exc;}
                if ((i%count[2] == count[2]-1)||(i%(count[2]*count[1]) < count[2])) {throw exc;}
                else {
                    if (data[i+addz+suby] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + subz + addy > totalpoints || i+subz+addy < 0) {throw exc;}
                if ((i%count[2] == 0)||(i%(count[2]*(count[1]-1)+(numplane*count[2]*count[1]))< count[2])) {throw exc;}
                else {
                    if (data[i+subz+addy] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + subz + suby > totalpoints || i+subz+suby < 0) {throw exc;}
                if ((i%count[2] == 0)||(i%(count[2]*count[1]) < count[2])) {throw exc;}
                else {
                    if (data[i+subz+suby] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + addz + addx > totalpoints || i+addz+addx < 0) {throw exc;}
                if ((i%count[2] == count[2]-1)||(i > count[2]*count[1]*(count[0]-1) && i < count[2]*count[1]*count[0])) {throw exc;}
                else {
                    if (data[i+addz+addx] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + addz + subx > totalpoints || i+addz+subx < 0) {throw exc;}
                if ((i%count[2] == count[2]-1)||(i >= 0 && i < count[2]*count[1])) {throw exc;}
                else {
                    if (data[i+addz+subx] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + subz + addx > totalpoints || i+subz+addx < 0) {throw exc;}
                if ((i%count[2] == 0)||(i > count[2]*count[1]*(count[0]-1) && i < count[2]*count[1]*count[0])) {throw exc;}
                else {
                    if (data[i+subz+addx] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + subz + subx > totalpoints || i+subz+subx < 0) {throw exc;}
                if ((i%count[2] == 0)||(i >= 0 && i < count[2]*count[1])) {throw exc;}
                else {
                    if (data[i+subz+subx] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + addy + addx > totalpoints || i+addy+addx < 0) {throw exc;}
                if ((i%(count[2]*(count[1]-1)+(numplane*count[2]*count[1]))< count[2])||(i > count[2]*count[1]*(count[0]-1) && i < count[2]*count[1]*count[0])) {throw exc;}
                else {
                    if (data[i+addy+addx] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + addy + subx > totalpoints || i+addy+subx < 0) {throw exc;}
                if ((i%(count[2]*(count[1]-1)+(numplane*count[2]*count[1]))< count[2])||(i >= 0 && i < count[2]*count[1])) {throw exc;}
                else {
                    if (data[i+addy+subx] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + suby + addx > totalpoints || i+suby+addx < 0) {throw exc;}
                if ((i%(count[2]*count[1]) < count[2])||(i > count[2]*count[1]*(count[0]-1) && i < count[2]*count[1]*count[0])) {throw exc;}
                else {
                    if (data[i+suby+addx] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}
            try {
                //if (i + suby + subx > totalpoints || i+suby+subx < 0) {throw exc;}
                if ((i%(count[2]*count[1]) < count[2])||(i >= 0 && i < count[2]*count[1])) {throw exc;}
                else {
                    if (data[i+suby+subx] < low) {
                        SASA.data.push_back(1); continue;
                    }
                }
            }
            catch (int exc) {}


            SASA.data.push_back(0);
            continue;

        }
        else {
            SASA.data.push_back(0);
            continue;
        }
    }


    return SASA;
}

dx dx::clearByOne(double cutoff, char flag) {
    //makes a copy of the current dx file but eliminates all values below or above a cutoff
    //flag here is designed to determine greater than or less than

    dx ONE;

    ONE.totalpoints = totalpoints;
    for (int i = 0; i < 3; i++) {
        ONE.count[i] = count[i];
        ONE.origin[i] = origin[i];
        ONE.delta[i] = delta[i];
    }

    if (flag == 'g') {
        //set to greater than
        cout << "check greater than\n";
        for (int i = 0; i < totalpoints; i++) {
            if (data[i] < cutoff) {
                ONE.data.push_back(0);
            }
            else {
                ONE.data.push_back(1);
            }
        }
    }
    else if (flag == 'l') {
        //set to less than
        for (int i = 0; i < totalpoints; i++) {
            if (data[i] > cutoff) {
                ONE.data.push_back(0);
            }
            else {
                ONE.data.push_back(1);
            }
        }
    }
    else {
        cerr << "\nERROR!:\n\n"
        << "inequality flag is set improperly\n";
        exit(0);
    }
    return ONE;
}

dx dx::clearByTwo(dx E, double cutoff1, double  cutoff2, char flag1, char flag2) {
    //makes a copy of the current dx file and sets all values to 0 if below or above a cutoff of this and another file
    //cutoff1 applies to current file, cutoff2 to the secondary file
    //flag1 determines greater than or less than for current file
    //flage 2 determines greater than or less than for secondary file
    dx TWO;

    bool check = checkSame(E);
    if (check == false) {
        cout << "The two dx files used for clearing by cutoff do not have the same headers, therefore spitting out an empty dx copy\n";
        return TWO;
    }

    TWO.totalpoints = totalpoints;
    for (int i = 0; i < 3; i++) {
        TWO.count[i] = count[i];
        TWO.origin[i] = origin[i];
        TWO.delta[i] = delta[i];
    }

    if (flag1 == 'g') {
        //looking for points greater than in the current file
        if (flag2 == 'g') {
            //looking for points greater than in the secondary file too
            for (int i = 0; i < totalpoints; i++) {
                if (data[i] < cutoff1) {
                    TWO.data.push_back(0);
                }
                else {
                    if (E.data[i] < cutoff2) {
                        TWO.data.push_back(0);
                    }
                    else {
                        TWO.data.push_back(1);
                    }
                }
            }
        }
        else if (flag2 == 'l') {
            for (int i = 0; i < totalpoints; i++) {
                if (data[i] < cutoff1) {
                    TWO.data.push_back(0);
                }
                else {
                    if (E.data[i] > cutoff2) {
                        TWO.data.push_back(0);
                    }
                    else {
                        TWO.data.push_back(1);
                    }
                }
            }
        }
        else {
            cerr << "\nERROR!:\n\n"
            << "inequality flag2 is set improperly\n";
            exit(0);
        }
    }
    else if (flag1 == 'l') {
        if (flag2 == 'g') {
            //looking for points greater than in the secondary file too
            for (int i = 0; i < totalpoints; i++) {
                if (data[i] > cutoff1) {
                    TWO.data.push_back(0);
                }
                else {
                    if (E.data[i] < cutoff2) {
                        TWO.data.push_back(0);
                    }
                    else {
                        TWO.data.push_back(1);
                    }
                }
            }
        }
        else if (flag2 == 'l') {
            for (int i = 0; i < totalpoints; i++) {
                if (data[i] > cutoff1) {
                    TWO.data.push_back(0);
                }
                else {
                    if (E.data[i] > cutoff2) {
                        TWO.data.push_back(0);
                    }
                    else {
                        TWO.data.push_back(1);
                    }
                }
            }
        }
        else {
            cerr << "\nERROR!:\n\n"
            << "inequality flag2 is set improperly\n";
            exit(0);
        }
    }
    else {
        cerr << "\nERROR!:\n\n"
        << "inequality flag1 is set improperly\n";
        exit(0);
    }
    return TWO;
}


void dx::contour(dx G, double cutoff, char flag) {
    //this method contours current dx file by values in the second, ie if energy isn't above x then g set to 0
    bool check = checkSame(G);
    if (check == false) {
        cout << "The two dx files used for creating contour do not have the same headers, therefore not applying contour\n";
    }
    else {
        if (flag == 'g') {
            for (int i = 0; i < totalpoints; i++) {
                if (G.data[i] < cutoff) {
                    data[i] = 0;
                }
            }
        }
        else if (flag == 'l') {
            for (int i = 0; i < totalpoints; i++) {
                if (G.data[i] > cutoff) {
                    data[i] = 0;
                }
            }
        }
        else {
            cerr << "\nERROR!:\n\n"
            << "inequality flag is set improperly\n";
            exit(0);
        }
    }
}

//int dx::checkNeighbor(int index, vector<int > &grpNum, vector<int > &grpCount, int crnum) {
int dx::checkNeighbor(int index, vector<int > &grpNum, int crnum) {
    int addz = 1;
    int subz = -1*1;
    int addy = count[2];
    int suby = -1*count[2];
    int addx = count[2]*count[1];
    int subx = -1*count[2]*count[1];
    int exc = 0; int dummy;
    int numplane = index/(count[2]*count[1]);
    numgrpatt++;
    if (numgrpatt > 10000) {return 0;}

    try {
        //if (index + addz > totalpoints || index+addz < 0) {throw exc;}
        if (index%count[2] == count[2]-1) {throw exc;}
        else {
            if (data[index+addz] != 0 && grpNum[index+addz] == 0) {
                grpNum[index+addz] =crnum; dummy = checkNeighbor(index+addz, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index + addy > totalpoints || index+addy < 0) {throw exc;}
        if (index%(count[2]*(count[1]-1)+(numplane*count[2]*count[1]))< count[2]) {throw exc;}
        else {
            if (data[index+addy] != 0 && grpNum[index+addy] == 0) {
                grpNum[index+addy] =crnum; dummy = checkNeighbor(index+addy, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index + addx > totalpoints || index+addx < 0) {throw exc;}
        if (index >= count[2]*count[1]*(count[0]-1) && index < count[2]*count[1]*count[0]) {throw exc;}
        else {
            if (data[index+addx] != 0 && grpNum[index+addx] == 0) {
                grpNum[index+addx] =crnum; dummy = checkNeighbor(index+addx, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index + subz > totalpoints || index+subz < 0) {throw exc;}
        if (index%count[2] == 0) {throw exc;}
        else {
            if (data[index+subz] != 0 && grpNum[index+subz] == 0) {
                grpNum[index+subz] =crnum; dummy = checkNeighbor(index+subz, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index + suby > totalpoints || index+suby < 0) {throw exc;}
        if (index%(count[2]*count[1]) < count[2]) {throw exc;}
        else {
            if (data[index+suby] != 0 && grpNum[index+suby] == 0) {
                grpNum[index+suby] =crnum; dummy = checkNeighbor(index+suby, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index + subx > totalpoints || index+subx < 0) {throw exc;}
        if (index >= 0 && index < count[2]*count[1]) {throw exc;}
        else {
            if (data[index+subx] != 0 && grpNum[index+subx] == 0) {
                grpNum[index+subx] =crnum; dummy = checkNeighbor(index+subx, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index+addz+addy > totalpoints || index+addz+addy < 0) {throw exc;}
        if ((index%count[2] == count[2]-1) || (index%(count[2]*(count[1]-1)+(numplane*count[2]*count[1])) < count[2])) {throw exc;}
        else {
            if (data[index+addz+addy] != 0 && grpNum[index+addz+addy] == 0) {
                grpNum[index+addz+addy] =crnum; dummy = checkNeighbor(index+addz+addy, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index+addz+suby > totalpoints || index+addz+suby < 0) {throw exc;}
        if ((index%count[2] == count[2]-1)||(index%(count[2]*count[1]) < count[2])) {throw exc;}
        else {
            if (data[index+addz+suby] != 0 && grpNum[index+addz+suby] == 0) {
                grpNum[index+addz+suby] =crnum; dummy = checkNeighbor(index+addz+suby, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index+subz+addy > totalpoints || index+subz+addy < 0) {throw exc;}
        if ((index%count[2] == 0)||(index%(count[2]*(count[1]-1)+(numplane*count[2]*count[1]))< count[2])) {throw exc;}
        else {
            if (data[index+subz+addy] != 0 && grpNum[index+subz+addy] == 0) {
                grpNum[index+subz+addy] =crnum; dummy = checkNeighbor(index+subz+addy, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index+subz+suby > totalpoints || index+subz+suby < 0) {throw exc;}
        if ((index%count[2] == 0)||(index%(count[2]*count[1]) < count[2])) {throw exc;}
        else {
            if (data[index+subz+suby] != 0 && grpNum[index+subz+suby] == 0) {
                grpNum[index+subz+suby] =crnum; dummy = checkNeighbor(index+subz+suby, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index+addz+addx > totalpoints || index+addz+addx < 0) {throw exc;}
        if ((index%count[2] == count[2]-1)||(index >= count[2]*count[1]*(count[0]-1) && index < count[2]*count[1]*count[0])) {throw exc;}
        else {
            if (data[index+addz+addx] != 0 && grpNum[index+addz+addx] == 0) {
                grpNum[index+addz+addx] =crnum; dummy = checkNeighbor(index+addz+addx, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index+addz+subx > totalpoints || index+addz+subx < 0) {throw exc;}
        if ((index%count[2] == count[2]-1)||(index >= 0 && index < count[2]*count[1])) {throw exc;}
        else {
            if (data[index+addz+subx] != 0 && grpNum[index+addz+subx] == 0) {
                grpNum[index+addz+subx] =crnum; dummy = checkNeighbor(index+addz+subx, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index+subz+addx > totalpoints || index+subz+addx < 0) {throw exc;}
        if ((index%count[2] == 0)||(index >= count[2]*count[1]*(count[0]-1) && index < count[2]*count[1]*count[0])) {throw exc;}
        else {
            if (data[index+subz+addx] != 0 && grpNum[index+subz+addx] == 0) {
                grpNum[index+subz+addx] =crnum; dummy = checkNeighbor(index+subz+addx, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index+subz+subx > totalpoints || index+subz+subx < 0) {throw exc;}
        if ((index%count[2] == 0)||(index >= 0 && index < count[2]*count[1])) {throw exc;}
        else {
            if (data[index+subz+subx] != 0 && grpNum[index+subz+subx] == 0) {
                grpNum[index+subz+subx] =crnum; dummy = checkNeighbor(index+subz+subx, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index+addx+addy > totalpoints || index+addx+addy < 0) {throw exc;}
        if ((index%(count[2]*(count[1]-1)+(numplane*count[2]*count[1]))< count[2])||(index >= count[2]*count[1]*(count[0]-1) && index < count[2]*count[1]*count[0])) {throw exc;}
        else {
            if (data[index+addx+addy] != 0 && grpNum[index+addx+addy] == 0) {
                grpNum[index+addx+addy] =crnum; dummy = checkNeighbor(index+addx+addy, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index+addx+suby > totalpoints || index+addx+suby < 0) {throw exc;}
        if ((index%(count[2]*count[1]) < count[2])||(index >= count[2]*count[1]*(count[0]-1) && index < count[2]*count[1]*count[0])) {throw exc;}
        else {
            if (data[index+addx+suby] != 0 && grpNum[index+addx+suby] == 0) {
                grpNum[index+addx+suby] =crnum; dummy = checkNeighbor(index+addx+suby, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index+subx+addy > totalpoints || index+subx+addy < 0) {throw exc;}
        if ((index%(count[2]*(count[1]-1)+(numplane*count[2]*count[1]))< count[2])||(index >= 0 && index < count[2]*count[1])) {throw exc;}
        else {
            if (data[index+subx+addy] != 0 && grpNum[index+subx+addy] == 0) {
                grpNum[index+subx+addy] =crnum; dummy = checkNeighbor(index+subx+addy, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}
    try {
        //if (index+subx+suby > totalpoints || index+subx+suby < 0) {throw exc;}
        if ((index%(count[2]*count[1]) < count[2])||(index >= 0 && index < count[2]*count[1])) {throw exc;}
        else {
            if (data[index+subx+suby] != 0 && grpNum[index+subx+suby] == 0) {
                grpNum[index+subx+suby] =crnum; dummy = checkNeighbor(index+subx+suby, grpNum, crnum);
            }
        }
    }
    catch (int exc) {}

    return 0;

}

void dx::makeGroups() {
    /*
        This code should only ever be used on an already cutoff dx file, in other words we're grouping non zeros together. Anything else is not implemented
    */
    vector<int > groupNum(totalpoints, 0); //make a vector of length totalpoints initialized to 0
    //vector<int > groupCount; //vector to contain counts in each group index 0 = group 1
    int currnum = 1;
    int num;
    //groupCount.push_back(0);
    for (int i = 0; i < totalpoints; i++) {
        if (data[i] != 0) {
            if (groupNum[i] == 0) {
                groupNum[i] = currnum;
                //num = checkNeighbor(i, groupNum, groupCount, currnum); //this should be a recursive call which will get all neighbors
                num = checkNeighbor(i, groupNum, currnum); //this should be a recursive call which will get all neighbors
                currnum++;
                numgrpatt = 0;
                //groupCount[i]++;
                //groupCount.push_back(0);
            }
        }
    }



    /*
        Begin code for writing the groups. Here we want to:
        print file with all areas which are in any group (master group.dx?)
        print files with each individual group (incrementing file names (group1.dx, group2.dx...groupn.dx)
    */



    char fileName[80];
    sprintf(fileName, "mastergroup.dx");



    ofstream output(fileName);
    output << "object 1 class gridpositions counts " << count[0] << " " << count[1] << " " << count[2] << "\n";
    output << "origin " << origin[0] << " " << origin[1] << " " << origin[2] << "\n";
    output << "delta " << delta[0] << " 0 0\n";
    output << "delta " << "0 " << delta[1] << " 0\n";
    output << "delta " << "0 0 " << delta[2] << "\n";
    output << "object 2 class gridconnections counts " << count[0] << " " << count[1] << " " << count[2] << "\n";
    output << "object 3 class array type float rank 0 items " << totalpoints << " data follows\n";
    int a = 0;
    output.precision(16);
    for (int i = 0; i < totalpoints; i++) {
        if (groupNum[i] != 0) {
            if (a < 2) {
                output << data[i] << " ";
                a++;
            }
            else {
                output << data[i] << "\n";
                a = 0;
            }
        }
        else {
            if (a < 2) {
                output << 0 << " ";
                a++;
            }
            else {
                output << 0 << "\n";
                a = 0;
            }
        }
    }
    output.close();

    sprintf(fileName, "mastergroup.pdb");
    double vals[3];
    int pos = 0;
    FILE * pFile;
    pFile = fopen(fileName, "w");
    //Define the pdb file stuff
    string name = "ATOM"; string atom = "H"; string resname = "T3P"; string chainid = "C"; int resseq = 1; double occupancy = 0.0; double T = 0.0;
    for (int i = 0; i < count[0]; i++) {
        for (int j = 0; j < count[1]; j++) {
            for (int k = 0; k < count[2]; k++) {
                if (groupNum[pos] != 0) {
                    vals[0] = i*delta[0] + origin[0];
                    vals[1] = j*delta[1] + origin[1];
                    vals[2] = k*delta[2] + origin[2];
                    fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), pos/1000, atom.c_str(), resname.c_str(), chainid.c_str(), groupNum[pos], vals[0], vals[1], vals[2], data[pos], T);
                }
                pos++;
            }
        }
    }

    fclose(pFile);

    vector<int > grcount(currnum, 0);

    int val = 0;
    for (int i = 0; i< currnum; i++) {
        val = i+1;
        sprintf(fileName, "group%i.dx", val);
        ofstream output(fileName);
        output << "object 1 class gridpositions counts " << count[0] << " " << count[1] << " " << count[2] << "\n";
        output << "origin " << origin[0] << " " << origin[1] << " " << origin[2] << "\n";
        output << "delta " << delta[0] << " 0 0\n";
        output << "delta " << "0 " << delta[1] << " 0\n";
        output << "delta " << "0 0 " << delta[2] << "\n";
        output << "object 2 class gridconnections counts " << count[0] << " " << count[1] << " " << count[2] << "\n";
        output << "object 3 class array type float rank 0 items " << totalpoints << " data follows\n";
        int a = 0;
        output.precision(16);
        for (int i = 0; i < totalpoints; i++) {
            if (groupNum[i] == val) {
                grcount[val-1]++;
                if (a < 2) {
                    output << data[i] << " ";
                    a++;
                }
                else {
                    output << data[i] << "\n";
                    a = 0;
                }
            }
            else {
                if (a < 2) {
                    output << 0 << " ";
                    a++;
                }
                else {
                    output << 0 << "\n";
                    a = 0;
                }
            }
        }
        output.close();
    }

    ofstream CNT("grcount.txt");
    for (int i = 0; i < currnum; i++) {
        CNT << "group " << i+1 << ": " <<grcount[i] << endl;
    }

}

void dx::mult(dx G) {
    bool check = checkSame(G);
    if (check == false) {
        cout << "The two dx files used for multiplication do not have the same headers, therefore not applying contour\n";
    }
    else {
        for (int i = 0; i < totalpoints; i++) {
            data[i] = data[i]*G.data[i];
        }
    }
}

/*void dx::writeGroups() {
    //to be done later?
}*/


double dx::sum() {
    double sumval = 0;
    for (int i = 0; i < totalpoints; i++) {
        sumval += data[i];
    }
    return sumval;
}

void dx::printPdb(string outfile) {
    double vals[3];
    int pos = 0;
    FILE * pFile;
    pFile = fopen(outfile.c_str(), "w");
    //Define the pdb file stuff
    string name = "ATOM"; string atom = "H"; string resname = "T3P"; string chainid = "C"; int resseq = 1; double occupancy = 0.0; double T = 0.0;
    for (int i = 0; i < count[0]; i++) {
        for (int j = 0; j < count[1]; j++) {
            for (int k = 0; k < count[2]; k++) {
                vals[0] = i*delta[0] + origin[0];
                vals[1] = j*delta[1] + origin[1];
                vals[2] = k*delta[2] + origin[2];
                fprintf (pFile, "%-6s%5i %-4s %3s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f\n", name.c_str(), pos, atom.c_str(), resname.c_str(), chainid.c_str(), resseq, vals[0], vals[1], vals[2], data[pos], T);
                //pos++;
            }
        }
    }

}

void dx::add(dx G) {
    bool check = checkSame(G);
    if (check == false) {
        cout << "The two dx files used for addition do not have the same headers, therefore not applying contour\n";
    }
    else {
        for (int i = 0; i < totalpoints; i++) {
            data[i] = data[i]+G.data[i];
        }
    }
}


void dx::sub(dx G) {
    bool check = checkSame(G);
    if (check == false) {
        cout << "The two dx files used for subtraction do not have the same headers, therefore not applying contour\n";
    }
    else {
        for (int i = 0; i < totalpoints; i++) {
            data[i] = data[i]-G.data[i];
        }
    }
}

void dx::div(dx G) {
    bool check = checkSame(G);
    if (check == false) {
        cout << "The two dx files used for division do not have the same headers, therefore not applying contour\n";
    }
    else {
        for (int i = 0; i < totalpoints; i++) {
            data[i] = data[i]/G.data[i];
        }
    }
}


void dx::addbyConst(double C) {
    for (int i = 0; i < totalpoints; i++) {
        data[i] += C;
    }
}

void dx::multbyConst(double C) {
    for (int i = 0; i < totalpoints; i++) {
        data[i] *= C;
    }
}


double dx::avg() {
    double sumval = 0;
    for (int i = 0; i < totalpoints; i++) {
        sumval += data[i];
    }
    return sumval/totalpoints;
}


pop::pop() {
    tempx = 0.0;
    average = 0.0;
    avgtnorm = 0.0;
    avgonorm = 0.0;
    avgEwwnorm = 0.0;
    avgEswnorm = 0.0;

}

void pop::calcpop(string gistfile) {
    ifstream input(gistfile.c_str());
    if (!input.is_open()) {
	cerr << "Could not find file: " << gistfile << " please check that it is available.\n";
	exit(0);
    }
    getline(input, temp); getline(input, temp);
    while (!input.eof()) {
        input >> tempx >> tempx >> tempx >> tempx >> tempx;
        /*if (tempx != 0) {
            maxx = max(tempx, maxx); minn = min(tempx, minn); pop.push_back(tempx);
        }*/
        N.push_back(tempx);
        input >> tempx >> tempx >> tempx >> tempx;
        transnorm.push_back(tempx);
        input >> tempx >> tempx;
        orientnorm.push_back(tempx);
        input >> tempx >> tempx;
        Eswnorm.push_back(tempx);
        input >> tempx >> tempx;
        Ewwnorm.push_back(tempx);
        input >> tempx >> tempx >> tempx >> tempx >> tempx >> tempx >> tempx;
    }

    for (int i = 0; i < N.size(); i++) {
        average += N[i];
        avgtnorm += transnorm[i];
        avgonorm += orientnorm[i];
        avgEswnorm += Eswnorm[i];
        avgEwwnorm += Ewwnorm[i];
    }

    average /= N.size(); avgtnorm /= N.size(); avgonorm /= N.size(); avgEwwnorm /= N.size(); avgEswnorm /= N.size();
    cout << "average pop: " << average << endl;
    cout << "average tnorm: " << avgtnorm << endl;
    cout << "average onorm: " << avgonorm << endl;
    cout << "average Ewwnorm: " << avgEwwnorm << endl;
    cout << "average Eswnorm: " << avgEswnorm << endl;
}

void dx::write_out_dx(string infile, int column) {
    /*
        This operation will take user specified gist text output file, a dx output from the same system, and the specified column of interest
        It will then print that columns data as a dx file
    */

    int C = column; string temp; double tempx;
    ifstream input(infile.c_str());
    if (!input.is_open()) {
   	cerr << "Could not open file: " << infile << " please ensure it is available.\n";
	exit(0);
    }
    getline(input, temp); //skip the header line
    getline(input, temp); 
    //there are 22 slots per line
	//false this approach is shyte. What happens is any time someone uses a different version or a change occurs this code eats dirt.
	// Instead we will implement a vector, yeah?

    cout << temp << endl;

    //double row[22];
    //vector <double > row;
    vector <double > vals;
    
	//double tempx;

    /*
    while (!input.eof()) {
        for (int i = 0; i < 22; i++) {
            input >> tempx;
            row[i] = tempx; //read in entire row into vector
            //cout << tempx << endl;
        }
        vals.push_back(row[C]); //save column of interest
        //cout << row[C] << endl;
    }

    for (int i = 0; i < totalpoints; i++) {
        data[i] = vals[i]; //reset original dx file stored values with the values from the column read in
        //cout << vals[i] << endl;
    }
    */
	//cout << "Starting to read the input file: " << infile << endl;
	while (!input.eof()) {
		input >> tempx;
		vals.push_back(tempx);
	}
	//cout << "Finished reading the input file: " << infile << endl;
	int rowcount = vals.size()/totalpoints;
	cout << "rowcount: " << rowcount << endl;
	int currrow = -1;	
	int j = 0;
	for (int i = 0; i < vals.size(); i++) {
		if (i%rowcount == C) {
			data[j] = vals[i];
			j++;
		}
	}
}


void dx::calcVdw(double vdw) {
    /*
        This operation will take the dxfile provided and the vdw radius provided and flag all voxels within that radial distance of any other flagged voxel
    */

    double minn[3], maxx[3];
    double x, y, z;

    for (int i = 0; i < 3; i++) {
        minn[i] = 10000.0;
        maxx[i] = 0.0;
    }

    for(int i = 0; i < count[0]; i++) {
        for (int j = 0; j < count[1]; j++) {
            for(int k = 0; k < count[2]; k++) {
                if (data[i+j+k] == 1) {
                    x = i*delta[0] + origin[0];
                    y = j*delta[1] + origin[1];
                    z = k*delta[2] + origin[2];
                    minn[0] = min(minn[0], x); maxx[0] = max(maxx[0], x);
                    minn[1] = min(minn[1], y); maxx[1] = max(maxx[1], y);
                    minn[2] = min(minn[2], z); maxx[2] = max(maxx[2], z);
                }
            }
        }
    }

    double dx, dn; double dt = vdw*vdw;

    for(int i = 0; i < count[0]; i++) {
        for (int j = 0; j < count[1]; j++) {
            for(int k = 0; k < count[2]; k++) {

                x = i*delta[0] + origin[0];
                y = j*delta[1] + origin[1];
                z = k*delta[2] + origin[2];

                dx = pow(maxx[0]-x, 2) + pow(maxx[1] - y, 2) + pow(maxx[2] - z, 2);
                dn = pow(minn[0]-x, 2) + pow(minn[1] - y, 2) + pow(minn[2]-z, 2);

                if (dn < dt || dx < dt) {
                    data[i+j+k] = 1; //set every voxel within our specified region to 1
                }
                else {
                    data[i+j+k] = 0;
                }
            }
        }
    }
}

void dx::histogram(string outfile) {
	double avg = 0;
	double total = 0;
	double stddev = 0;
	double stot = 0;
	double maxx = -10000000;
	double minn = 100000;
	int temp = 0;
	double h = 0; //bin size
	int numbin = 0;
	int adjpts = totalpoints;
	for (int i = 0; i < totalpoints; i++) {
		if (data[i] == 0) {
			adjpts--;			
			continue;
		}		
		else {
			total += data[i];
			maxx = max(maxx, data[i]);
			minn = min(minn, data[i]);
		}		
	}
	avg = total/totalpoints;
	cout << "totalpoints: " << totalpoints << endl;
	cout << "made average and min/max\n";
	cout << "avg: " << avg <<endl;
	cout << "min: " << minn << endl;
	cout << "max: " << maxx << endl;
	for (int i = 0; i < totalpoints; i++) {
		if (data[i] == 0) {
			continue;
		}
		else {
			stot += pow((data[i]-avg), 2);
		}
	}
	stddev = sqrt(stot/(totalpoints-1));
	cout << "made stddev: " << stddev << endl;
	h = (3.5*stddev)/(pow(totalpoints,1./3.));
	cout << "made binsize: " << h << endl;
	numbin = (maxx-minn)/h;
	cout << "made binnum: " << numbin <<endl;
	vector<int> colh;
	colh.resize(numbin);
		

	for (int i = 0; i < totalpoints; i++) {
		if (data[i] == 0) {
			continue;
		}
		else {
			temp = (data[i]-minn)/h;
			temp = floor(temp);
			//cout << temp << endl;
			colh[temp]++;
		}
	}

	ofstream output(outfile.c_str());
	output << "#Average: " << avg << " Stddev: " << stddev << " Bin Width: " << h << " Number bins: " << numbin << endl;

	output << "\n\n\n\n";
	
	for (int i = 0; i < numbin; i++) {
		output << minn+(i*h) << "\t\t" << colh[i] << endl;
	}
	
}

void dx::printcol(string outfile) {
	ofstream output(outfile.c_str());	
	for (int i = 0; i < totalpoints; i++) {
		output << data[i] << endl;
	}
}

dx dx::cat(dx N) {
	/*
		This function will take two dx files that overlap and concatenate them into one new large dx file. This requires that both dx grids line up appropriately, we will not be handling cases where the two .dx files do not line up perfectly
		Generally speaking the incorporated algorithm will determine what edge the two .dx files share in common, then determine the overlap, then finally create a new .dx file that will contain both data sets. Overlapping voxels should be identical, therefore they will be set from the first file only. 
	*/

	

}

dx dx::zeros() {
	dx zero;
	zero.totalpoints = totalpoints;
	for (int i = 0; i < 3; i++ ) {
		zero.count[i] = count[i];
		zero.delta[i] = delta[i];
		zero.origin[i] = origin[i];
	}
	for (int i = 0; i < totalpoints; i++) {
		zero.data.push_back(0.0);
	}
	return zero;
}

void dx::citation() {
	cout << "***********************************************************************************************************************************" << endl;
	cout << "***********************************************************************************************************************************" << endl;
	cout << "If data created through the use of GISTPP is used in a publication please cite the following source: " << endl;
	cout << " Ramsey, S., Nguyen, C., Salomon‐Ferrer, R., Walker, R. C., Gilson, M. K., & " << endl;
	cout << " Kurtzman, T. (2016). Solvation thermodynamic mapping of molecular " << endl;
	cout << " surfaces in AmberTools: GIST. Journal of computational chemistry, 37(21), " << endl;
	cout << " 2029-2037. " << endl;	
	cout << "***********************************************************************************************************************************" << endl;
	cout << "***********************************************************************************************************************************" << endl;	
}
