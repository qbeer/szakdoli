#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>

int parallerRuns;
int subsequentRuns;

struct DataPoints{
    
        std::vector<int> hadron;
        std::vector<int> charge;
        std::vector<double> mass;
        std::vector< std::array< double, dim > > posAndMom;
        std::vector<int> num_of_coll;
    
};

int main(int argc, char* argv[]){

    DataPoints data;

    // first param is the input file name
    std::ifstream fin( argv[1] );

    parallerRuns = atoi(argv[2]);
    subsequentRuns = atoi(argv[3]);

    fin >> data;

    return 0;

}

std::istream& operator>>(std::istream& is, DataPoints& points){
    
    std::array<double, dim> tempPosAndMom;
        
    std::string line;

    for(int cntr = 0; cntr < subsequentRuns; cntr++){

        std::getline(is, line);
        std::getline(is, line);
    
        for(int inner = 0; inner < parallerRuns; inner++){

            std::getline(is, line);
            std::istringstream rowCount(line);

            int row;

            rowCount >> row;

            std::ofstream fout( std::to_string(inner + counter*inner) );

            for(int data = 0; data < row; data++){
    
                std::getline(is, line);

                line >> fout;

    
            }

            fout.close();

        }

    }
    
    return is;
    
}