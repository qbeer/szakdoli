#include <iostream>
#include <vector>
#include <utility> 
#include <set>
#include <array>
#include <string>
#include <fstream>
#include <numeric> // for std::accumulate
#include <sstream> // for std::istringstream
#include <algorithm> // for std::for_each, std::find, etc.
#include <cmath> // for std::sqrt
#include <list>
#include <chrono> // to measure time

/*
 * @Author: Olar Alex
 * @From: GeeksForGeeks
 */

class Graph{

private:

    int numberOfVertices;
    std::list<int>* adjecencyList;
    std::vector<std::vector<int>> BFSroute;

public:

    // constructor of the graph in adjecency list representation
    Graph(int numberOfVertices){
        // this represents a pointer to the current object
        this->numberOfVertices = numberOfVertices;
        this->adjecencyList = new std::list<int> [this->numberOfVertices];
    }

    void addEdge(int from, int to){
        // it is a UNDIRECTED graph 
        adjecencyList[from].push_back(to);
        adjecencyList[to].push_back(from);
    }

    void BFS(){
        std::vector<bool> visitedVertices(this->numberOfVertices, false);
        int startingVertex = 0;
        while( std::count(visitedVertices.begin(), visitedVertices.end(), true) != 
                                                            this->numberOfVertices ){
            std::vector<int> route;
            visitedVertices[startingVertex] = true;
            std::list<int> enqueuedVertices;
            enqueuedVertices.push_back(startingVertex);
            int currentVertex;
            while(!enqueuedVertices.empty()){
                currentVertex = enqueuedVertices.front();
                route.push_back(currentVertex);
                enqueuedVertices.pop_front();

                std::for_each( adjecencyList[currentVertex].begin(),
                adjecencyList[currentVertex].end(),
                [&](const int& currentNeighbor){ if(!visitedVertices[currentNeighbor]){
                    visitedVertices[currentNeighbor] = true;
                    enqueuedVertices.push_back(currentNeighbor);
                }} 
                );
            }
            BFSroute.push_back(route);
            startingVertex = std::find(visitedVertices.begin(), visitedVertices.end(), false) 
                                                                    - visitedVertices.begin();
        }
    }

    void printBFS(){
        if(BFSroute.size()!=0){
            for(auto route : BFSroute){
                if(route.size() > 1){
                for(auto vertex : route){
                    std::cout << vertex << " ";
                }
                std::cout << std::endl;
                }
            }
            std::cout << "Number of clusters: " << std::count_if(BFSroute.begin(),
            BFSroute.end(), [&](std::vector<int>& route){ return route.size() > 1; })
            << std::endl;
        }else{
            std::cout << "BFS algorithm has not been executed yet." << std::endl;
        }
    }

    std::vector< std::vector<int> > getRoutes(){
        std::vector<std::vector<int>> routes;
        if(BFSroute.size()!=0){
            for(auto route : BFSroute){
                if(route.size() > 1){
                    routes.push_back(route);
                }
            }
        }
        return routes;
    }

};

double infinity = 10000.; // handled as infinity in prim's algorithm
constexpr int dim = 6;
constexpr double scalingFactor = 0.0000000000000; // for distance calculations


struct DataPoints{
    
        std::vector<int> hadron;
        std::vector<int> charge;
        std::vector<double> mass;
        std::vector< std::array< double, dim > > posAndMom;
        std::vector<int> num_of_coll;
    
};

std::istream& operator>>(std::istream& is, DataPoints& points);

double distance( const std::array<double, dim>& a, const std::array<double, dim>& b);

int main(int argc, char* argv[]){

    std::ifstream fin(argv[1]);
    infinity = atof(argv[2]);

    DataPoints data;
    fin >> data;

    int numberOfVertices = data.mass.size();

    Graph* g = new Graph(numberOfVertices);

    int P,Q;
    double distanceBetweenVertices;

    auto t0 = std::chrono::high_resolution_clock::now();
    
    for( unsigned int i = 0; i < numberOfVertices; i++ ){
    
        P = i;
    
        for( unsigned int j = 0; j < i; j++ ){
    
            Q = j;

            distanceBetweenVertices = distance( data.posAndMom[i],data.posAndMom[j] );
    
            if(distanceBetweenVertices < infinity){
    
                g->addEdge(P,Q);
    
            }
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();

 //   std::cout << "Graph construction took: " << 
  //  std::chrono::duration_cast< std::chrono::microseconds >(t1-t0).count() << " microseconds\n";

    auto t2 = std::chrono::high_resolution_clock::now();

    g->BFS();

    auto t3 = std::chrono::high_resolution_clock::now();

//std::cout << "Clustering took: " <<
  //  std::chrono::duration_cast< std::chrono::microseconds >(t3-t2).count() << " microseconds\n";

    //g->printBFS();

    std::vector<std::vector<int>> routes = g->getRoutes();

    for(auto route : routes){
        int numberOfProtons = std::count_if( route.begin(), route.end(), [&](const int& routeId){
            return data.charge[routeId] == 1;
        });
        for(int routeId : route){
            for(int i = 0; i < 3; i++){
                std::cout << data.posAndMom[routeId][i] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "A : " << route.size() << "\t" << "Z : " << numberOfProtons << std::endl;
        std::cout << std::endl;
        int polarAngle = std::
    }

    return 0;

}

std::istream& operator>>(std::istream& is, DataPoints& points){
    
    std::array<double, dim> tempPosAndMom;
        
    std::string line;
    
    while(std::getline(is, line)){
    
        double tempData;
    
        std::istringstream rawData(line);
    
        rawData >> tempData;
        points.hadron.push_back(int(tempData));
        rawData >> tempData;
        points.charge.push_back(int(tempData));
        rawData >> tempData;
        points.mass.push_back(tempData);
    
        for(unsigned int i=0; i < dim; i++){
            rawData >> tempPosAndMom[i];
        }
    
        points.posAndMom.push_back(tempPosAndMom);
    
        rawData >> tempData;
        points.num_of_coll.push_back(int(tempData));
    
    }
    
    return is;
    
}


double distance( const std::array<double, dim>& a, const std::array<double, dim>& b){
    
    double distPos = 0.0;
    double distMom = 0.0;
        
    for(unsigned int i = 0; i < dim/2; i++){
        
        distPos += (a[i] - b[i])*(a[i] - b[i]);
        distMom += scalingFactor*(a[dim/2+i] - b[dim/2+i])*(a[dim/2+i] - b[dim/2+i]);    
         
    }
        
    return std::sqrt(distPos+distMom);
        
}