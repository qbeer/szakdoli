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
#include <math.h>
#include <iomanip>

/*
 * @Author: Olar Alex
 * @From: GeeksForGeeks
 */

constexpr double PI = 3.14159;

template<typename T>
T SQR(T x){
     return x*x;
}

 struct Cluster{

     Cluster(){}

     Cluster(std::array<double, 3> com,
             std::array<double, 3> momentum,
             int charge,
             double numberOfNucleons,
             double velocity,
             double energy,
             double polarAngle){
         this->com = com;
         this->momentum = momentum;
         this->charge = charge;
         this->numberOfNucleons = numberOfNucleons;
         this->velocity = velocity;
         this->energy = energy;
         this->polarAngle = polarAngle;
     }

     std::array<double, 3> com;
     std::array<double, 3> momentum;
     int charge;
     double numberOfNucleons;
     double velocity;
     double energy;
     double polarAngle;

 };

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

std::vector<Cluster> getClusters(const DataPoints& data, const std::vector<std::vector<int>>& routes);

int main(int argc, char* argv[]){

    std::vector<Cluster> allClusters;
    infinity = atof(argv[1]);
    std::string directory = std::string(argv[2]);

    for(int filename=0; filename < atoi(argv[3]); filename++){

        std::ifstream fin( directory + "/" + std::to_string(filename) );

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

        std::vector<Cluster> clusters = getClusters(data, g->getRoutes());

        allClusters.insert(allClusters.end(), clusters.begin(), clusters.end());

    }

    std::cout << std::right <<  std::setw(10) << "p_x" << "\t" << std::setw(10) << "p_y" << "\t" 
        << std::setw(10) << "p_z" << "\t" << std::setw(8) << "v" <<  "\t" << std::setw(8) << 
        "E" << "\t" << "A" << "\t" << "Z"  << "\t" << "polar angle" <<"\n\n";
    
    std::for_each(allClusters.begin(), allClusters.end(), [](Cluster& cluster){
        std::cout << std::right <<  std::setw(10) << cluster.momentum[0] << "\t" << std::setw(10) << cluster.momentum[1] << "\t" 
        << std::setw(10) << cluster.momentum[2] << "\t" << std::setw(8) << cluster.velocity <<  "\t" << std::setw(8) << 
        cluster.energy << "\t" << cluster.numberOfNucleons << "\t" << cluster.charge << "\t" << cluster.polarAngle << std::endl;
    });

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
        
        distPos += scalingFactor*(a[i] - b[i])*(a[i] - b[i]);
        distMom += (a[dim/2+i] - b[dim/2+i])*(a[dim/2+i] - b[dim/2+i]);    
         
    }
        
    return std::sqrt(distPos+distMom);
        
}

std::vector<Cluster> getClusters(const DataPoints& data, const std::vector<std::vector<int>>& routes){

        std::vector<Cluster> clusters; 

        // cluster
        for(const auto& route : routes){
            
            int numberOfProtons = std::count_if( route.begin(), route.end(), [&](const int& routeId){
                return data.charge[routeId] == 1;});

            int numberOfNucleons = route.size();

            std::array<double, 3> com;
            std::array<double, 3> momentum;

            com.fill(0);
            momentum.fill(0);

            // cluster element
            for(const int& id : route){
                for(int i = 0; i < 3; i++){
                    com[i] += data.posAndMom[id][i];
                    momentum[i] += data.posAndMom[id][i+3];
                }
            }

            for( double& elm : com ){
                elm /= (double)numberOfNucleons;
            }

            double momentumLength = std::sqrt(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2]);

            double clusterMass = numberOfProtons * 0.938 + (numberOfNucleons - numberOfProtons)* 0.939;

            double clusterEnergy = std::sqrt(
                SQR(momentumLength) + 
                SQR(clusterMass)
            );

            double clusterVelocity = momentumLength / clusterEnergy;

            double polarAngle = std::atan( std::sqrt(SQR(momentum[1])+SQR(momentum[0])) / momentum[3]) * 180/PI;

            Cluster elem( com, momentum,  numberOfProtons, numberOfNucleons, clusterVelocity,
                            clusterEnergy, polarAngle);

            clusters.push_back(elem);

        }

        return clusters;

};