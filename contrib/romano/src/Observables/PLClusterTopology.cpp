//
// Created by josh on 8/14/25.
//

/*
 * PLClusterTopology.cpp
 *
 *  Created on: Feb 12, 2013
 *      Author: petr
 */

#include "PLClusterTopology.h"
#include "../Interactions/PatchyShapeInteraction.h"

#include <sstream>
#include <map>
#include <iostream>
#include <vector>


using namespace std;
//
//// Data structure to store a graph edge
//struct Edge {
//    int src, dest;
//};
//
//
//// A class to represent a graph object
//class Graph
//{
//public:
//    // a vector of vectors to represent an adjacency list
//    vector<vector<int>> adjList;
//
//    // Graph Constructor
//    Graph(vector<Edge> const &edges, int n)
//    {
//        // resize the vector to hold `n` elements of type `vector<int>`
//        adjList.resize(n);
//
//        // add edges to the directed graph
//        for (auto &edge: edges)
//        {
//            // insert at the end
//            adjList[edge.src].push_back(edge.dest);
//
//            // uncomment the following code for undirected graph
//            // adjList[edge.dest].push_back(edge.src);
//        }
//    }
//};
//
//// Function to print adjacency list representation of a graph
//void printGraph(Graph const &graph, int n)
//{
//    for (int i = 0; i < n; i++)
//    {
//        // print the current vertex number
//        cout << i << " ——> ";
//
//        // print all neighboring vertices of a vertex `i`
//        for (int v: graph.adjList[i]) {
//            cout << v << " ";
//        }
//        cout << endl;
//    }
//}


PLClusterTopology::PLClusterTopology() {
    _show_types = true;
}

PLClusterTopology::~PLClusterTopology() = default;

void PLClusterTopology::get_settings(input_file &my_inp, input_file &sim_inp) {
    bool show_types = false;

    if( getInputBool(&my_inp,"show_types",&show_types,1) == KEY_FOUND)
    {
        this->_show_types = show_types;
    }

}

std::string PLClusterTopology::get_output_string(llint curr_step) {

    PatchyShapeInteraction *interaction = dynamic_cast<PatchyShapeInteraction * >(CONFIG_INFO->interaction);
    interaction->check_patchy_locks();

    BaseParticle *p;
    int N = CONFIG_INFO->N();
    int *cluster_membership = new int[N];
    for (int i = 0; i < N; i++)
        cluster_membership[i] = -1;
    //BaseParticle *q;
    number pair_energy;
    //number total_energy;
    int cluster_index = 0;
    std::stringstream output_str;
    int cluster_count = 0;


    std::vector<std::vector<int>> adjList;

    for (int i = 0; i < CONFIG_INFO->N(); i++) {
        p = CONFIG_INFO->particles()[i];
        std::vector<BaseParticle *> neighs =
                CONFIG_INFO->lists->get_all_neighbours(p);

        vector<int> i_neighbors;

        //printf("Particle %d has %d neighbors \n",i,neighs.size());
        for (unsigned int j = 0; j < neighs.size(); j++) {
            pair_energy = CONFIG_INFO->interaction->pair_interaction_term(
                    PatchyShapeInteraction::PATCHY, p, neighs[j]
            );
            if (pair_energy < 0 && i < neighs[j]->index) {
                i_neighbors.push_back(neighs[j]->index);
                if (
                        cluster_membership[p->index] == -1
                        && cluster_membership[neighs[j]->index] == -1
                        ) {
                    cluster_membership[p->index] =
                    cluster_membership[neighs[j]->index] =
                            cluster_index;
                    cluster_index++;
                } else if (
                        cluster_membership[p->index] == -1
                        && cluster_membership[neighs[j]->index] > -1
                        ) {
                    cluster_membership[p->index] =
                            cluster_membership[neighs[j]->index];
                } else if (
                        cluster_membership[p->index] > -1
                        && cluster_membership[neighs[j]->index] == -1
                        ) {
                    cluster_membership[neighs[j]->index] =
                            cluster_membership[p->index];
                } else if (
                        cluster_membership[p->index] > -1
                        && cluster_membership[neighs[j]->index] > -1
                        && cluster_membership[neighs[j]->index]
                           != cluster_membership[p->index]
                        ) {
                    int oldindex = cluster_membership[neighs[j]->index];
                    cluster_membership[neighs[j]->index] =
                            cluster_membership[p->index];
                    for (int k = 0; k < N; k++) {
                        if (cluster_membership[k] == oldindex) {
                            cluster_membership[k] =
                                    cluster_membership[p->index];
                        }
                    }
                }
            }
        }
        adjList.push_back(i_neighbors);
    }

    if (cluster_index > 0) {
        vector<int> cluster_members;
        for (int j = 0; j < cluster_index; j++) {
            cluster_members.clear();
            std::stringstream clusout;
            std::stringstream topologyout;
            clusout << "( ";
            bool isempty = true;
            for (int k = 0; k < N; k++) {
                if (cluster_membership[k] == j) {
                    vector<int> mylist = adjList[k];

                    if (mylist.size() > 0)
                    {
                        topologyout << k << " -> ("  ;
                        for (auto neigh = mylist.begin(); neigh != mylist.end(); ++neigh)
                        {
                            topologyout << *neigh << " ";
                        }
                        // Remove last space
                        topologyout.seekp(-1, topologyout.cur);

                        topologyout << "), ";
                    }

                    if(! this->_show_types)
                    {
                        cluster_members.push_back(k);
                    }
                    else
                    {
                        cluster_members.push_back(CONFIG_INFO->particles()[k]->type );
                    }
                    isempty = false;
                }
            }

            if (!isempty) {
                cluster_count++;
                //std::sort(cluster_members.begin(),cluster_members.end());
                for(vector<int>::iterator k = cluster_members.begin(); k != cluster_members.end(); ++k)
                {
                    clusout << *k << " ";
                }
                // Remove last space
                topologyout.seekp(-1, topologyout.cur);

                clusout << ") ";
                output_str << clusout.str() << " [" << topologyout.str();

                // Remove last comma and space
                output_str.seekp(-2, output_str.cur);

                output_str << "]";
            }
        }
    }


    delete [] cluster_membership;
    std::stringstream out;
    out << cluster_count << " " << output_str.str();
    return out.str();
}

