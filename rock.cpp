#include "rock.hpp"

#include <cmath>
#include <climits>
#include <iostream>
#include <omp.h>

//#include "utils.hpp"

//utils
inline double euclidean_distance_sqrt(const std::vector<double> * const point1, const std::vector<double> * const point2) {
    double distance = 0.0;
    /* assert(point1->size() != point1->size()); */
    for (unsigned int dimension = 0; dimension < point1->size(); dimension++) {
        double difference = (point1->data()[dimension] - point2->data()[dimension]);
        distance += difference * difference;
    }

    return distance;
}

//
namespace container {

const double adjacency_matrix::DEFAULT_EXISTANCE_CONNECTION_VALUE = 1.0;
const double adjacency_matrix::DEFAULT_NON_EXISTANCE_CONNECTION_VALUE = 0.0;


adjacency_matrix::adjacency_matrix(void) {
    m_adjacency = adjacency_matrix_container();
}


adjacency_matrix::adjacency_matrix(const adjacency_matrix & another_matrix) {
    m_adjacency = another_matrix.m_adjacency;
}


/*adjacency_matrix::adjacency_matrix(adjacency_matrix && another_matrix) {
    m_adjacency = std::move(another_matrix.m_adjacency);
}*/


adjacency_matrix::adjacency_matrix(const size_t node_amount) : m_adjacency(node_amount, std::vector<double>(node_amount, DEFAULT_NON_EXISTANCE_CONNECTION_VALUE)) { }


adjacency_matrix::~adjacency_matrix(void) { }


size_t adjacency_matrix::size(void) const { return m_adjacency.size(); }


void adjacency_matrix::set_connection(const size_t node_index1, const size_t node_index2) {
    m_adjacency[node_index1][node_index2] = DEFAULT_EXISTANCE_CONNECTION_VALUE;
}


bool adjacency_matrix::has_connection(const size_t node_index1, const size_t node_index2) const {
    return (m_adjacency[node_index1][node_index2] != DEFAULT_NON_EXISTANCE_CONNECTION_VALUE);
}


void adjacency_matrix::erase_connection(const size_t node_index1, const size_t node_index2) {
    m_adjacency[node_index1][node_index2] = DEFAULT_NON_EXISTANCE_CONNECTION_VALUE;
}


void adjacency_matrix::get_neighbors(const size_t node_index, std::vector<size_t> & node_neighbors) const {
    node_neighbors.clear();

    const std::vector<double> & node_neighbor_connections = m_adjacency[node_index];
    for (size_t neighbor_index = 0; neighbor_index != node_neighbor_connections.size(); neighbor_index++) {
        if (node_neighbor_connections[neighbor_index] != DEFAULT_NON_EXISTANCE_CONNECTION_VALUE) {
            node_neighbors.push_back(neighbor_index);
        }
    }
}


void adjacency_matrix::set_connection_weight(const size_t node_index1, const size_t node_index2, const double weight_connection) {
    m_adjacency[node_index1][node_index2] = weight_connection;
}


double adjacency_matrix::get_connection_weight(const size_t node_index1, const size_t node_index2) const {
    return m_adjacency[node_index1][node_index2];
}


void adjacency_matrix::clear(void) {
    m_adjacency.clear();
}


adjacency_matrix & adjacency_matrix::operator=(const adjacency_matrix & another_collection) {
    if (this != &another_collection) {
        m_adjacency = another_collection.m_adjacency;
    }

    return *this;
}


/*adjacency_matrix & adjacency_matrix::operator=(adjacency_matrix && another_collection) {
    if (this != &another_collection) {
        m_adjacency = std::move(another_collection.m_adjacency);
    }

    return *this;
}*/

}

namespace cluster_analysis {

cluster_algorithm::~cluster_algorithm(void) { }

cluster_data::cluster_data(void) : m_clusters(new cluster_sequence()) { }


cluster_data::cluster_data(const cluster_data & p_other) : m_clusters(p_other.m_clusters) { }


//cluster_data::cluster_data(cluster_data && p_other) : m_clusters(std::move(p_other.m_clusters)) { }


cluster_data::~cluster_data(void) { }


cluster_sequence_ptr cluster_data::clusters(void) { return m_clusters; }


size_t cluster_data::size(void) const { return m_clusters->size(); }


cluster & cluster_data::operator[](const size_t p_index) { return (*m_clusters)[p_index]; }


const cluster & cluster_data::operator[](const size_t p_index) const { return (*m_clusters)[p_index]; }


cluster_data & cluster_data::operator=(const cluster_data & p_other) {
    if (this != &p_other) {
        m_clusters = p_other.m_clusters;
    }

    return *this;
}


/*cluster_data & cluster_data::operator=(cluster_data && p_other) {
    if (this != &p_other) {
        m_clusters = std::move(p_other.m_clusters);
    }

    return *this;
}*/


bool cluster_data::operator==(const cluster_data & p_other) const {
    return (m_clusters == p_other.m_clusters);
}


bool cluster_data::operator!=(const cluster_data & p_other) const {
    return !(*this == p_other);
}

rock::rock(void) :
    m_adjacency_matrix(adjacency_matrix()),
    m_radius(0.0),
    m_degree_normalization(0.0),
    m_number_clusters(0)
{ }


rock::rock(const double radius, const size_t num_clusters, const double threshold) :
    m_adjacency_matrix(adjacency_matrix()),
    m_radius(radius * radius),
    m_degree_normalization(1.0 + 2.0 * ( (1.0 - threshold) / (1.0 + threshold) )),
    m_number_clusters(num_clusters)
{ }


rock::~rock(void) { }


void rock::process(const dataset & p_data, cluster_data & p_result) {
    std::cout << "Create matrix" << std::endl;
    create_adjacency_matrix(p_data);
    std::cout << "Create matrix ok" << std::endl;
    /* initialize first version of clusters */
    for (size_t index = 0; index < p_data.size(); index++) {
        m_clusters.push_back(cluster(1, index));
    }

    while( (m_number_clusters < m_clusters.size()) && (merge_cluster()) ) { }

    /* copy results to the output result (it much more optimal to store in list representation for ROCK algorithm) */
    p_result = rock_data();
    p_result.clusters()->insert(p_result.clusters()->begin(), m_clusters.begin(), m_clusters.end());

    m_clusters.clear();         /* no need it anymore - clear to save memory */
    m_adjacency_matrix.clear(); /* no need it anymore - clear to save memory */
}


void rock::create_adjacency_matrix(const dataset & p_data) {
    m_adjacency_matrix = adjacency_matrix(p_data.size());
    #pragma omp parallel for
    for (size_t i = 0; i < m_adjacency_matrix.size(); i++) {
        for (size_t j = i + 1; j < m_adjacency_matrix.size(); j++) {
            double distance = euclidean_distance_sqrt(&p_data[i], &p_data[j]);
            if (distance < m_radius) {
                m_adjacency_matrix.set_connection(i, j);
                m_adjacency_matrix.set_connection(j, i);
            }
        }
    }
}


bool rock::merge_cluster(void) {
    rock_cluster_sequence::iterator cluster1 = m_clusters.end();
    rock_cluster_sequence::iterator cluster2 = m_clusters.end();
    rock_cluster_sequence::iterator iter_i = m_clusters.begin(), iter_j = m_clusters.begin();
    double maximum_goodness = 0;

    omp_set_num_threads(1);
    #pragma omp parallel for
    for( int i = 0; i < m_clusters.size(); i++ ) {
        for ( int j = i + 1; j < m_clusters.size(); j++ ) {
            double goodness = calculate_goodness(m_clusters[i], m_clusters[j]);
            if (goodness > maximum_goodness) {
                maximum_goodness = goodness;

                cluster1 = m_clusters.begin() + i;
                cluster2 = m_clusters.begin() + j;
            }
        }
    }

    if (cluster1 == cluster2) {
        return false;   /* clusters are totally separated (no links between them), it's impossible to made a desicion which of them should be merged */
    }

    (*cluster1).insert((*cluster1).end(), (*cluster2).begin(), (*cluster2).end());
    m_clusters.erase(cluster2);

    return true;
}

size_t rock::calculate_links(const cluster & cluster1, const cluster & cluster2) const {
size_t number_links = 0;

    std::vector<size_t> links;
    links.resize(cluster1.size());

    for( int i = 0; i < links.size(); i++ ) {
        links[i] = 0;
    }

    //std::cout <<"\nThread "<<omp_get_thread_num()<<" on cpu "<<sched_getcpu()<<std::endl;
    //#pragma omp parallel for
    for( int i = 0; i < cluster1.size(); i++) {
        for( int j = 0; j < cluster2.size(); j++) {
            links[i] += (size_t) m_adjacency_matrix.has_connection(cluster1[i], cluster2[j]);
        }
    }

    for( int i = 0; i < cluster1.size(); i++ ) {
        number_links += links[i];
    }

    return number_links;
}

double rock::calculate_goodness(const cluster & cluster1, const cluster & cluster2) const {
    const double number_links = calculate_links(cluster1, cluster2);

    const double size_cluster1 = (double) cluster1.size();
    const double size_cluster2 = (double) cluster2.size();

    return number_links / ( std::pow( size_cluster1 + size_cluster2, m_degree_normalization ) -
        std::pow( size_cluster1, m_degree_normalization ) -
        std::pow( size_cluster2, m_degree_normalization ) );
}


}

long ReadSample(const std::string &file_name, dataset &set)
{
    std::fstream pFile;
    pFile.open(file_name.c_str());

    if( !pFile.is_open() ) {
        std::cout << "Catton open sample" << std::endl;
        return -1;
    }

    std::string line, coord;
    int space_pos = 0;

    point temp_point;
    while( std::getline(pFile, line) ) {
        space_pos = line.find(" ");
        temp_point.clear();

        // X
        coord = line.substr(0, space_pos);
        temp_point.push_back(std::strtod(coord.c_str(), NULL));

        // Y
        coord = line.substr(space_pos + 1);
        temp_point.push_back(std::strtod(coord.c_str(), NULL));

        set.push_back(temp_point);
    }

    return 0;
}
