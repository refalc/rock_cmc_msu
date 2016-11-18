#ifndef _ROCK_H_
#define _ROCK_H_

#include <vector>
#include <list>

//#include "container/adjacency_matrix.hpp"

//#include "cluster/cluster_algorithm.hpp"

//#include "definitions.hpp"


#include <memory>
#include <cstdlib>
#include <stdexcept>

#include <string>
#include <fstream>
#include <sstream>

#include <stack>
#include <cmath>
#include <algorithm>

 

typedef std::vector<double> pattern;

typedef std::vector<double> point;

typedef std::vector<point> dataset;


 
typedef struct data_representation {
public:
    unsigned int            size;
    unsigned int            dimension;
    double                  ** objects;
} data_representation;



namespace container {

 
class adjacency_collection {
public:
     
    virtual ~adjacency_collection(void) { }

public:
     
    virtual size_t size(void) const = 0;

     
    virtual void set_connection(const size_t node_index1, const size_t node_index2) = 0;

     
    virtual void erase_connection(const size_t node_index1, const size_t node_index2) = 0;

     
    virtual bool has_connection(const size_t node_index1, const size_t node_index2) const = 0;

     
    virtual void get_neighbors(const size_t node_index, std::vector<size_t> & node_neighbors) const = 0;

     
    virtual void clear(void) = 0;
};



 
class adjacency_weight_collection : public adjacency_collection {
public:
     
    virtual void set_connection_weight(const size_t node_index1, const size_t node_index2, const double weight) = 0;

     
    virtual double get_connection_weight(const size_t node_index1, const size_t node_index2) const = 0;
};

 
class adjacency_matrix : public adjacency_weight_collection {
private:
    typedef std::vector<std::vector<double> >  adjacency_matrix_container;

protected:
    adjacency_matrix_container  m_adjacency;

public:
     
    adjacency_matrix(void);

     
    adjacency_matrix(const adjacency_matrix & another_matrix);

     
    //adjacency_matrix(adjacency_matrix && another_matrix);

     
    adjacency_matrix(const size_t node_amount);

     
    virtual ~adjacency_matrix(void);


private:
     
    static const double DEFAULT_EXISTANCE_CONNECTION_VALUE;

     
    static const double DEFAULT_NON_EXISTANCE_CONNECTION_VALUE;


public:
     
    virtual size_t size(void) const;

     
    virtual void set_connection(const size_t node_index1, const size_t node_index2);

     
    virtual void erase_connection(const size_t node_index1, const size_t node_index2);

     
    virtual bool has_connection(const size_t node_index1, const size_t node_index2) const;

     
    virtual void get_neighbors(const size_t node_index, std::vector<size_t> & node_neighbors) const;

     
    virtual void set_connection_weight(const size_t node_index1, const size_t node_index2, const double weight);

     
    virtual double get_connection_weight(const size_t node_index1, const size_t node_index2) const;

     
    virtual void clear(void);

public:
    adjacency_matrix & operator=(const adjacency_matrix & another_collection);

    //adjacency_matrix & operator=(adjacency_matrix && another_collection);
};

}

using namespace container;


namespace cluster_analysis {

typedef std::vector<size_t> noise ;
typedef std::vector<size_t> cluster  ;
typedef std::vector<cluster> cluster_sequence  ;
typedef cluster_sequence* cluster_sequence_ptr  ;


 
class cluster_data {
protected:
    cluster_sequence_ptr    m_clusters;

public:
     
    cluster_data(void);

     
    cluster_data(const cluster_data & p_other);

     
    //cluster_data(cluster_data && p_other);

     
    virtual ~cluster_data(void);

public:
     
    cluster_sequence_ptr clusters(void);

     
    void resize_clusters(std::size_t p_count_clusters, std::size_t p_count_data);

     
    size_t size(void) const;

public:
     
    cluster & operator[](const size_t p_index);

     
    const cluster & operator[](const size_t p_index) const;

     
    cluster_data & operator=(const cluster_data & p_other);

     
    //cluster_data & operator=(cluster_data && p_other);

     
    bool operator==(const cluster_data & p_other) const;

     
    bool operator!=(const cluster_data & p_other) const;
};

typedef cluster_data rock_data ;

 
class cluster_algorithm {
public:
     
    virtual ~cluster_algorithm(void);

public:
     
    virtual void process(const dataset & p_data, cluster_data & p_result) = 0;
};


class rock : public cluster_algorithm {
private:
     
    typedef std::vector<cluster> rock_cluster_sequence;

private:
    adjacency_matrix        m_adjacency_matrix;

    double                  m_radius;

    double                  m_degree_normalization;

    size_t                  m_number_clusters;

    rock_cluster_sequence   m_clusters;

public:
     
    rock(void);

     
    rock(const double radius, const size_t number_clusters, const double threshold);

     
    virtual ~rock(void);

public:
     
    virtual void process(const dataset & p_data, cluster_data & p_result);

private:
     
    void create_adjacency_matrix(const dataset & p_data);

     
    bool merge_cluster(void);

     
    size_t calculate_links(const cluster & cluster1, const cluster & cluster2) const;

     
    double calculate_goodness(const cluster & cluster1, const cluster & cluster2) const;
};


}

// next utils and other

long ReadSample(const std::string &file_name, dataset &set);


#endif // ROCK_H
