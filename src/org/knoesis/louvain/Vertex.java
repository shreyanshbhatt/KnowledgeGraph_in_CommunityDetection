package org.knoesis.louvain;

import java.util.ArrayList;

class EdgeWeights {
  double total_edge_weight;
  ArrayList<Double> feature_edge_weight_derivative_helper;
}


public class Vertex {
  int id;
  double total_vertex_degree;
  ArrayList<Neighbor> neighbors;
  int prev_iter_community_label;
  int community_label;
  boolean should_skip;
  ArrayList<Double> feature_vertex_derivate_helper;
}
