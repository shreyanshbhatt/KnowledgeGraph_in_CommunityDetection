package org.knoesis.louvain;

import java.util.ArrayList;

public class Neighbor {
  int neighbour_label;
  ArrayList<Double> similarity;
  ArrayList<Double> weighted_similarity_der;
  double edge_weight;
}
