package org.knoesis.louvain;

import java.util.ArrayList;
import java.util.HashMap;

public class GradientHelper {
  private static GradientHelper gHelper = null;
  Graph g;

  public static GradientHelper getInstance(int total_comms, int total_features, Graph g) {
    if (gHelper != null || true) { // we are creating new structure every time
                                   // since we do not want to initialize
      gHelper = new GradientHelper(total_comms, total_features, g);
    }
    return gHelper;
  }

  HashMap<Integer, ArrayList<Double>> communities_of_interest;

  private GradientHelper(int total_comms, int total_features, Graph g) {
    this.g = g;
    communities_of_interest = new HashMap<Integer, ArrayList<Double>>();
    for (int commId : g.community_to_vertices_index.keySet()) {
      ArrayList<Double> perFeatureWeight = new ArrayList<Double>();
      communities_of_interest.put(commId, perFeatureWeight);
      for (int j = 0; j < total_features; j++) {
        perFeatureWeight.add(0.0);
      }
      communities_of_interest.put(commId, perFeatureWeight);
    }
  }

}
