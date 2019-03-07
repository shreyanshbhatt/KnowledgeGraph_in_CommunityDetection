package org.knoesis.louvain;

import java.io.IOException;
import java.util.ArrayList;

public class GradientDescent {
  private static int GRAD_MAX = 90;
  int numCommunities;
  int numTopics;

  public GradientDescent(int numCommunities, int numTopics) {
    this.numCommunities = numCommunities;
    this.numTopics = numTopics;
  }

  private void copyPreviousIterationWeights(Graph g, GradientHelper newWeightsAll) {
    for (int commId : g.community_to_vertices_index.keySet()) {
      for (int k = 0; k < g.communities[commId].topic_weights.size(); k++) {
        newWeightsAll.communities_of_interest.get(commId).set(k,
            g.communities[commId].topic_weights.get(k));
      }
    }
  }

  // public void updateGh(GradientHelper gh, int k, Graph g, double step, double
  // old_comm_weights[]) {
  // CommunityHelper ch[] = gh.communities_of_interest[k];
  //
  // for (int i = 0; i < ch.length; i++) {
  // old_comm_weights[i] = g.communities[k].topic_weights.get(i);
  // double update_val = (double)((double)g.communities[k].topic_weights.get(i)
  // + (double)(ch[i].gradient * step));
  // update_val = (update_val < 0)? 0 : update_val;
  // g.communities[k].topic_weights.set(i, update_val);
  // }
  // }

  // public void replaceWeights(Graph g, double old_comm_weights[], int k,
  // GradientHelper gh) {
  // for (int i = 0; i < old_comm_weights.length; i++) {
  // g.communities[k].topic_weights.set(i, old_comm_weights[i]);
  // gh.communities_of_interest[k][i].actual_weight = old_comm_weights[i];
  // }
  // }

  // public boolean getNormSatisfying(ArrayList<Double> theta_new,
  // CommunityHelper []grad, double step) {
  // double squared_norm_2_theta_new = 0.0;
  // double squared_norm_2_grad = 0.0;
  // for (int i = 0; i < theta_new.size(); i++) {
  // squared_norm_2_theta_new += (double) (theta_new.get(i) * theta_new.get(i));
  // squared_norm_2_grad += (double) (grad[i].gradient * grad[i].gradient);
  // }
  // squared_norm_2_theta_new = Math.sqrt(squared_norm_2_theta_new);
  // squared_norm_2_grad = Math.sqrt(squared_norm_2_grad);
  // return (squared_norm_2_grad * step < 0.05 * squared_norm_2_theta_new);
  // }

  private double[] getCopyOfCurWeights(GradientHelper gh, int k) {
    double[] old_weights = new double[gh.communities_of_interest.get(k).size()];
    for (int i = 0; i < gh.communities_of_interest.get(k).size(); i++) {
      old_weights[i] = gh.communities_of_interest.get(k).get(i);
    }
    return old_weights;
  }

  private void updateCurrentCommunityWeightVector(ArrayList<Double> gradients, Graph g,
      GradientHelper newWeightsAll, int comm, double step) {
    // update the new weight to newWeightsAll and to the graph
    for (int i = 0; i < gradients.size(); i++) {
      double existing_weight = newWeightsAll.communities_of_interest.get(comm).get(i);
      existing_weight += step * gradients.get(i);
      if (existing_weight < 0) {
        existing_weight = 0.0;
      }
      newWeightsAll.communities_of_interest.get(comm).set(i, existing_weight);
      g.communities[comm].topic_weights.set(i, existing_weight);
    }
  }

  private void replaceWeights(Graph g, GradientHelper newWeightsAll, double[] copyCurrentWeights,
      int k) {
    for (int i = 0; i < copyCurrentWeights.length; i++) {
      newWeightsAll.communities_of_interest.get(k).set(i, copyCurrentWeights[i]);
      g.communities[k].topic_weights.set(i, copyCurrentWeights[i]);
    }
  }

  private void replaceGraphWeights(Graph g, double[] copyOfOldWeights, int k) {
    for (int i = 0; i < copyOfOldWeights.length; i++) {
      g.communities[k].topic_weights.set(i, copyOfOldWeights[i]);
    }
  }

  private void replaceNewWeightVectorGraph(GradientHelper newWeightsAll, Graph g) {
    for (int comm : newWeightsAll.communities_of_interest.keySet()) {
      for (int topic = 0; topic < newWeightsAll.communities_of_interest.get(comm).size(); topic++) {
        g.communities[comm].topic_weights.set(topic,
            newWeightsAll.communities_of_interest.get(comm).get(topic));
      }
    }
  }

  public void performGradientDescent(Graph g) throws Exception {
    double old_gain = g.computeGainAll();

    int count = 0;
    int converge = 0;
    GradientHelper newWeightsAll =
        GradientHelper.getInstance(g.community_to_vertices_index.size(), numTopics, g);
    copyPreviousIterationWeights(g, newWeightsAll);

    // Currently we are doing it in a very inefficient way where we optimize for
    // each community
    // at once and for each community we traverse whole graph.
    for (int k : g.community_to_vertices_index.keySet()) {
      double step = 0.05;
      // now optimize according to each feature
      if (!g.communities[k].isActive) {
        continue;
      }
      if (g.communities[k].total_vertices < 2) {
        continue;
      }
      count = 0;
      converge = 0;
      double[] copyOfOldWeights = getCopyOfCurWeights(newWeightsAll, k);
      double old_gain_temp = old_gain;
      while (count < GRAD_MAX && converge < 5) {
        ArrayList<Double> gradients = g.computeGradientOfCommunity(k);
        // update community weight vector
        double[] copyCurrentWeights = getCopyOfCurWeights(newWeightsAll, k);
        updateCurrentCommunityWeightVector(gradients, g, newWeightsAll, k, step);
        g.updateWeightsDegreeGradient(k);
        // compute new gain
        double new_gain = g.computeGainOfCommunity(k);
        // System.out.println("old = "+old_gain[k]+" == "+new_gain);
        if (Double.compare(new_gain, old_gain_temp) > 0) {
          old_gain_temp = new_gain;
          count++;
          // if (getNormSatisfying(g.communities[k].topic_weights,
          // gh.communities_of_interest[k], step)) {
          // converge++;
          // } else {
          // converge = 0;
          // }
          step *= (double) 1.2;
        } else {
          count++;
          converge = 0;
          step = step * 0.5;
          replaceWeights(g, newWeightsAll, copyCurrentWeights, k);
          g.updateWeightsDegreeGradient(k);
        }
      }
      replaceGraphWeights(g, copyOfOldWeights, k);
    }
    replaceNewWeightVectorGraph(newWeightsAll, g);
    g.printCommAffiliations();
  }

}
