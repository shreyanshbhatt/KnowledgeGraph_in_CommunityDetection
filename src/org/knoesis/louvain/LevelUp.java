package org.knoesis.louvain;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

public class LevelUp {
  private void addNewCommEdge(HashMap<Integer, HashMap<Integer, Double>> newGraphInfo, int own_comm,
      int neigh_comm, double w, boolean should_divide) {
    HashMap<Integer, Double> edgeInfo = null;
    if (newGraphInfo.containsKey(own_comm)) {
      edgeInfo = newGraphInfo.get(own_comm);
    } else {
      edgeInfo = new HashMap<Integer, Double>();
    }
    if (should_divide) // bug fix: after level 1 of louvain iterations, the self
                       // loop community edges won't appear more than once, so
                       // we must not divide those by 2.
      w = (double) ((double) w / (double) 2);
    if (edgeInfo.containsKey(neigh_comm)) {
      w += edgeInfo.get(neigh_comm);
    }
    edgeInfo.put(neigh_comm, w);
    newGraphInfo.put(own_comm, edgeInfo);
  }

  private void printNewGraph(HashMap<Integer, HashMap<Integer, Double>> newGraphInfo,
      String outputFile) throws IOException {
    BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
    for (int source : newGraphInfo.keySet()) {
      HashMap<Integer, Double> edgeInfo = newGraphInfo.get(source);
      for (int target : edgeInfo.keySet()) {
        bw.write(source + " " + edgeInfo.get(target) + " " + target);
        bw.newLine();
      }
    }
    bw.close();
  }

  private HashMap<Integer, HashMap<Integer, Double>> orderAndMapIds(
      HashMap<Integer, HashMap<Integer, Double>> newGraphInfo, HashMap<Integer, Integer> idMap)
      throws IOException {
    HashMap<Integer, HashMap<Integer, Double>> orderedGraphInfo =
        new HashMap<Integer, HashMap<Integer, Double>>();
    int idIncrementer = -1;
    for (int source : newGraphInfo.keySet()) {
      int newSource = 0;
      if (idMap.containsKey(source)) {
        newSource = idMap.get(source);
      } else {
        idIncrementer++;
        newSource = idIncrementer;
      }
      HashMap<Integer, Double> orderedEdgeInfo = new HashMap<Integer, Double>();
      idMap.put(source, newSource);
      HashMap<Integer, Double> existingEdgeInfo = newGraphInfo.get(source);
      for (int dest : existingEdgeInfo.keySet()) {
        int newDest = 0;
        if (idMap.containsKey(dest)) {
          newDest = idMap.get(dest);
        } else {
          idIncrementer++;
          newDest = idIncrementer;
        }
        idMap.put(dest, newDest);
        orderedEdgeInfo.put(newDest, existingEdgeInfo.get(dest));
      }
      orderedGraphInfo.put(newSource, orderedEdgeInfo);
    }

    return orderedGraphInfo;
  }

  // This must be called after prepareForLouvain iteration
  public int prepareForNextLevel(Vertex vertices[], Community[] communities, String outputFile,
      HashMap<Integer, Integer> idMapper, HashMap<Integer, Integer> vertices_to_comm)
      throws IOException {
    HashMap<Integer, HashMap<Integer, Double>> newGraphInfo =
        new HashMap<Integer, HashMap<Integer, Double>>();

    for (Vertex v : vertices) {
      vertices_to_comm.put(v.id, v.community_label);
      for (Neighbor n : v.neighbors) {
        int own_comm = v.community_label;
        int neigh_comm = vertices[n.neighbour_label].community_label;
        if (neigh_comm < own_comm) {
          int temp = neigh_comm;
          neigh_comm = own_comm;
          own_comm = temp;
        }
        addNewCommEdge(newGraphInfo, own_comm, neigh_comm, n.edge_weight,
            (v.id != n.neighbour_label));
      }
    }

    newGraphInfo = orderAndMapIds(newGraphInfo, idMapper);
    printNewGraph(newGraphInfo, outputFile);
    int new_num_vertices = -100000000;

    for (int s : newGraphInfo.keySet()) {
      if (new_num_vertices < s) {
        new_num_vertices = s;
      }
      for (int d : newGraphInfo.get(s).keySet()) {
        if (new_num_vertices < d) {
          new_num_vertices = d;
        }
      }
    }

    return (new_num_vertices + 1);
  }
}
