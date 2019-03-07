package org.knoesis.results;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

class JaccardC {

  private void readFileInDS(String fileName, HashMap<Integer, ArrayList<String>> ds)
      throws Exception {
    List<String> linesGt = Files.readAllLines(Paths.get(fileName));
    int comm_no = 0;
    for (String line : linesGt) {
      String splits[] = line.split("\t");
      ArrayList<String> nodesBelonging = new ArrayList<String>();
      for (int i = 1; i < splits.length; i++) {
        nodesBelonging.add((splits[i]));
      }
      ds.put(comm_no++, nodesBelonging);
    }
  }

  private double getF1Score(ArrayList<String> detected, ArrayList<String> groundTruth) {
    double true_pos = 0.0;
    HashSet<String> union = new HashSet<String>();
    // [1, 2, 3, 4] -- [3, 6, 7]

    for (String det : detected) {
      for (String gt : groundTruth) {
        union.add(gt);
        union.add(det);
        if (det.equals(gt)) {
          true_pos++;
        }
      }
    }

    return (double) ((double) true_pos / (double) (union.size()));
  }

  public double computeFMeasure(String groundTruth, String detected) throws Exception {
    HashMap<Integer, ArrayList<String>> gtCNMembership = new HashMap<Integer, ArrayList<String>>();
    HashMap<Integer, ArrayList<String>> dtCNMembership = new HashMap<Integer, ArrayList<String>>();

    readFileInDS(groundTruth, gtCNMembership);
    readFileInDS(detected, dtCNMembership);

    double total_f1_score = 0.0;
    for (int dt : dtCNMembership.keySet()) {
      ArrayList<String> detectedListOfNodes = dtCNMembership.get(dt);
      double max_f1_score = 0.0;
      for (int gt : gtCNMembership.keySet()) {
        double f1_score_for_this = getF1Score(detectedListOfNodes, gtCNMembership.get(gt));
        if (f1_score_for_this > max_f1_score) {
          max_f1_score = f1_score_for_this;
        }
      }
      total_f1_score +=
          (double) ((double) (max_f1_score) / (double) (2.0 * (double) dtCNMembership.size()));
    }
    return (total_f1_score);
  }
}


public class CommunityJaccard {
  private static String groundTruth = "lawyer_circle_val.txt";
  private static String detected = "fb_rws_911_updated";

  public double getJaccard(int id, String detected) throws Exception {
    JaccardC fmc = new JaccardC();
    double one_way = fmc.computeFMeasure(detected, "newCircles" + id);
    double other = fmc.computeFMeasure("newCircles" + id, detected);
    double anss = (double) ((double) (one_way + other));
    return anss;
  }

  public double getJaccardGplus(String detected, String groundtruth) throws Exception {
    JaccardC fmc = new JaccardC();
    double one_way = fmc.computeFMeasure(detected, groundtruth);
    double other = fmc.computeFMeasure(groundtruth, detected);
    double anss = (double) ((double) (one_way + other));
    return anss;
  }

  public static void main(String[] args) throws Exception {
    JaccardC fmc = new JaccardC();
    double one_way = fmc.computeFMeasure(detected, groundTruth);
    double other = fmc.computeFMeasure(groundTruth, detected);
    double anss = (double) ((double) (one_way + other));
    System.out.println(anss);
  }
}
