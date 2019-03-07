package org.knoesis.results;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

class FMeasureCesna {
  public void readFileInDS(String fileName, HashMap<Integer, ArrayList<String>> ds)
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

  private double getPrecision(double true_pos, double false_pos) {
    return (double) ((double) true_pos / (double) (true_pos + false_pos));
  }

  private double getRecall(double true_pos, double false_neg) {
    return (double) ((double) true_pos / (double) (true_pos + false_neg));
  }

  private double getF1Score(ArrayList<String> detected, ArrayList<String> groundTruth) {
    double true_pos = 0.0;
    double false_pos = (double) detected.size();
    double false_neg = (double) groundTruth.size();
    // [1, 2, 3, 4] -- [3, 6, 7]

    for (String det : detected) {
      for (String gt : groundTruth) {
        if (det.equals(gt)) {
          true_pos++;
          false_pos--;
          false_neg--;
        }
      }
    }

    if (true_pos == 0) {
      return 0.0;
    }
    double prec = getPrecision(true_pos, false_pos);
    double rec = getRecall(true_pos, false_neg);
    double precRecProd = (double) ((double) prec * (double) rec);
    double precRecSum = (double) ((double) prec + (double) rec);
    return (double) (((double) 2 * precRecProd) / precRecSum);
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
        if (Double.compare(f1_score_for_this, max_f1_score) > 0) {
          max_f1_score = f1_score_for_this;
        }
      }
      // System.out.println("\t\t\t--"+max_f1_score);
      // total_f1_score += (double)((max_f1_score * dtCNMembership.size())/
      // (double)((double)total_nodes));
      total_f1_score += (double) ((max_f1_score) / (double) (2.0 * (double) dtCNMembership.size()));
    }
    return (total_f1_score);
  }
}


public class CommunityFMeasureCesna {
  private static String groundTruth = "lawyer_circle_val.txt";
  private static String detected = "fb_rws_911_updated";

  public double measureFScore(int id, String detected) throws Exception {
    FMeasureCesna fmc = new FMeasureCesna();
    double one_way = fmc.computeFMeasure("newCircles" + id, detected);
    double other = fmc.computeFMeasure(detected, "newCircles" + id);
    double anss = (double) ((double) (one_way + other));
    return anss;
  }

  public double measureFScoreGplus(String detected, String groundtruth) throws Exception {
    FMeasureCesna fmc = new FMeasureCesna();
    double one_way = fmc.computeFMeasure(detected, groundtruth);
    double other = fmc.computeFMeasure(groundtruth, detected);
    double anss = (double) ((double) (one_way + other));
    return anss;
  }

  public int diffNum(String detected, String groundTruth) throws Exception {
    FMeasureCesna fmc = new FMeasureCesna();
    HashMap<Integer, ArrayList<String>> gtCNMembership = new HashMap<Integer, ArrayList<String>>();
    HashMap<Integer, ArrayList<String>> dtCNMembership = new HashMap<Integer, ArrayList<String>>();
    fmc.readFileInDS(groundTruth, gtCNMembership);
    fmc.readFileInDS(detected, dtCNMembership);
    return (gtCNMembership.size() - dtCNMembership.size());
  }

  public static void main(String[] args) throws Exception {
    FMeasureCesna fmc = new FMeasureCesna();
    double one_way = fmc.computeFMeasure(groundTruth, detected);
    double other = fmc.computeFMeasure(detected, groundTruth);
    double anss = (double) ((double) (one_way + other));
    System.out.println(anss);
  }
}
