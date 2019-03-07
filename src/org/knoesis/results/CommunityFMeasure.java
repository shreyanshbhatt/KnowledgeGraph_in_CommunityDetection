package org.knoesis.results;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;

class FMeasureC {

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

  private int getSizes(HashMap<Integer, ArrayList<String>> map) {
    int total = 0;
    for (int dt : map.keySet()) {
      total += map.get(dt).size();
    }
    return total;
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

    int total_size = 0;
    total_size = getSizes(dtCNMembership);

    double total_f1_score = 0.0;
    for (int dt : dtCNMembership.keySet()) {
      ArrayList<String> detectedListOfNodes = dtCNMembership.get(dt);
      double size_ratio = (double) ((double) detectedListOfNodes.size() / (double) (total_size));
      double max_f1_score = 0.0;
      for (int gt : gtCNMembership.keySet()) {
        double f1_score_for_this = getF1Score(detectedListOfNodes, gtCNMembership.get(gt));
        if (f1_score_for_this > max_f1_score) {
          max_f1_score = f1_score_for_this;
        }
      }
      // System.out.println("got max = "+max_f1_score);
      total_f1_score += (double) (size_ratio * max_f1_score);
    }
    return total_f1_score;
  }
}


public class CommunityFMeasure {
  private static String id = "3980";

  public double measureFScore(int id, String detected) throws Exception {
    FMeasureC fmc = new FMeasureC();
    double one_way = fmc.computeFMeasure("newCircles" + id, detected);
    // double other = fmc.computeFMeasure("foundLabelsedges_"+id,
    // "newCircles"+id);
    // double anss = (double)((double)(one_way + other));
    return one_way;
  }

  public double measureFScoreGplus(String detected, String groundtruth) throws Exception {
    FMeasureC fmc = new FMeasureC();
    double one_way = fmc.computeFMeasure(groundtruth, detected);
    // double other = fmc.computeFMeasure("foundLabelsedges_"+id,
    // "newCircles"+id);
    // double anss = (double)((double)(one_way + other));
    return one_way;
  }

  // private static String groundTruth = "newCircles"+id;
  // private static String detected = "foundLabelsedges_"+id;
  private static String groundTruth = "lawyer_circle_val.txt";
  private static String detected = "fb_rws_911_updated";

  public static void main(String[] args) throws Exception {

    // String groundTruth = "newCirclesCora";
    //
    // for (final File fileEntry : new File("./").listFiles()) {
    // String detected = fileEntry.getName();
    // if (!detected.startsWith("foundLabelsedges_cora")) {
    // continue;
    // }
    FMeasureC fmc = new FMeasureC();
    // System.out.println(detected+" : "+fmc.computeFMeasure(groundTruth,
    // detected));
    System.out.println(detected + " : " + fmc.computeFMeasure(groundTruth, detected));

  }
}
