package org.knoesis.results;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import JavaMI.MutualInformation;

class NMIG {
  public double getNMI(String file1, String file2) throws Exception {
    List<String> file1lines = Files.readAllLines(Paths.get(file1));
    List<String> file2lines = Files.readAllLines(Paths.get(file2));
    ArrayList<Double> vec1 = new ArrayList<Double>();
    ArrayList<Double> vec2 = new ArrayList<Double>();
    for (String file1line : file1lines) {
      String splits[] = file1line.split(" ");
      for (String split : splits) {
        vec1.add(Double.parseDouble(split));
      }
    }
    for (String file2line : file2lines) {
      String splits[] = file2line.split(" ");
      for (String split : splits) {
        vec2.add(Double.parseDouble(split));
      }
    }
    ArrayList<Double> largerVec = vec1, smallerVec = vec2;
    if (vec1.size() != vec2.size()) {
      largerVec = (vec1.size() > vec2.size()) ? vec1 : vec2;
      smallerVec = (vec1.size() < vec2.size()) ? vec1 : vec2;
      double prev = 0;
      for (int i = 0; i < largerVec.size(); i++) {
        if (i < smallerVec.size()) {
          prev = smallerVec.get(i);
          continue;
        }
        smallerVec.add(prev);
      }
    }
    if (largerVec.size() != smallerVec.size()) {
      System.out.println("Issueeee");
    }
    double[] vecarr1 = new double[largerVec.size()];
    double[] vecarr2 = new double[smallerVec.size()];
    for (int i = 0; i < largerVec.size(); i++) {
      vecarr1[i] = largerVec.get(i);
      vecarr2[i] = smallerVec.get(i);
    }
    return MutualInformation.calculateMutualInformation(vecarr1, vecarr2);
  }
}


public class NMIComputer {
  public double computeNMI(int id) throws Exception {
    NMIG nmig = new NMIG();
    double total = nmig.getNMI("for_nmi_" + id, "for_nmi_" + id);
    double for_this = nmig.getNMI("for_nmi_" + id, "for_nmi_veredges_" + id);
    return (for_this / total);
  }

  public double computeNMIGplus(String gt, String det) throws Exception {
    NMIG nmig = new NMIG();
    double total = nmig.getNMI(gt, gt);
    double for_this = nmig.getNMI(det, gt);
    return (for_this / total);
  }

  public static void main(String[] args) throws Exception {
    NMIG nmig = new NMIG();
    int id = 686;
    System.out.println(nmig.getNMI("for_nmi_" + id, "for_nmi_veredges_" + id));
  }
}
