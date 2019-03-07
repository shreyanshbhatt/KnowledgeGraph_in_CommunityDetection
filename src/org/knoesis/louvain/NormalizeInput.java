package org.knoesis.louvain;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import scala.Tuple3;

public class NormalizeInput {
  public String normalizeInput(String inputFile) throws IOException {
    List<String> lines = Files.readAllLines(Paths.get(inputFile));
    ArrayList<Tuple3<Integer, ArrayList<Double>, Integer>> data =
        new ArrayList<Tuple3<Integer, ArrayList<Double>, Integer>>();

    double divisor = 0.0;
    for (String line : lines) {
      String splits[] = line.split(" ");
      int src = Integer.parseInt(splits[0]);
      int target = Integer.parseInt(splits[splits.length - 1]);
      ArrayList<Double> vals = new ArrayList<Double>();
      for (int i = 1; i < splits.length - 1; i++) {
        double div_val = 0.0;
        if (!splits[i].equals("NaN")) {
          div_val = Double.parseDouble(splits[i]);
        }
        divisor += div_val * div_val;
        vals.add(div_val);
      }
      data.add(new Tuple3<Integer, ArrayList<Double>, Integer>(src, vals, target));
    }
    divisor = Math.sqrt(divisor);
    divisor = (divisor == 0.0) ? 1.0 : divisor;
    String newFileName = inputFile + "_normalized_input";
    BufferedWriter bw = new BufferedWriter(new FileWriter(newFileName));
    for (Tuple3<Integer, ArrayList<Double>, Integer> tpl : data) {
      bw.write("" + tpl._1());
      for (double d : tpl._2()) {
        bw.write(" " + (double) (d / divisor));
      }
      bw.write(" " + tpl._3());
      bw.newLine();
    }
    bw.close();
    return newFileName;
  }

  public static void main(String[] args) throws Exception {
    NormalizeInput ni = new NormalizeInput();
    ni.normalizeInput("testEdgesWeights.gr");
  }
}
