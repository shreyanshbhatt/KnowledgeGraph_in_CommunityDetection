package org.knoesis.louvain;

import java.util.*;

public class Solution {
  private static void addTransformation(ArrayList<HashSet<Integer>> map, int x, int y) {
    map.get(x).add(y);
    for (int i = 0; i < map.size(); i++) {
      if (map.get(i).contains(x)) {
        map.get(i).add(x);
      }
      if (y == i) {
        map.get(x).add(i);
      }
    }
  }

  private static boolean areEqual(int x, int y, ArrayList<HashSet<Integer>> index) {
    HashSet<Integer> possiblities;

    possiblities = index.get(x);

    if (possiblities.contains(y)) {
      return true;
    }

    for (int transPoss : possiblities) {
      if (index.get(transPoss).contains(y)) {
        return true;
      }
    }

    return false;
  }

  public static void main(String[] args) {
    Scanner in = new Scanner(System.in);
    int n = in.nextInt();
    int k = in.nextInt();
    int m = in.nextInt();
    n++;
    ArrayList<HashSet<Integer>> possibleTrans = new ArrayList<HashSet<Integer>>();

    for (int i = 0; i < n; i++) {
      possibleTrans.add(new HashSet<Integer>());
    }
    for (int a0 = 0; a0 < k; a0++) {
      int x = in.nextInt();
      int y = in.nextInt();
      addTransformation(possibleTrans, x, y);
    }
    int[] a = new int[m];
    for (int a_i = 0; a_i < m; a_i++) {
      a[a_i] = in.nextInt();
    }
    in.close();
    int dynamicHelper[][] = new int[m][m];
    int maxSoFar = 1;
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < m - i; j++) {
        if (i == 0) {
          dynamicHelper[j][i] = 1;
          continue;
        }
        if (a[j] == a[i + j] || areEqual(a[j], a[i + j], possibleTrans)
            || areEqual(a[i + j], a[j], possibleTrans)) {
          dynamicHelper[j][i] = dynamicHelper[j + 1][i - 1] + 2;
        } else {
          dynamicHelper[j][i] = Math.max(dynamicHelper[j + 1][i], dynamicHelper[j][i - 1]);
        }
        if (dynamicHelper[j][i] > maxSoFar) {
          maxSoFar = dynamicHelper[j][i];
        }
      }
    }
    System.out.println(maxSoFar);
  }
}
