/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package mfruzan.algorithm;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 *
 * @author mario.fruzangohar
 */
public class MapUtil<T> {
    public  void merge(Map<T, Integer> m1, Map<T, Integer> m2){
        // we add elemnts of m2 to m1 
        for(T k: m2.keySet()){
           Integer i = m1.get(k);
           if (i == null)
               m1.put(k, m2.get(k));
           else
               m1.put(k, i + m2.get(k));
        }
    }
    
     public  Map<T, Integer> merge_return(Map<T, Integer> m1, Map<T, Integer> m2){
        // we add elements of m2 to m1 and return the result as a new Map
        Map<T, Integer> map = new HashMap();
        merge(map, m1);
        merge(map, m2);
        
        return map;
        
    }

    public void increment(Map<T, Integer> map, Set<T> set){
        // map holds number of times an object of type T has been visited. We update the map based on element of set
        for(T k : set){
            Integer i = map.get(k);
            if (i == null)
                map.put(k, 1);
            else
                map.put(k, i+1);
        }
    }
    public void increment(Map<T, Integer> map, T k){
        // map holds number of times an object of type T has been visited. We update the map based on element of set
            Integer i = map.get(k);
            if (i == null)
                map.put(k, 1);
            else
                map.put(k, i+1);
    }
    
        public void increment(Map<T, Double> map, T k, double d){
        // map holds number of times an object of type T has been visited. We update the map based on element of set
            Double i = map.get(k);
            if (i == null)
                map.put(k, d);
            else
                map.put(k, i+d);
    }
    
    public int totalCount(Map<T, Integer> map){
        int sum = 0;
        for (T st : map.keySet()){
            sum += map.get(st);
        }
        return sum;
    }
    
    public boolean overlap(Map<Integer, Integer> map, int start, int end){
        // each element of map contains start and end position of multiple genomic regions
        // checks to see if resion specified by start and end ovelpas with any region in the map.
        
        for (int st : map.keySet()){
            int en = map.get(st);
            if ((start >= st && start <= en) || (end >= st && end <=en))
                return true;
        }
        
        return false;
    }
    
    public int notCovered(TreeMap<Integer, Integer> map){
        //map contains genomics regions start of region is the key and end of region value of the map 
        // we assume map is sorted by keys (TreeMap implementation can do that)
        
        int total_not_covered = 0;
        Iterator<Map.Entry<Integer, Integer>> it = map.entrySet().iterator();
        
        int right_end = map.firstEntry().getKey(); // end position of the last entry
        //int left_start = map.firstEntry().getKey();
        while(it.hasNext()){
            Map.Entry<Integer,Integer> entry = it.next();
            int start = entry.getKey();
            int end = entry.getValue();
            if (start > right_end)
                total_not_covered += start-right_end-1;
            if(end>right_end)
                right_end=end;
        }// while in posMap
        
        return total_not_covered;
    }
    
    public int right_end(TreeMap<Integer, Integer> map){
        // returns right most position
        Iterator<Map.Entry<Integer, Integer>> it = map.entrySet().iterator();
        int right_end = 0; 
        //int left_start = map.firstEntry().getKey();
        while(it.hasNext()){
            Map.Entry<Integer,Integer> entry = it.next();
            int end = entry.getValue();
            if (end > right_end)
                right_end=end;
        }// while in posMap
        
        return right_end;        
    }
}
