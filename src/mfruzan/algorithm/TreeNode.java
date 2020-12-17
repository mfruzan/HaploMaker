/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package mfruzan.algorithm;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;

/**
 *
 * @author mario.fruzangohar
 */
public class TreeNode<T extends Comparable> {
    TreeNode leftChild = null;
    TreeNode rightChild = null;
    Map<TreeNode, Double> children = new HashMap<TreeNode, Double>();  
    TreeNode  parent = null;
    boolean hasParent = false;
    T content = null; 
    double total_evidence = 0; 
    TreeNode backref= null; 
    double min_weight = Double.MAX_VALUE; 
    TreeNode backref2= null; 
    double min_weight2 = Double.MAX_VALUE; 
    
    public TreeNode(T content){
        this.content = content;
        
    }
     public TreeNode(){
        
    }
     
    
     public double getWeight(TreeNode child){         
          double i = children.get(child).doubleValue(); 
          double pr = 0;
          if (i==0 && total_evidence>0){
              pr = (double)1/(total_evidence*2);
          }else if (total_evidence == 0){
              pr = (double)1/children.size();  
          }else{
              pr = (double)i/total_evidence;
          }
          
          return -Math.log10(pr);
     }
     
     public int hashCode() {
        return new HashCodeBuilder(11, 73). // two randomly chosen prime numbers
            append(content). // content itself overrides             
            toHashCode();
        
    }
    
    public boolean equals(Object other){
        if (other == null) return false;
        if (other == this) return true;
        if (!(other instanceof TreeNode))return false;
        TreeNode node = (TreeNode)other;
        return new EqualsBuilder().append(content, node.content).isEquals();
        
    }
}
