/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package mfruzan.algorithm;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Stack;
import mfruzan.common.Helper;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.Collections;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author mario.fruzangohar
 */

public class HaplotypeDAG {
    private Map<Integer, List<TreeNode>> graphNodes = new HashMap(); // Level -> List of TreeNodes
    private int depth = 0;
    
    private int start_level = 1; //added 8/5/2019, we start looking for haplotype from this level
    public List<SNP> far_evidence_from = new ArrayList(); 
    public List<SNP> far_evidence_to = new ArrayList();
    public List<FarEvidence> far_evidence = new ArrayList();
    public List<Boolean> processed = new ArrayList();// this is a list with equal size of far_evidence
    List<GraphEdge> new_edges = new ArrayList();
    public HaplotypeDAG(int depth){    
        this.depth = depth;
        for (int i=0; i<=depth; i++)
            graphNodes.put(i, new ArrayList());
        TreeNode root = new TreeNode(); // empty node, root has level zero, (its content also remains null)
        graphNodes.get(0).add(root);  // first Map (key=0) has got just one Node and that node is the root
    }
   
   
    public void clear(){
        //root = null;
        //pointer = null;
    }
    
    public void insert(List<SNP> snpChain, double weight){
         
        TreeNode parent = null;
        for (SNP snp: snpChain){
            if (snp.allele==null)
                continue;// if a snp is not resolved we dont add it to the DAG

            if (snp.orderID == 1)
                parent = graphNodes.get(0).get(0); // root node should be the parent node
            
            TreeNode node = getNode(snp);
            if(node == null){
                node= new TreeNode(snp); // content of TreeNode is snp
                graphNodes.get(snp.orderID).add(node);
            }
            if(parent != null){
                new MapUtil().increment(parent.children, node, weight);
                parent.total_evidence = weight;
                node.hasParent = true;
            }           
        } // for
    }
    
    public void resolve_far_evidence() throws Exception{
        for(int i=0; i<far_evidence_from.size(); i++){
            if (far_evidence_from.get(i).contig.equals("chr3") && far_evidence_from.get(i).position>=27201554 && far_evidence_from.get(i).position<=27204109)
                System.out.println(String.format("%d %s To %d %s", far_evidence_from.get(i).position, far_evidence_from.get(i).allele, far_evidence_to.get(i).position, far_evidence_to.get(i).allele));
            if (far_evidence_to.get(i).orderID - far_evidence_from.get(i).orderID <= 30 ){
                traverse_partial(getNode(far_evidence_from.get(i)), new Stack(), getNode(far_evidence_to.get(i)), 1);
            }else
                System.out.println("More than 30 levels far evidence, This fragment will be skipped.");
        }
        
    }
    public void initFarEvidence(){
        Collections.sort(far_evidence);
        processed.clear();
        for(int i=0; i<far_evidence.size(); i++){
            processed.add(false) ; 
            
        }        
    }
    public int resolve_far_evidence_Basic() throws Exception{
        int cc = 0;
        boolean new_path_added = false;
            do{
                new_path_added = false;            
                for(int i=0; i<far_evidence.size(); i++){

                    if(processed.get(i))
                        continue;
                    List<List<SNP>> paths = new ArrayList();
                    List<SNP> deadend = new ArrayList();
                    traverse_partial_V2(getNode(far_evidence.get(i).from), new Stack(), getNode(far_evidence.get(i).to),paths, deadend);
                    if(paths.size()==1){
                        insert(paths.get(0), far_evidence.get(i).weight);
                        
                        processed.set(i, true);
                        
                        new_path_added = true;                        
                        cc++;
                    }

                }//for each far evidence
            }while(new_path_added );
        
        return cc;
    }
   
    
    public int  finalConnect(){
        int cc = 0;
        
        return cc;
    }
    
    private TreeNode getNode(SNP snp){
        for(TreeNode node : graphNodes.get(snp.orderID))
            if (node.content!=null && node.content.compareTo(snp)==0)
                return node;
        
        return null;
    }
    
    public SNP getAlternateSNP(SNP snp){
        //For example if genotype at a SNP is GD and snp has got allele G, the allele D will be returned
        for(TreeNode node : graphNodes.get(snp.orderID))
           if (((SNP)node.content).allele!=null && node.content.compareTo(snp)!=0)
                return (SNP)node.content; 
        // if we have only one node at this level return itself
        return snp;
    }
    public TreeNode getAlternateNode(TreeNode aNode){
        //For example if genotype at a SNP is GD and snp has got allele G, the allele D will be returned
        SNP aSNP = (SNP)aNode.content;
        for(TreeNode node : graphNodes.get(aSNP.orderID)){
            SNP snp = (SNP) node.content;
            if (snp!=null && snp.compareTo(aSNP)!=0)
                return node;
        }
 
        
        return null;
    }
    
    public List<SNP>  getAlternateSNPList(List<SNP> list){
        List<SNP> out = new ArrayList();
        for(SNP snp : list){
            SNP otherSNP = snp.clone();
            String[] als = snp.genotype.split("/");
            if (als.length>1){
                if (als[0].equals(snp.allele))
                    otherSNP.allele = als[1];
                else
                    otherSNP.allele = als[0];
            }
            otherSNP.setVCFAlleleCode();
            out.add(otherSNP);
        }
        
        return out;
    }
   
 
    
    

    

   
   
    public int getEmptyLevel(){
        int cc = 0;
        // check for empty levels returns 0, if not such a Node was found
        for(int i=1; i<=depth; i++)
            if(graphNodes.get(i).size() == 0)
                cc++;
        return cc;
    }
    public int getOneLevel(){
        // check for levels have only one node returns 0, if not such a Node was found
        int cc = 0;    
        for(int i=1; i<=depth; i++)
            if(graphNodes.get(i).size() == 1)
                cc++;
        return cc;
    }
    
    public void writeFullLevels(String file){
        Writer writer = null;
        try{
            writer = new FileWriter(file, true);//append mode is true
            for(int i=1; i<=depth; i++)
            if(graphNodes.get(i).size() == 2){
                SNP snp = (SNP)graphNodes.get(i).get(0).content;
                
                writer.write(String.format( "%s\t%d\n", snp.contig, snp.position));
            }

        }catch(Exception ex){
        }finally{
            try {
                writer.close();
            } catch (IOException ex) {}
        }
        
    }
    public List<Integer> connectBrokenLevels(){
        List<Integer> brokenLevels = new ArrayList();
        for(int i=1; i<depth; i++){
            boolean hasChild = false;
            for(TreeNode node : graphNodes.get(i)){
                if (node.total_evidence>0)// so if children map has a key-value but total_evidence is zero, it still considered as broken (zero children not counted)
                    hasChild = true;
            }
            if (graphNodes.get(i).size()>0 && !hasChild){
                brokenLevels.add(i);
                for(TreeNode node_i : graphNodes.get(i)){
                    for(TreeNode node_i_1 : graphNodes.get(i+1)){
                        node_i.children.put(node_i_1, 0);
                        node_i_1.hasParent = true;
                    }
                }
            }
                
        } // for
        
        return brokenLevels;
    }

 
   
    
    private boolean isLevelBroken(int level){
        if (graphNodes.get(level).size()==0)
            return false;
        if (level<depth && graphNodes.get(level+1).size()==0)
            return false;
        for(TreeNode node : graphNodes.get(level))
            if(node.total_evidence > 0)
                return false;
        
        return true;
        
    }
    
    private void connectBrokenLevels(List<Integer> levels, boolean ref_to_ref){
        //developed on 30/4/2020
        // This method looks between levels start and end and if finds any broken level, connect ref snp to ref/non_ref snp based on ref_to_ref value
        for(int lv : levels){
            if(ref_to_ref ){
                // find reference SNP at level and reference SNP at level+1
                TreeNode up_ref = null;
                TreeNode up_non_ref = null;
                TreeNode down_ref = null;
                TreeNode down_non_ref = null;
                
                for(TreeNode node : graphNodes.get(lv))
                    if (((SNP)node.content).isVCFReference()){
                        up_ref = node; 
                        break;
                    }
                for(TreeNode node : graphNodes.get(lv+1))
                    if (((SNP)node.content).isVCFReference()){
                        down_ref = node; 
                        break;
                    }
                
                //now get alternate alleles
                if (up_ref != null)
                   up_non_ref = getAlternateNode(up_ref);
                if(down_ref!=null)
                    down_non_ref = getAlternateNode(down_ref);
                
                // now connect edges
                if(up_ref != null && down_ref != null){
                    new MapUtil().increment(up_ref.children, down_ref); 
                    up_ref.total_evidence++; 
                    down_ref.hasParent = true;
                }
                if(up_non_ref != null && down_non_ref != null){
                    new MapUtil().increment(up_non_ref.children, down_non_ref); 
                    up_non_ref.total_evidence++; 
                    down_non_ref.hasParent = true;                    
                }
                
            }
            
        }//for each level
        
    }
    
   
    
    
    
    
    
    
    public List<Integer> countConnectedLevels(List<Integer> levels){
        List<Integer> list = new ArrayList();
        for(int level : levels){
            //System.out.println("Level " + level);
            for(TreeNode node : graphNodes.get(level)){
              if (node.total_evidence>0){
                  list.add(level);
                  break;
              }
                  
            }  // for each node in the level             
        } // for each broken level
        
        return list;
    }
    
    
    
    
    
    
   
    
   
    
    
    public void updateHasParent(){
        for(int i=1; i<depth; i++){
             // we could use has_parent feature instead
            for(TreeNode node : graphNodes.get(i)){
                for(Object child: node.children.keySet()){
                    TreeNode c_child = (TreeNode)child;
                    if((int)node.children.get(c_child) > 0)
                        c_child.hasParent = true;
                }
            }                
        } // for
    }
    public int removeIncompleteNodes(){
        //This method removes any node that has not children (incompelte); start from bottom
         int cnt = 0; // number of incomlete paths
        for(int i=depth-1; i>=1; i--){
            Iterator<TreeNode> iterator = graphNodes.get(i).iterator();
            while(iterator.hasNext()){
                TreeNode node = iterator.next();
                if (node.children.isEmpty()){
                   //first remove all the references from this nodes in its parents
                    List<TreeNode> potential_parents =  graphNodes.get(i-1);
                    for(TreeNode parent : potential_parents){
                        if (parent.children.containsKey(node)){
                            parent.total_evidence -= (Integer)parent.children.get(node);
                            parent.children.remove(node);
                        }
                    }
                    //now remove the node itself
                    iterator.remove();
                    cnt++;
                }
            }//while
        }
        
        return cnt;
    }
    
    
    public int getBrokenPointsNum(){
        int cc = 0;
        for(int i=1; i<depth; i++){
            boolean hasChild = false;
            for(TreeNode node : graphNodes.get(i)){
                if (node.total_evidence>0)// so if children map has a key-value but total_evidence is zero, it still considered as broken 
                    hasChild = true;
            }
            if (!hasChild)
                cc++;
        } // for
        return cc;
    }
    public void fillLevels(List<SNP> list){
        if (list.size() != depth)
            return;
        for(int i=1; i<=depth; i++){
            if(graphNodes.get(i).size()==1){
                SNP snp = (SNP)graphNodes.get(i).get(0).content;
                SNP otherSNP = snp.clone();
                String[] als = snp.genotype.split("/");
                if (als.length>1){
                    if (als[0].equals(snp.allele))
                        otherSNP.allele = als[1];
                    else
                        otherSNP.allele = als[0];
                }
                otherSNP.setVCFAlleleCode();
                graphNodes.get(i).add(new TreeNode(otherSNP));
            }else if(graphNodes.get(i).size()==0){
              SNP snp = list.get(i-1);
              String[] als = snp.genotype.split("/");
              if(als.length>1){
                  SNP snp1 = snp.clone();
                  SNP snp2 = snp.clone();
                  snp1.allele = als[0];
                  snp2.allele = als[1];
                  snp1.setVCFAlleleCode();
                  snp2.setVCFAlleleCode();
                   graphNodes.get(i).add(new TreeNode(snp1));
                   graphNodes.get(i).add(new TreeNode(snp2));
                  
              }
            }
        }
    }
    
    public List<Integer> getBrokenLevels(){
        List<Integer> list  = new ArrayList();
        for(int i=1; i<depth; i++){
             // we could use has_parent feature instead           
            if (isLevelBroken(i))
                list.add(i);
        } // for
        return list;
    }
    
    public int markBrokenPoints(List<SNP> list) throws Exception{
        int cc = 0;
        if (list.size()!= depth)
            throw new Exception("Exception inside detectBrokenPoints");
        for(int i=1; i<depth; i++){
            boolean hasChild = false;
            for(TreeNode node : graphNodes.get(i)){
                if (node.total_evidence>0)// so if children map has a key-value but total_evidence is zero, it still considered as broken 
                    hasChild = true;
            }
            if (!hasChild){
                list.get(i-1).broken = true;
                cc++;
            }
        } // for
        return cc;
    }
    
    
    
    public void printStructure(){
        for(int i =0; i<=depth; i++){
           List<TreeNode> nodes = graphNodes.get(i);
           System.out.println("level " + i + " of graph has got " + nodes.size() + " Nodes.");
           for(TreeNode node: nodes){
               String allele = "";
               if (node.content==null)
                   allele = "root";
               else 
                   allele = ((SNP)node.content).allele;
               
               System.out.print("Node:"+ allele + " has got " +  node.children.size() + " children.");
               for(Object child : node.children.keySet()) System.out.print(  ((SNP)((TreeNode)child).content).allele  + ":" +  node.children.get(child) + " ");
               
               System.out.println();
           }
        }
        
    }
    
   public void printStructureSimple(){
        for(int i =1; i<depth; i++){
           List<TreeNode> nodes = graphNodes.get(i);
           List<TreeNode> next_nodes = graphNodes.get(i+1);
           
           TreeNode nodes_left = nodes.get(0);
           
           TreeNode nodes_right = null;
           if (nodes.size()>1)
              nodes_right = nodes.get(1);
           
           TreeNode next_nodes_right = null;
           TreeNode next_nodes_left = next_nodes.get(0);
           if(next_nodes.size()>1)
              next_nodes_right = next_nodes.get(1);
           
           StringBuffer labels = new StringBuffer();
           labels.append(((SNP)nodes_left.content).position+"\t");
           labels.append(((SNP)nodes_left.content).allele);
           labels.append("  ");
           if(nodes_right!=null)
               labels.append(((SNP)nodes_right.content).allele);
           
           System.out.println(labels.toString());
           
           StringBuffer edges = new StringBuffer();
           edges.append(((SNP)nodes_left.content).position+"\t");
           
            if (nodes_left.children.containsKey(next_nodes_left))
                edges.append('|');
            else 
                edges.append(' ');
            
            if (next_nodes_right!=null && nodes_left.children.containsKey(next_nodes_right))
                edges.append('\\');
            else
                edges.append(' ');
            
           
            
            if(nodes_right!=null && nodes_right.children.containsKey(next_nodes_left)){
                edges.append('/');
            }else{
                 edges.append(' ');
            }
            
           if(nodes_right!=null && next_nodes_right!=null && nodes_right.children.containsKey(next_nodes_right)){
                edges.append('|');
            } 
           
           System.out.println(edges.toString());
           
        } // for each level
        //now we need to print labels of last level 
        List<TreeNode> nodes = graphNodes.get(depth);
        TreeNode nodes_left = nodes.get(0);
        TreeNode nodes_right = null;
        if (nodes.size()>1)
           nodes_right = nodes.get(1);
        
        StringBuffer labels = new StringBuffer();
        labels.append(((SNP)nodes_left.content).allele);
        labels.append("  ");
        if(nodes_right!=null)
            labels.append(((SNP)nodes_right.content).allele);           
        System.out.println(labels.toString());
        
        
    }
   public void printStructureSimple(String cntg, int start, int end){
        for(int i =1; i<depth; i++){
           List<TreeNode> nodes = graphNodes.get(i);
           List<TreeNode> next_nodes = graphNodes.get(i+1);
           
           TreeNode nodes_left = nodes.get(0);
           if ( !((SNP)nodes_left.content).contig.equals(cntg) || ((SNP)nodes_left.content).position<start || ((SNP)nodes_left.content).position>end)
               continue;
           TreeNode nodes_right = null;
           if (nodes.size()>1)
              nodes_right = nodes.get(1);
           
           TreeNode next_nodes_right = null;
           TreeNode next_nodes_left = next_nodes.get(0);
           if(next_nodes.size()>1)
              next_nodes_right = next_nodes.get(1);
           
           StringBuffer labels = new StringBuffer();
           
           labels.append(((SNP)nodes_left.content).position+"\t");
           labels.append(((SNP)nodes_left.content).allele);
           labels.append("  ");
           if(nodes_right!=null)
               labels.append(((SNP)nodes_right.content).allele);
           
           System.out.println(labels.toString());
           
           StringBuffer edges = new StringBuffer();
           StringBuffer counts = new StringBuffer();
           edges.append(((SNP)nodes_left.content).position+"\t");
           
            if (nodes_left.children.containsKey(next_nodes_left)){
                edges.append('|');
                counts.append(String.format("%.2f,", (double)nodes_left.children.get(next_nodes_left)));
            }else 
                edges.append(' ');
            
            if (next_nodes_right!=null && nodes_left.children.containsKey(next_nodes_right)){
                edges.append('\\');
                counts.append(String.format("%.2f,", (double)nodes_left.children.get(next_nodes_right)));
            }else
                edges.append(' ');
            
           
            
            if(nodes_right!=null && nodes_right.children.containsKey(next_nodes_left)){
                edges.append('/');
                //counts.append((int)nodes_right.children.get(next_nodes_left) + ",");
                counts.append(String.format("%.2f,", (double)nodes_right.children.get(next_nodes_left)));
            }else{
                 edges.append(' ');
            }
            
           if(nodes_right!=null && next_nodes_right!=null && nodes_right.children.containsKey(next_nodes_right)){
               
                edges.append('|');
                //counts.append((int)nodes_right.children.get(next_nodes_right) + ",");
                counts.append(String.format("%.2f,", (double)nodes_right.children.get(next_nodes_right)));
            } 
           
           System.out.println(edges.toString() + "\t" + counts.toString());
           
        } // for each level
       
        
        
    }
   
   public Stack<SNP> findNextHaplotype(){
       int step = 1;
       Stack<SNP> hap1 = new Stack();
       
       if (start_level >= depth-1)//we already have reached to the end
           return null;
       step = 2;
       int end_level = findNextShortestPath();
       step = 3;
       if(end_level>start_level){
            TreeNode closestNode = null;
            double min_dist = Double.MAX_VALUE;
            for(TreeNode node : graphNodes.get(end_level)){
                //System.out.println("Node " + ((SNP)node.content).allele + " : " + node.min_weight);
                if (node.min_weight < min_dist){
                    min_dist = node.min_weight;
                    closestNode = node;
                }
            }
            TreeNode prevNode = closestNode;
            while (prevNode.content != null){
                hap1.push((SNP)prevNode.content);
                prevNode = prevNode.backref;           
            }           
       }
      
       step = 4;
       if (end_level == depth)
          start_level = depth;
       else{
           start_level = end_level + 1;
           for(TreeNode node : graphNodes.get(end_level))
               if (node.children.size()>0){
                   start_level = end_level;
                   break;
               }
           
       }
       
       return hap1;
       
   }
   
   
    
    
    private double sumDistance(Stack<Double> dists){
        double sum = 0;
        for(double d : dists)
            sum += d;
        
        return sum;
    }
    
    
    
    private int findNextShortestPath(){
 int step = 1;
       
       while(graphNodes.get(start_level).size() == 0)
           start_level++;
step = 2;       
        if(start_level == 1)
            graphNodes.get(0).get(0).min_weight = 1;
        else if(start_level>1){
            TreeNode newRoot = new TreeNode(); // empty node, (its content also remains null)
            newRoot.min_weight = 1;
            double total = 0;
            for(TreeNode node : graphNodes.get(start_level)){
                node.min_weight = Double.MAX_VALUE;
                
                node.hasParent = true;
                
                total++;
            }
            newRoot.total_evidence = total; //it must be 2
            graphNodes.get(start_level+1).clear();
            graphNodes.get(start_level+1).add(newRoot);   
            
        }
step = 3;        
        for(int i=start_level-1; i<=depth; i++){
            //System.out.println("level " + i);
            boolean stuck = true;
            for(TreeNode node : graphNodes.get(i)){
                if (node.min_weight != Double.MAX_VALUE){
                    stuck = false;
                    for(Object child : node.children.keySet()){
                        double sum_weight = node.min_weight + node.getWeight((TreeNode)child);
                        if (sum_weight < ((TreeNode)child).min_weight){
                            ((TreeNode)child).min_weight = sum_weight;
                            
                        }
step = 5;                        
                    }
                }
            }// for each node at level i
            
            if (stuck)
                return i+1;
        }
       
        return depth;
    }
    
    

    private void traverse_partial(TreeNode node, Stack<SNP> snps, TreeNode finalNode, double weight) throws Exception{
          try{
          
          if (node.content != null){ // only root content is null
              SNP snp = (SNP)node.content;
              snps.push(snp);
          }
          
          if(node.content.compareTo(finalNode.content)==0){
              
              insert(snps.subList(0, snps.size()), weight);
              return;
          }
          
          if (((SNP)node.content).orderID>=((SNP)finalNode.content).orderID)
                return;
          
          if (node.children.size()==0)  // we should not use node.total_evidence==0 , because a node with zero total_evidence can still have children with 0 counts
              return;

          for(Object child : node.children.keySet()){

              if (node.content != null){

                 int idx = snps.indexOf(node.content);
                 
                 for (int i=snps.size()-1; i>idx; i--){
                     snps.pop(); // removing from end of Stack or ArrayList is very efficient (done in O(1) time)
                 }
              }else{ // if we get back to the root we need to empty whole stack again
                  // in this method, this condition never met, because start node is always after level one
                  snps.clear(); // empty whole stack
              } 

              traverse_partial((TreeNode)child, snps, finalNode, weight);

          } //for 
          
          
       }catch(Exception ex){
         throw new Exception (" Stack size :" + snps.size() + " last allele:" + snps.peek().allele);
       }
    }
    
 
    private List<TreeNode> getPotentialChildrenFlexible(TreeNode node){
        List<TreeNode> out = new ArrayList();
        int node_level = ((SNP)node.content).orderID;
        if (node_level==depth)
            return out;// no potential children
        TreeNode sibiling = getAlternateNode(node);
        if (sibiling != null && sibiling.children.size()>1 ){
            // both nodes at the next level could be  potential children; the logic behind this is if sibling has got 2 children (one of them cause by wrong mapping of reads) so why not me?
            for(TreeNode child : graphNodes.get(node_level+1)){
                out.add(child);
            }

        }else{
            // if there is no sibling, or is a sibling with 0 or 1 childs; 
            // if there is a sibling with 0 children then both nodes at the next level are potential children
            for(TreeNode child : graphNodes.get(node_level+1)){
                if(!child.hasParent)
                    out.add(child);
            }
            
        }
        
        return out;
    }
    
    public int  addAlternativeEdges(){
        int cc = 0;
        for(int i=1; i<depth; i++){
            for(TreeNode node : graphNodes.get(i)){
                TreeNode sibling = getAlternateNode(node);
                if (node.children.size() == 1 && sibling != null){                                     
                    Entry<TreeNode, Integer> entry = (Entry<TreeNode, Integer>)node.children.entrySet().iterator().next(); 
                    if ((int)entry.getValue()>1 ){                    
                       TreeNode nephew = getAlternateNode(entry.getKey());
                       //check if there is an edge from sibling to nephew
                       if (nephew != null && !sibling.children.containsKey(nephew)){
                           // now we add alternative edge
                           sibling.children.put(nephew, entry.getValue());
                           cc++;
                       }

                    }
                }
            }//for each TreeNode at this level
        }//for each level of DAG
        return cc;
    }
    
    
    
    
    
    
    private void traverse_partial_V2(TreeNode node, Stack<SNP> snps, TreeNode finalNode, List<List<SNP>> solutions, List<SNP> deadend_path) throws Exception{
          try{
          
          if (node.content != null){ 
              SNP snp = (SNP)node.content;
              snps.push(snp);
          }
          if(node.content.compareTo(finalNode.content)==0){
              solutions.add(new ArrayList(snps.subList(0, snps.size())));
              return;
          }
          
          if (((SNP)node.content).orderID>=((SNP)finalNode.content).orderID){
              // If we are at the level of finalNode but pointiong to its sibiling, this means we have come to a deadend
              deadend_path.addAll(new ArrayList(snps.subList(0, snps.size())));
                return;
          }
        if (node.children.size()==0){ 
          for(TreeNode pchild : getPotentialChildrenFlexible(node)){
              if (node.content != null){

                 int idx = snps.indexOf(node.content);

                 for (int i=snps.size()-1; i>idx; i--){
                     snps.pop(); // removing from end of Stack or ArrayList is very efficient (done in O(1) time)
                 }
              }
              traverse_partial_V2((TreeNode)pchild, snps, finalNode, solutions, deadend_path);
          }

      }else{
          for(Object child : node.children.keySet()){

              if (node.content != null){

                 int idx = snps.indexOf(node.content);

                 for (int i=snps.size()-1; i>idx; i--){
                     snps.pop(); // removing from end of Stack or ArrayList is very efficient (done in O(1) time)
                 }
              }

              traverse_partial_V2((TreeNode)child, snps, finalNode, solutions, deadend_path);

        } //for 
      }//else
          
       }catch(Exception ex){
         throw new Exception (" Stack size :" + snps.size() + " last allele:" + snps.peek().allele);
       }
    }
    
   
}
