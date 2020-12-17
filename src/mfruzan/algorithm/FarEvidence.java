/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mfruzan.algorithm;

/**
 *
 * @author a1195806
 */
public class FarEvidence implements Comparable<FarEvidence>{
    //Used in DAGHaplotype to sort far evidences
    public SNP from = null;
    public SNP to = null;
    public double weight = 0;
     public FarEvidence(SNP f, SNP t, double w){
         from = f;
         to = t;
         weight = w;
     }
    @Override
   
    public int compareTo(FarEvidence anotherFE) {
        if (from.contig.equals(anotherFE.from.contig) && to.contig.equals(anotherFE.to.contig) && from.position==anotherFE.from.position && to.position==anotherFE.to.position)
            return 0;
        if (from.contig.equals(anotherFE.from.contig) && to.contig.equals(anotherFE.to.contig) && to.position<anotherFE.to.position)
              return -1;
        if (from.contig.equals(anotherFE.from.contig) && to.contig.equals(anotherFE.to.contig) &&  to.position==anotherFE.to.position && from.position>anotherFE.from.position )
            return -1;
        
        return 1;
    }
    @Override
   
    public String toString() {
        if (from==null || to==null)
            return "";
        return (String.format("From %d %s To %d %s", from.position, from.allele, to.position, to.allele));
    }
    
    
    
    
}
