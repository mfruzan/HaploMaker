/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package mfruzan.algorithm;


import mfruzan.common.Helper;

import java.util.ArrayList;
import java.util.List;
import org.apache.commons.lang.builder.EqualsBuilder;
import org.apache.commons.lang.builder.HashCodeBuilder;



/**
 *
 * @author mario.fruzangohar
 */
public class SNP implements Comparable<SNP> , java.io.Serializable{
    
    public byte type = Helper.SNP_MultiBase;  
    public String ref = null; 
    public String vcf_alt = null; 
    public String contig = null;
    public int position = 0; 
    //public int position_end = 0; // only used in vcf
    public int orderID = 0; 
    public String genotype = null; 
    public String allele = null;
    public byte allele_code = -1; 
    public boolean broken = false; 
    public boolean deleted = false; 
    
    public boolean flagged = false; 
    public boolean isMutation = false; 
    //public boolean isReference = false; //this SNP carries reference allele (is used in CVF file format)
    List<Byte> mapq = new ArrayList();
        
    public SNP(){}
    
    public SNP(byte tp, int order, int pos){
        type = tp;
        orderID = order;
        position = pos;
        
    }
   public SNP(byte tp, int order, int pos, String gt){
        type = tp;
        orderID = order;
        position = pos;
        genotype = gt;
        
    }
   public SNP(byte tp, int order, String cntg, int pos, String gt){
        
        type = tp;
        orderID = order;
        contig = cntg;
        position = pos;
        genotype = gt;
        
    }
   public SNP(byte tp, int order, String cntg, int pos, String gt, String alt){
        
        type = tp;
        orderID = order;
        contig = cntg;
        position = pos;
        genotype = gt;
        vcf_alt = alt;
        
    }   
    public SNP(byte tp, int order, String allele){
        type = tp;
        orderID = order;        
        this.allele = allele;
        
    }
    
    public SNP clone(){
        SNP aSNP = new SNP();
        aSNP.type = this.type;
        aSNP.orderID = this.orderID;
        aSNP.contig = this.contig;
        aSNP.position = this.position;
        aSNP.ref = this.ref;
        aSNP.vcf_alt = this.vcf_alt;
        aSNP.genotype = this.genotype;
        aSNP.allele = this.allele;
        aSNP.allele_code = this.allele_code;
        //aSNP.isReference = this.isReference;
        aSNP.deleted = this.deleted;
        aSNP.broken = this.broken;
        aSNP.flagged = this.flagged;
        aSNP.isMutation = this.isMutation;
        return aSNP;
    }
    public String[] getAllVariaons(){
        String[] arr = null;
        if (type==Helper.SNP_MultiBase){
            arr = new String[genotype.length()];
            for (int i =0; i<genotype.length(); i++)
                arr[i] = genotype.substring(i, i+1);                
        }else if (type==Helper.SNP_MultiInsert){
            arr = genotype.split("/");
        }
        
        return arr;
    }
    
    public boolean hasAlleleLength(int len){
        String[] alleles = genotype.split("/");
        for(String al : alleles){
            if (al.length() == len){
                return true;
            }
        }
            
        return false;
                
    }
    public int maxAlleleLength(){
       int max = 0;
        String[] alleles = genotype.split("/");
        for(String al : alleles){
            if (al.length() > max){
                max = al.length();
            }
        }
            
        return max;                
    }
    public boolean isReference(){
        if (type == Helper.SNP_MultiInsert){
            if (allele.equals("+"))
                    return true;
        }else if (type == Helper.SNP_MultiBase){                       
            if (allele.equals(ref))
                return true;
        }
        return false;
    }
    public boolean isVCFReference(){
        if (allele.equals(ref))
            return true;
        return false;
    }
    public boolean isVCFDeletion(){
           if(allele==null || ref==null)
               return false;
           if(allele.length() < ref.length())
               return true;
           
           return false;
    }
    public boolean isVCFInsertion(){
           if(allele==null || ref==null)
               return false;
           if(allele.length() > ref.length())
               return true;
           
           return false;
    }

    public boolean validate(){
        if (genotype == null)
            return true;
        if (type == Helper.SNP_MultiBase && genotype.indexOf(allele)>=0)
            return true;
        if (type==Helper.SNP_MultiInsert  ){
            for(String ins : genotype.split("/"))
                if (ins.equals(allele))
                    return true;
        }
        return false;
    }
    public boolean isAlleleValid(String alll){
        if (genotype == null)
            return true;
        String[] als = genotype.split("/");
        for(String al : als)
            if (al.equals(alll))
                return true;
        if (als.length == 1)
            return false;
        if (alll.indexOf('N')>=0){
            if (alll.length()== als[0].length() && alll.length()!= als[1].length()){
                return true;
            }else if (alll.length()== als[1].length() && alll.length()!= als[0].length()){
                return true;                
            }
        }
        return false;
        
    }
    public boolean validateVCF(){
        if (genotype == null)
            return true;
        
        String[] als = genotype.split("/");
        for(String al : als)
            if (al.equals(allele))
                return true;
        if (als.length == 1)
            return false;
        if (allele.indexOf('N')>=0){
            if (allele.length()== als[0].length() && allele.length()!= als[1].length()){
                allele = als[0];
                return true;
            }else if (allele.length()== als[1].length() && allele.length()!= als[0].length()){
                allele = als[1];
                return true;                
            }
        }
        return false;
    }
    public void setVCFAlleleCode(){
        if (allele==null)
            return;
       
       if (allele.equals(ref))        
           allele_code = 0;
       if (vcf_alt == null)
           return;
       String[] alts = vcf_alt.split("/");
       for(int i=0; i<alts.length; i++)
           if (allele.equals(alts[i]))
               allele_code =(byte)(i + 1);
    }
    
    public String printAsFull(){
        return String.format("%s\t%d\t%s\t%s", contig, position, ref, allele);
    }
    
    public String printAsTriplet(){
        return String.format("(%d,%d,%s,%s)", orderID, position, genotype, allele);
    }
    
    @Override
    public int compareTo(SNP aSNP) {
        if (aSNP.type==type && aSNP.position==position && aSNP.orderID==orderID && aSNP.allele.equals(allele) )
            return 0;
       return -1;
    }
    
    @Override
    public String toString() { 
        return "SNP at position:"+position+" orderID:"+orderID+" genotype:"+genotype;
    }
    
    @Override
    public int hashCode() {
        return new HashCodeBuilder(23, 67). 
            append(type).
            append(position).
            append(orderID).            
            append(allele).            
            toHashCode();
    }
    @Override
    public boolean equals(Object other){
        if (other == null) return false;
        if (other == this) return true;
        if (!(other instanceof SNP))return false;
        SNP snp = (SNP)other;
        return new EqualsBuilder().append(type, snp.type).append(position, snp.position).append(orderID, snp.orderID).append(allele, snp.allele).isEquals();
        
    }
    
}
