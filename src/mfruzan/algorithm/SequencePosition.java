/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package mfruzan.algorithm;

/**
 *
 * @author mario.fruzangohar
 * This class encapsulates one position of a sequence
 */
public class SequencePosition implements Comparable<SequencePosition>{
    public String contig = null;
    public int pos = 0;
    public int count = 0; // This is number of reads/evidences that confirm this position 
    
    //The following variable is used for multi-sample situation, we want to know which samples have this position set in them
    public boolean[] samples;
    
    public SequencePosition(String cont, int posss){
        this.pos = posss;
        this.contig = cont;
    }
    
    public SequencePosition(String cont, int posss, int samples_num){
        this.pos = posss;
        this.contig = cont;
        this.samples = new boolean[samples_num];
        //initialize the boolean array to false
        for (int i=0; i<samples.length; i++)
            this.samples[i]=false;
    }
    
    @Override
    public int compareTo(SequencePosition anotherSequence) {
        if(contig.compareTo(anotherSequence.contig)!=0)
            return contig.compareTo(anotherSequence.contig);
        //so we get here if 2 positions have the same contig
        if (pos<anotherSequence.pos)
            return -1;
        else if (pos<anotherSequence.pos)
            return 1;
        else
            return 0;
    }
    //Fuzzy comparison
    public int compareWith(SequencePosition seq, int radious){
        if(contig.compareTo(seq.contig)<0)
            return -1;
        if (contig.compareTo(seq.contig)>0)
            return 1;
        // it returns 0 if current seq can be joined to the seq, -1 if current seq comes before seq and +1 if current seq comes after
        if (Math.abs(pos-seq.pos)<=radious)
            return 0;
        if (pos>seq.pos)
            return 1;
        // we get here when current seq comes before seq
        return -1;
    }
    
    public int getSetSamples(){
        int cnt = 0;
        for(int i=0;i<samples.length;i++)
            if (samples[i])
                cnt++;
        
        return cnt;
    }
}
