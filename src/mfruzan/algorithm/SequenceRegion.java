/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package mfruzan.algorithm;

import java.util.LinkedList;

/**
 *
 * @author mario.fruzangohar
 */
public class SequenceRegion {
    public String contig = null;
    public long start = 0;  //start is always less than end
    public long end = 0;
    private LinkedList<SequenceRegion> regions = new LinkedList();
    
    public SequenceRegion(){
        
    }
    
    public SequenceRegion(String cntg, long s, long e){
        this.contig = cntg;
        this.start = s;
        this.end = e;        
    }
    
    public boolean isClose(SequenceRegion reg, int max_gap){
        if (!reg.contig.equals(this.contig))
            return false;
        //so we get here if contigs are the same
        // see if there is overlap
        if(reg.start<=this.end && reg.start>=this.start)
            return true;
        if(this.start<=reg.end && this.start>=reg.start)
            return true;
        // now if there is no overlap, see if they are close enough
        if(this.start>reg.end){
            if (this.start-reg.end <=max_gap)
                return true;
        }else if (reg.start>this.end){
            if (reg.start-this.end <=max_gap)
                return true;            
        }
        return false;
    }
    private boolean comesBefore(SequenceRegion reg){
        if(!reg.contig.equals(contig))
            return false;
        if(start < reg.start)
            return true;
        
        return false;
    }
    private boolean comesAfter(SequenceRegion reg){
        if(!reg.contig.equals(contig))
            return false;
        if(start > reg.start)
            return true;
        
        return false;
    } 
    private boolean inside(SequenceRegion reg){
        // returns true if reg is inside current Region
        if(!reg.contig.equals(contig))
            return false;
        if(start >= reg.start && end <= reg.end)
            return true;
        
        return false;
    }
    
    public boolean add(SequenceRegion region){
        // finds proper position in the list for region and insert or replace it in the list
        int i=0;
        for(; i<regions.size(); i++){
            SequenceRegion reg = regions.get(i);
            if (i==0){
                if (region.comesBefore(reg)){
                    regions.add(0, region);
                    return true;
                }
            }else{
                // i>0, then aHSP should be between element i-1 and i to be able to be inserted into i
                if (region.comesBefore(reg) && region.comesAfter(regions.get(i-1))){
                    regions.add(i, region);
                    return true;
                }
            }
        }
        // we get here for 3 reasons
        //1: list is empty
        if (regions.size()==0){
            regions.add(region);
            return true;
        }
        //2: aHSP comes after last element of the list
        if (region.comesAfter(regions.get(regions.size()-1))){
            regions.add(region);// add it 
            return true;
        }
        //3: aHSP has significant overlap with existing HSPs in the list
        
        return false;
    }
    
    public boolean contains(SequenceRegion region){
        // if regions have overlap then this method is not working properly and needs to be re-written
        for(SequenceRegion reg : regions){
            if(region.comesAfter(reg))
                break;
            if(region.inside(reg))
                return true;
        }
        
        return false;
            
    }
}
