/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package mfruzan.common;

import htsjdk.samtools.BAMRecord;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;
import mfruzan.algorithm.FarEvidence;
import mfruzan.algorithm.HaplotypeDAG;
import mfruzan.algorithm.MapUtil;
import mfruzan.algorithm.SNP;
import mfruzan.algorithm.SequencePosition;
import mfruzan.algorithm.SequenceRegion;

/**
 *
 * @author mario.fruzangohar
 */
public class PileupProcessor {
    
   
    
    
   
 
    
    
   
    public void phaseGenomeHaplotypesPairedFastVCF_V2(FileWriter writer, FastaSequenceIndex fsi, IndexedFastaSequenceFile ifsi, SAMRecordIterator it, String ref_contig,  List<SNP> snpList,   int max_insert_length, int min_mapq, int min_hap_len) {
 
       

               Map<String, BAMRecord> map = new HashMap(); 
int step = 1;
        int reads_have_unexpected_SNP_positions = 0;
        int snps_have_unequal_allels = 0;
        int snps_have_unequal_indel_allels = 0;
        int indel_in_critical_regions = 0;
        try{
        int contigLength = (int)fsi.getContigSize(ref_contig);
        if (snpList.size() == 1){
            System.out.println("In whole region there is only one SNP and no need for haplotype.");
            return ;
            
        }

            HaplotypeDAG dag = new HaplotypeDAG(snpList.size()) ;
            int cnt_inserted =0; 
            int cnt_null_allele =0; 
            int cnt_null_read = 0;
            int start_search = 0;
            start_search = snpList.get(0).position - (max_insert_length * 2);
            if (start_search <=0)
                start_search = 1;

            int end_search = 0; 
            end_search = snpList.get(snpList.size()-1).position + (max_insert_length * 2);
            if (end_search > contigLength)
                end_search = contigLength;


            System.out.println("Read fragments from position " + start_search + " to " + end_search);
            
            BAMRecord reco = (BAMRecord)it.next();
            while(it.hasNext() && !reco.getReferenceName().equals(ref_contig))
                    reco = (BAMRecord)it.next();
            if (!it.hasNext()){
                System.out.println("No reads were found in bam file for this region");
                return;
            }
            while(reco.getReferenceName().equals(ref_contig) && reco.getAlignmentStart()<start_search)
                reco = (BAMRecord)it.next();            
            if(!reco.getReferenceName().equals(ref_contig)){
                System.out.println("No reads were found in bam file for this region");
                return;                
            }
            
            int paired_cnt = 0;
            int snp_list_search_idx = 0; //index of snpList that we start search 
            while (reco.getReferenceName().equals(ref_contig) &&  reco.getAlignmentStart() <= end_search){
                if (reco.getMappingQuality()>=min_mapq){
step = -4;                    
                    BAMRecord pair1 = map.get(reco.getReadName());
step = -3;                    
                    if (pair1!=null){
                        paired_cnt++;
                        if (paired_cnt%10000 == 0)
                           System.out.println(paired_cnt/10000 + "0000 paired-end reads were loaded"); 
                        BAMRecord pair2 = reco;
                        int fragment_start = pair1.getAlignmentStart(); 
                        int fragment_end = pair2.getAlignmentEnd();

                        double weight =(double)(pair1.getMappingQuality()+pair2.getMappingQuality())/(40*2);
//---------------------------------------------------------------------------------------
                        List<SNP> subList = Helper.getSNPSubList(snp_list_search_idx, snpList, fragment_start, fragment_end) ;
                       if (subList!= null && subList.size()>=1  /*&& !has_unexpected_SNP_position*/){
                            
                           snp_list_search_idx = subList.get(0).orderID -1;
                           
                           if(snp_list_search_idx>=10)
                               snp_list_search_idx -= 10;
                           else
                               snp_list_search_idx = 0;
                            boolean has_unexpected_allele = false;
                            
                            for(SNP snp : subList){    
                                   boolean allele_confirmed_by_both_reads = false;

                                    String pair1_allele = null;
                                    String pair2_allele = null;
                                    if (snp.position >= pair1.getAlignmentStart() && snp.position<= pair1.getAlignmentEnd()){
                                        int left = pair1.getAlignmentStart()-pair1.getUnclippedStart();
                                        int right = pair1.getUnclippedEnd()-pair1.getAlignmentEnd();
                                        String read = pair1.getReadString();
                                        read = read.substring(left, read.length()-right);
                                        Cigar cigar = pair1.getCigar();
                                        Cigar ncg = new Cigar();
                                        for (CigarElement ce : cigar.getCigarElements()){
                                            if (!ce.getOperator().equals(CigarOperator.S)) 
                                                ncg.add(ce);
                                        }
                                        snp.allele = Helper.getSNPAlleleFromRead(read, ncg, pair1.getAlignmentStart(), snp.position, snp.position + snp.ref.length()-1);
                                        pair1_allele = snp.allele;
                                    }
                                    if (snp.position >= pair2.getAlignmentStart() && snp.position<= pair2.getAlignmentEnd()){
                                        int left = pair2.getAlignmentStart()-pair2.getUnclippedStart();
                                        int right = pair2.getUnclippedEnd()-pair2.getAlignmentEnd();
                                        String read = pair2.getReadString();
                                        read = read.substring(left, read.length()-right);
                                                       
                                        Cigar cigar = pair2.getCigar();
                                        Cigar ncg = new Cigar();
                                        for (CigarElement ce : cigar.getCigarElements()){
                                            if (!ce.getOperator().equals(CigarOperator.S)) 
                                                ncg.add(ce);
                                        }
 
                                        snp.allele = Helper.getSNPAlleleFromRead(read, ncg, pair2.getAlignmentStart(), snp.position, snp.position + snp.ref.length()-1);
                                        pair2_allele = snp.allele;  
                                    }
                                    if (pair1_allele != null && pair2_allele != null){

                                        if (!pair2_allele.equals(pair1_allele)){
                                            snps_have_unequal_allels++;
                                            if(pair1_allele.length() != pair2_allele.length())
                                                snps_have_unequal_indel_allels++;

                                            if(snp.isAlleleValid(pair1_allele) && !snp.isAlleleValid(pair1_allele))
                                                snp.allele = pair1_allele;
                                            else if (!snp.isAlleleValid(pair1_allele) && snp.isAlleleValid(pair1_allele))
                                                snp.allele = pair2_allele;
                                            else{
                                                int pair1_distance_to_middle = Math.abs(pair1.getReadPositionAtReferencePosition(snp.position, true) - (pair1.getReadLength()/2));
                                                int pair2_distance_to_middle = Math.abs(pair2.getReadPositionAtReferencePosition(snp.position, true) - (pair2.getReadLength()/2));
                                                if (pair2_distance_to_middle < pair1_distance_to_middle){
                                                    snp.allele = pair2_allele;
                                                }else{
                                                    snp.allele = pair1_allele;
                                                }
                                            }
                                            System.out.println(snp.toString() + " was resolved by 2 pairs differently, pair 1: " + pair1_allele + " Pair 2:" + pair2_allele + " Finally resolvd by " + snp.allele);

                                        } else{
                                            allele_confirmed_by_both_reads = true;
                                            snp.allele = pair1_allele;
                                        }
                                    }else if (pair1_allele != null){
                                        snp.allele = pair1_allele;
                                    }else if (pair2_allele != null){
                                       snp.allele = pair2_allele; 
                                    }else{
                                        cnt_null_allele++;
                                    }                                

                                if (snp.allele!=null && !snp.validateVCF() ){
                                    snp.allele = null; //unexpected allele
                                    has_unexpected_allele = true;

                                }
                               if( (snp.isVCFDeletion() || snp.isVCFInsertion()) && !allele_confirmed_by_both_reads){
                                    if ((snp.position>=pair1.getAlignmentStart() && snp.position <= pair1.getAlignmentStart() + 2) || (snp.position>=(pair1.getAlignmentEnd() - 2) && snp.position <= pair1.getAlignmentEnd() ) || (snp.position>=pair2.getAlignmentStart() && snp.position <= pair2.getAlignmentStart() + 2) || (snp.position>=(pair2.getAlignmentEnd() - 2) && snp.position <= pair2.getAlignmentEnd() )){
                                        snp.allele = null;
                                        indel_in_critical_regions++;
                                    }
                                    
                                }
                                snp.setVCFAlleleCode();
step = 11;
                            }
                            int null_chunck_start_idx = -1; 
                            int snp_chunck_start_idx = -1;  
                            for (int k=0; k<subList.size() ;k++){
                                if (subList.get(k).allele==null){
                                    if (null_chunck_start_idx==-1){
                                        null_chunck_start_idx = k;
                                    }
                                    if(snp_chunck_start_idx>-1){
                                        if(k - snp_chunck_start_idx>0){
                                            // insert it into dag
                                            dag.insert(subList.subList(snp_chunck_start_idx, k), weight);
                                            dag.insert(dag.getAlternateSNPList(subList.subList(snp_chunck_start_idx, k)), weight);
                                            cnt_inserted++;
                                        }
                                        
                                        snp_chunck_start_idx = -1;
                                    }
                                }else{
                                    if (snp_chunck_start_idx == -1)
                                        snp_chunck_start_idx = k;
                                    if (null_chunck_start_idx>-1){
                                        if (null_chunck_start_idx>0 ){
                                            dag.far_evidence.add(new FarEvidence(subList.get(null_chunck_start_idx-1), subList.get(k), weight)); 
                                                                                        
                                        }
                                        null_chunck_start_idx = -1;
                                    }
                                }
                            }//for k
                            if (snp_chunck_start_idx > -1){
 
                                cnt_inserted++;
                            }
                       }
                       else 
                           cnt_null_read++;
                        
                    }else{
                        map.put(reco.getReadName(), reco);
                    }
                }
                if (it.hasNext())
                   reco = (BAMRecord)it.next();
                else
                    break;
            } 
            
            dag.writeFullLevels("accepted_SNP.bed");
            
            System.out.println("Number of paired reads processed " + paired_cnt);
            System.out.println("Number of single end reads processed " + map.size());
            System.out.println("Number of fragments contain no SNP " + cnt_null_read);
            System.out.println("Number of fragments have SNP at unexpected Positions " + reads_have_unexpected_SNP_positions);
            System.out.println("Number of null alleles " + cnt_null_allele);
            System.out.println("Number of times a SNP was resolved by R1 and R2 differently " + snps_have_unequal_allels);
            System.out.println("Number of times a SNP was resolved by R1 and R2 differently having diffrent length " + snps_have_unequal_indel_allels);
            System.out.println("Number of times an InDel SNP was in critical regions " + indel_in_critical_regions);
            map.clear();
            map = null;
            
            System.out.println(cnt_inserted + " fragments inserted into the Graph");
            if (cnt_inserted==0){
                System.out.println("No Fragments have been added.");
                return ;
            }
            int emptylevel = dag.getEmptyLevel();
            if (emptylevel>0){
                System.out.println("Number of Levels with no node " + emptylevel + " because no read has been mapped, program has missed some reads, you probabily need to change parameters.");
            }

            int onelevel = dag.getOneLevel();
            if (onelevel>0){
                System.out.println("Number of Levels with 1 node " + onelevel  );
            }
            
            if(emptylevel>0 || onelevel>0){
                dag.fillLevels(snpList);
            }

            List<Integer> brokenLevels = dag.getBrokenLevels();
            System.out.println("number of broken levels " + dag.getBrokenPointsNum());
            if(dag.far_evidence.size() > 0){
               System.out.println("Total Far Evidence : " + dag.far_evidence.size()); 
               dag.initFarEvidence();
               int added_far_ev = dag.resolve_far_evidence_Basic(); 
               System.out.println("Added Far Evidence with unique paths : " + added_far_ev); 
               List<Integer> connected_levels = dag.countConnectedLevels(brokenLevels);
               System.out.println("Using broken evidence strategy helped us to connect " + connected_levels.size()  + " broken levels.");
            }
 

            System.out.println("Estimating Haplotypes..." );
            List<SNP> master_hap1 = new ArrayList();
            List<SNP> master_hap2 = new ArrayList();          
            
            step = 5;
            Stack<SNP> hap = dag.findNextHaplotype();
            int debug = 1;
            while(hap != null){
               
 step = 24;                
                if (hap.size()>= min_hap_len){
                    List<SNP> hap1 = new ArrayList();
                    List<SNP> hap2 = new ArrayList();          
                    while(!hap.empty()){                              
                          hap1.add(hap.pop());                              
                    }
                    for(int j=0; j< hap1.size(); j++){   
                        SNP dSNP = hap1.get(j);
                        step = 51;
                        SNP other_SNP = dag.getAlternateSNP(dSNP);
                        hap2.add(other_SNP);
                    }
                    if(master_hap1.size()>0 && master_hap1.get(master_hap1.size()-1).orderID != hap1.get(0).orderID){
                        writer.write(String.format("%s\t%d\t%d\t%d\n",ref_contig, master_hap1.get(0).position, master_hap1.get(master_hap1.size()-1).position, master_hap1.size())); // each block is separated by empty line
                        for(int j=0; j< master_hap1.size(); j++){   
                             SNP dSNP = master_hap1.get(j);
                             SNP pSNP = master_hap2.get(j);
                            step = 53;         
                            writer.write(String.format("%s\t%s\t%d|%d\n", dSNP.printAsFull(), pSNP.allele, dSNP.allele_code, pSNP.allele_code));
                        } 
                        master_hap1.clear();
                        master_hap2.clear();
                        master_hap1.addAll(hap1);
                        master_hap2.addAll(hap2);
                        
                        
                    }else{
                        if(master_hap1.size()>0){
                            if(hap1.get(0).equals(master_hap1.get(master_hap1.size()-1))){
                                master_hap1.addAll(hap1.subList(1, hap1.size()));
                                master_hap2.addAll(hap2.subList(1, hap2.size()));
                            }else{
                                master_hap1.addAll(hap2.subList(1, hap2.size()));
                                master_hap2.addAll(hap1.subList(1, hap1.size()));                            
                            }
                        }else{
                            master_hap1.addAll(hap1);
                            master_hap2.addAll(hap2);                            
                        }
                        
                    }                    
                    
                }
                debug++;
                hap = dag.findNextHaplotype();
            }
            if(master_hap1.size()>0){
                writer.write(String.format("%s\t%d\t%d\t%d\n",ref_contig, master_hap1.get(0).position, master_hap1.get(master_hap1.size()-1).position, master_hap1.size())); // each block is separated by empty line
                for(int j=0; j< master_hap1.size(); j++){   
                     SNP dSNP = master_hap1.get(j);
                     SNP pSNP = master_hap2.get(j);
                    writer.write(String.format("%s\t%s\t%d|%d\n", dSNP.printAsFull(), pSNP.allele, dSNP.allele_code, pSNP.allele_code));
                } 

            }
            System.out.println("Done");
                       
            
        } catch (Exception ex) {System.out.println("Exception inside phaseGenomeHaplotypesPairedFastVCF contig:" + ref_contig + " at step:" + step + ex.getMessage()); }
        
    }
    
    public void phaseGenomeHaplotypesSEFastVCF(FileWriter writer, FastaSequenceIndex fsi, IndexedFastaSequenceFile ifsi, SAMRecordIterator it, String ref_contig,  List<SNP> snpList,   int max_insert_length, int min_mapq, int min_hap_len, Set<String> read_names) {
       
        int step = 1;
        try{
        int indel_in_critical_regions = 0;
        int contigLength = (int)fsi.getContigSize(ref_contig);
        if (snpList.size() == 1){
            // haplotype is already determined!            
            System.out.println("In whole region there is only one SNP and no need for haplotype.");
            return ;
            
        }
        // we get here when there are at least 2 SNP

            HaplotypeDAG dag = new HaplotypeDAG(snpList.size()) ;
            int cnt_inserted =0; // number of inserted fragments into hapDag
            int cnt_null_allele =0; // number of null alleles
            int cnt_null_read = 0; //number of fragments contain no SNP
            int start_search = 0;
            start_search = snpList.get(0).position - (max_insert_length * 2);
            if (start_search <=0)
                start_search = 1;

            int end_search = 0; 
            end_search = snpList.get(snpList.size()-1).position + (max_insert_length * 2);
            if (end_search > contigLength)
                end_search = contigLength;

            System.out.println("Read fragments from position " + start_search + " to " + end_search);
            
            int unexpected_allele_count = 0;  
            BAMRecord reco = (BAMRecord)it.next();
            while(it.hasNext() && !reco.getReferenceName().equals(ref_contig))
                    reco = (BAMRecord)it.next();
            if (!it.hasNext()){
                System.out.println("No reads were found in bam file for this region");
                return;
            }
            while(reco.getReferenceName().equals(ref_contig) && reco.getAlignmentStart()<start_search)
                reco = (BAMRecord)it.next();            
            if(!reco.getReferenceName().equals(ref_contig)){
                System.out.println("No reads were found in bam file for this region");
                return;                
            }
            
            int SE_cnt = 0;
            int snp_list_search_idx = 0; 
            while (reco.getReferenceName().equals(ref_contig) &&  reco.getAlignmentStart() <= end_search){
                
                if (reco.getMappingQuality()>=min_mapq &&  !reco.getReadUnmappedFlag() ){
                    
                        // process the read  
                        SE_cnt++;
                        if (SE_cnt%10000 == 0)
                           System.out.println(SE_cnt/10000 + "0000 single-end reads were loaded"); 
                        int fragment_start = reco.getAlignmentStart();
                        int fragment_end = reco.getAlignmentEnd();
//---------------------------------------------------------------------------------------
//50 for blasr
                        double weight =(double)(reco.getMappingQuality()/50);
                        List<SNP> subList = Helper.getSNPSubList(snp_list_search_idx, snpList, fragment_start, fragment_end) ;
                        
                        
                       if (subList!= null && subList.size()>=1){
                            
                           snp_list_search_idx = subList.get(0).orderID -1;
                           if(snp_list_search_idx>=10)
                               snp_list_search_idx -= 10;
                           else
                               snp_list_search_idx = 0;
                           
                          
                           
                            for(SNP snp : subList){  
                                int left = reco.getAlignmentStart()-reco.getUnclippedStart();
                                int right = reco.getUnclippedEnd()-reco.getAlignmentEnd();
                                String read = reco.getReadString();
                                read = read.substring(left, read.length()-right);
                                               
                                Cigar cigar = reco.getCigar();
                                Cigar ncg = new Cigar();
                                for (CigarElement ce : cigar.getCigarElements()){
                                    if (!ce.getOperator().equals(CigarOperator.S))
                                        ncg.add(ce);
                                }          
                                snp.allele = Helper.getSNPAlleleFromRead(read, ncg, reco.getAlignmentStart(), snp.position, snp.position + snp.ref.length()-1);

                                if (snp.allele!=null && !snp.validateVCF() ){
                                    snp.allele = null; //unexpected allele                                    
                                    unexpected_allele_count++;
 

                                }
                               if( (snp.isVCFDeletion() || snp.isVCFInsertion()) ){
                                    
                                    if ((snp.position>=reco.getAlignmentStart() && snp.position <= reco.getAlignmentStart() + 2) || (snp.position>=(reco.getAlignmentEnd() - 2) && snp.position <= reco.getAlignmentEnd() ) ){
                                        snp.allele = null;
                                        indel_in_critical_regions++;
                                    }
                                    
                                }                                
                                snp.setVCFAlleleCode();

                            }
 
                            int null_chunck_start_idx = -1; 
                            int snp_chunck_start_idx = -1;  
                            for (int k=0; k<subList.size() ;k++){
                                if (subList.get(k).allele==null){
                                    if (null_chunck_start_idx==-1){
                                        null_chunck_start_idx = k;
                                    }
                                    if(snp_chunck_start_idx>-1){
                                        //process snp chunck
                                        if(k - snp_chunck_start_idx>0){
                                            
                                            dag.insert(subList.subList(snp_chunck_start_idx, k), weight);
                                            dag.insert(dag.getAlternateSNPList(subList.subList(snp_chunck_start_idx, k)), weight);
                                            cnt_inserted++;
                                        }
                                        
                                        
                                        snp_chunck_start_idx = -1;
                                    }
                                }else{
                                    if (snp_chunck_start_idx == -1)
                                        snp_chunck_start_idx = k;
                                    if (null_chunck_start_idx>-1){
                                        if (null_chunck_start_idx>0 ){
                                            dag.far_evidence.add(new FarEvidence(subList.get(null_chunck_start_idx-1), subList.get(k), weight)); // not a good idea to discount weights of far evidence
                                               
                                        }
                                       
                                        null_chunck_start_idx = -1;
                                    }
                                }
                            }//for k


                           
                           
                       }
                       else 
                           cnt_null_read++;

                }
                
                if (it.hasNext())
                   reco = (BAMRecord)it.next();
                else
                    break;
            } 
            
            System.out.println("Number of single-end reads processed " + SE_cnt);
            System.out.println("Number of fragments contain no SNP " + cnt_null_read);
            System.out.println("Number of indels in critical regions  " + indel_in_critical_regions);
             System.out.println("Number of unxpected alleles  " + unexpected_allele_count);
           

            System.out.println(cnt_inserted + " reads inserted into the Graph");
            if (cnt_inserted==0){
                System.out.println("No read have been added.");
                return ;
            }
            
           
            int emptylevel = dag.getEmptyLevel();
            if (emptylevel>0){
           
                System.out.println("Number of Levels with no node " + emptylevel + " because no read has been mapped, program has missed some reads, you probably need to change parameters.");
               
            }

            int onelevel = dag.getOneLevel();
            if (onelevel>0){
                System.out.println("Number of Levels with 1 node " + onelevel  );
            }

            System.out.println("number of broken levels " + dag.getBrokenPointsNum());

            if(emptylevel>0 || onelevel>0){
                dag.fillLevels(snpList);
            }            

            List<Integer> brokenLevels = dag.getBrokenLevels();
            
            
            if(dag.far_evidence.size() > 0){
               System.out.println("Total Far Evidence : " + dag.far_evidence.size());
               dag.initFarEvidence();
               int added_far_ev = dag.resolve_far_evidence_Basic(); 
             
               List<Integer> connected_levels = dag.countConnectedLevels(brokenLevels);
               System.out.println("Using broken evidence strategy helped us to connect " + connected_levels.size()  + " broken levels.");
               
            }// if we have any far evidence

         

        
            System.out.println("Estimating Haplotypes..." );
            List<SNP> master_hap1 = new ArrayList();
            List<SNP> master_hap2 = new ArrayList();          
            
            step = 5;
            Stack<SNP> hap = dag.findNextHaplotype();
            int debug = 1;
            while(hap != null){
               
 step = 24;                
                if (hap.size()>= min_hap_len){

step = 25;
                    List<SNP> hap1 = new ArrayList();
                    List<SNP> hap2 = new ArrayList();          
                   
  
                    while(!hap.empty()){                              
                          hap1.add(hap.pop());                              
                    }
                    for(int j=0; j< hap1.size(); j++){   
                        SNP dSNP = hap1.get(j);
                        step = 51;
                        SNP other_SNP = dag.getAlternateSNP(dSNP);
                        hap2.add(other_SNP);
                    }
                    if(master_hap1.size()>0 && master_hap1.get(master_hap1.size()-1).orderID != hap1.get(0).orderID){
                      
                        
                        writer.write(String.format("%s\t%d\t%d\t%d\n",ref_contig, master_hap1.get(0).position, master_hap1.get(master_hap1.size()-1).position, master_hap1.size())); // each block is separated by empty line
 
                        for(int j=0; j< master_hap1.size(); j++){   
                             SNP dSNP = master_hap1.get(j);
                             SNP pSNP = master_hap2.get(j);
                            step = 53;    
               
                            writer.write(String.format("%s\t%s\t%d|%d\n", dSNP.printAsFull(), pSNP.allele, dSNP.allele_code, pSNP.allele_code));
 
                        } 
                        
                        //empty master haplotypes and new one into it
                        master_hap1.clear();
                        master_hap2.clear();
                        master_hap1.addAll(hap1);
                        master_hap2.addAll(hap2);
                        
                        
                    }else{
                      
                        if(master_hap1.size()>0){
                            if(hap1.get(0).equals(master_hap1.get(master_hap1.size()-1))){
                               
                                master_hap1.addAll(hap1.subList(1, hap1.size()));
                                master_hap2.addAll(hap2.subList(1, hap2.size()));
                            }else{
                                master_hap1.addAll(hap2.subList(1, hap2.size()));
                                master_hap2.addAll(hap1.subList(1, hap1.size()));                            
                            }
                        }else{
                            master_hap1.addAll(hap1);
                            master_hap2.addAll(hap2);                            
                        }
                        
                    }                    
                    
step = 26;

                     
                            
                }//if haplotype has minimum length
step = 27;

                hap = dag.findNextHaplotype();
step = 28;                
            }   
            if(master_hap1.size()>0){
                writer.write(String.format("%s\t%d\t%d\t%d\n",ref_contig, master_hap1.get(0).position, master_hap1.get(master_hap1.size()-1).position, master_hap1.size())); // each block is separated by empty line
                // then rows of haplotype
                for(int j=0; j< master_hap1.size(); j++){   
                     SNP dSNP = master_hap1.get(j);
                     SNP pSNP = master_hap2.get(j);

  
                    writer.write(String.format("%s\t%s\t%d|%d\n", dSNP.printAsFull(), pSNP.allele, dSNP.allele_code, pSNP.allele_code));

                } 

            }
            
            System.out.println("Done");
                       
            
        } catch (Exception ex) {System.out.println("Exception inside phaseGenomeHaplotypesSEFastVCF contig:" + ref_contig + " at step:" + step + ex.getMessage()); }
        
    }
    
    
}
