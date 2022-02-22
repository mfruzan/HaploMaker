/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package mfruzan.common;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
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
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import mfruzan.algorithm.SNP;

/**
 *
 * @author mario.fruzangohar
 * This class encapsulate output of genotype task
 * each line of file : contig id, position, ref_base, sampe_1 ,......sample_n  
 * each sample genotype can be something like A, T , AD, D, A/+CC and ...
 */
public class GenotypeFileProcessor {
    
    public String loadSNPListFromVCF(String genotypeFile, String ref_seq, String contig, int start_pos, int end_pos, List<SNP> list ){
        // This method  build a list of SNPs from sample that is mapped to reference contig from start_pos up to end_pos (both inclusive), 
        int step = 1;
           
                       
       
        VCFFileReader gfr = null;
        
        try { 
            
            gfr =  new VCFFileReader(new File(genotypeFile), false); 
            // we know that positions in genotypeFile are ordered
            Iterator<VariantContext> it = gfr.iterator();
            
            //go to the begining of the contig in genotype file
            
           step = 2;
           VariantContext vc = null;
            while(it.hasNext()) { 
                vc = it.next();
                if (vc.getContig().equals(contig))
                    break;
            }
            step = 3;
            if(!it.hasNext())
                return null; // there is no entry for this region in vcf file
            //so we get here when we are at the first entry for this contig
            step = 4;
            // here we advance on genotype file until we get to the first position inside the range [start_pos, end_pos]
            while(it.hasNext()){                
                if (!vc.getContig().equals(contig) || vc.getStart() >= start_pos)
                    break;
                
                vc = it.next();
            }
            step = 5;
            // so we get here if contig does not match or snp.position>=start_pos
            if(!vc.getContig().equals(contig)) 
                return null; // there is no entry for this region in genotype file
            
            // so we get here if snp.pos>=start_pos (we are inside the desired region)
            int last_pos = -1; // last index in ref_seq (0-based)
            int snp_order = 1;
            step = 6;
            while (it.hasNext() && vc.getContig().equals(contig) && vc.getStart() <= end_pos){
                // from last_pos to current postion add reference bases to buffer 
                step = 61;

                    step = 62; 
                    String ref_str = vc.getReference().getBaseString().toUpperCase();                 
                    String alleles = Helper.getVCFAllelesStr(vc, 0); // first sample in VCF file has index 0, it returns something like A/AT for het or T for homozygous
                    String alt_str = Helper.getVCFAlternateStr(vc);// returns something like A or A/AT
                    
                    String[] parts = alleles.split("/");
                    step = 63;
                    if(parts.length == 2){
                        // it is het 
                        byte snp_type = 0;                      
                        // now define a SNP 
                        if (ref_str.length()==parts[0].length() && ref_str.length()==parts[1].length()){
                            // it is multibase 
                            snp_type = Helper.SNP_MultiBase;
                        }else if (ref_str.length() < parts[0].length() || ref_str.length() < parts[1].length()){
                            // it is an insertion het position
                            snp_type = Helper.SNP_MultiInsert;
                        }else if (ref_str.length() > parts[0].length() || ref_str.length() > parts[1].length()){
                            // it is a deletion het position
                            snp_type = Helper.SNP_MultiDelete;
                        }

                        SNP nSNP = new SNP(snp_type, snp_order++, vc.getContig(), vc.getStart(), alleles, alt_str);
                        nSNP.ref = ref_str;
                        list.add(nSNP);

                    }
                    last_pos = vc.getEnd() - start_pos;
                //}// if (vc.getStart()- start_pos>last_pos)
                vc = it.next();
            } // while
           step = 7;
            
            
        } catch (Exception ex) {System.out.println("exception inside  loadSNPListFromVCF at step "+ step + " " +ex.getMessage());}
        finally{
            try {               
                gfr.close();
            } catch (Exception ex) {System.out.println("exception closing vcf file" + ex.getMessage() );}
        }    
        
        return null;

    }
    
   
    public void constructAllHaplotypeBlocksVCF(String ref_file, String refIndex, String variantFile, String outputFile, String bamFile, String read_names_file, int insert_size, boolean paired, int min_mapq, int min_hap_len, int max_mapq, boolean highError){
        // This method is good when there is more than one DNA strand or haplotype, mapped to the reference 
        //(In human, where 2 DNA strands are mapped the the same reference, so we use this method. 
        // In wheat if our reference is not heomeolog specific or wheat is heterozygous then we use this method)        
        // bam file is used to extract its original reads for haplotype phasing (paired-end bam file is better, because is longer and can connect 2 consecutive SNPs)
        
        FileWriter writer = null;        
        FastaSequenceIndex fsi = null;
        IndexedFastaSequenceFile ifsi = null;
        SamReader sr = null;
        boolean has_reached = false; // this variable is just used for debugging
        SAMRecordIterator it = null;
        try {
            
            final SamReaderFactory factory =  SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS).validationStringency(ValidationStringency.LENIENT);
            //final SamInputResource resource = SamInputResource.of(new File(bamFile)).index(new File(bamIndex));  
            
            
            sr = factory.open(new File(bamFile));
            it = sr.iterator();
        
            fsi = new FastaSequenceIndex(new File(refIndex));
            ifsi =new IndexedFastaSequenceFile(new File(ref_file) , fsi);
             
            //for each contig in VCF file find first and last het positions 
            
            writer = new FileWriter(outputFile);
            VCFFileReader reader = new VCFFileReader(new File(variantFile), false);
            Iterator<VariantContext> vcit = reader.iterator();           
            
            String last_contig = null;
            int start_pos = 0;
            int last_pos = 0;
            while (vcit.hasNext()){
                VariantContext snp = vcit.next();
                if (!snp.getContig().equals(last_contig)){
                    //contig hast changed
                    if(last_contig !=null){
                        //process it 
                        if(last_pos>start_pos){
                            ReferenceSequence rs = ifsi.getSubsequenceAt(last_contig, start_pos, last_pos);   // length of rs will be = start_idx - end_idx + 1               
                            String ref_frag = new String(rs.getBases(), "UTF-8");
                            System.out.println("----------------------------------------------------------");
                            System.out.println("Reading " + last_contig+ " from " + start_pos + " until " + last_pos);
                            List<SNP> snpList = new ArrayList();  // to be filled in the next command
                            String polymorphic_frag = loadSNPListFromVCF(variantFile, ref_frag, last_contig, start_pos  , last_pos , snpList);
                            System.out.println("Total Number of SNPs:" + snpList.size());
                            if (snpList.size()>0){
                                for (SNP sn : snpList){  
                                    System.out.print(String.format("(%d,%d,%s)", sn.orderID, sn.position, sn.genotype)); 
                                }

                                System.out.println();
                                 // Now for each SNP chain (continous SNPs that dont contain broken point), phase haplotye separately                              
                                 PileupProcessor pp = new PileupProcessor();
                                if (paired){
                                    //try{
                                      pp.phaseGenomeHaplotypesPairedFastVCF_V2(writer, fsi, ifsi, it, last_contig , snpList,  insert_size, min_mapq, min_hap_len);// 155 is maximum insert size
                                    //}catch (Exception ex){System.out.println("Exception outside phaseGenomeHaplotypesPairedFast " + ex.getMessage());}
                                } else 
                                    pp.phaseGenomeHaplotypesSEFast(writer, fsi, ifsi, it, last_contig , snpList,  insert_size, min_mapq, min_hap_len);
                            } // if we have any snp
                            
                        }//if there are at least 2 SNP
                    }
                    start_pos = snp.getStart();
                }
                last_pos = snp.getStart();
                last_contig = snp.getContig();
            }//while there is more line in Variant file
            // now process the last bit 
            if(last_pos>start_pos){
                ReferenceSequence rs = ifsi.getSubsequenceAt(last_contig, start_pos, last_pos);   // length of rs will be = start_idx - end_idx + 1               
                String ref_frag = new String(rs.getBases(), "UTF-8");
                System.out.println("----------------------------------------------------------");
                System.out.println("Reading " + last_contig+ " from " + start_pos + " until " + last_pos);
                List<SNP> snpList = new ArrayList();  // to be filled in the next command
                String polymorphic_frag = loadSNPListFromVCF(variantFile, ref_frag, last_contig, start_pos  , last_pos , snpList);
                System.out.println("Total Number of SNPs:" + snpList.size());
                if (snpList.size()>0){
                    for (SNP sn : snpList){  
                        System.out.print(String.format("(%d,%d,%s)", sn.orderID, sn.position, sn.genotype)); 
                    }

                    System.out.println();
                    PileupProcessor pp = new PileupProcessor();
                    if (paired){
                        //try{
                          pp.phaseGenomeHaplotypesPairedFastVCF_V2(writer, fsi, ifsi, it, last_contig , snpList,  insert_size, min_mapq, min_hap_len);
                        //}catch (Exception ex){System.out.println("Exception outside phaseGenomeHaplotypesPairedFast " + ex.getMessage());}
                    } else 
                        pp.phaseGenomeHaplotypesSEFast(writer, fsi, ifsi, it, last_contig , snpList,  insert_size, min_mapq, min_hap_len);
                        
                } // if we have any snp                
            }
        } catch (Exception ex) {System.out.println("Exception inside constructAllHaplotypeBlocksVCF " + ex.getMessage());}
        finally{
            try {                
                writer.close();
                ifsi.close();
                it.close();
                sr.close();                                
            } catch (IOException ex) {System.out.println("Stream is already closed!!");}
        }    

    }
    
    
}
