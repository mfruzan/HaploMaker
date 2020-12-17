/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package mfruzan.common;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeMap;
import mfruzan.algorithm.SNP;
import static mfruzan.common.Helper.T;
import htsjdk.variant.variantcontext.VariantContext;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import mfruzan.algorithm.MapUtil;
import org.biojava3.core.sequence.AccessionID;
import org.biojava3.core.sequence.ChromosomeSequence;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.GeneSequence;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.Strand;
import org.biojava3.core.sequence.TranscriptSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;

/**
 *
 
 */
public class Helper {

    // type of mutation
   public static final byte MT_DELETION = 1; // base  deletion
   public static final byte MT_INSERTION = 2; // base  insertion
   public static final byte MT_Substitution = 3; // base  change

    // type of Variation
   public static final byte VR_DELETION = 1; // base  deletion
   public static final byte VR_INSERTION = 2; // base  insertion
   public static final byte VR_SNP = 3; // base  change
   public static final byte VR_Default = 4; // base  change
   
   // type of SNP or polymorphic position
   public static final byte SNP_MultiBase = 1; // different bases, or Deletion
   public static final byte SNP_MultiInsert = 2; // different insertions, or presence or absence of an Insetion
   public static final byte SNP_MultiDelete = 3; // different insertions, or presence or absence of an Insetion
   
   //relative positions of 2 sequences by each other
   public static final byte POS_Q_Extended_S = 1; // Query sequence is extended by the Subject one from right hand side
   public static final byte POS_S_Extended_Q = 2; // one sequence is extended by the second one from right hand side
   public static final byte POS_Q_Inside_S = 3; // one seq is inside the other one
   public static final byte POS_S_Inside_Q = 4; // one seq is inside the other one
   public static final byte POS_Distinct = 5; // 2 sequences have no overlap

   public static final byte MT_Default = 4; // default mutation type
   public static final byte MT_Hidden = 5; // This is subtype of Transition mutation, when ratio test will reveal mutation
  
   
   // number neucleotie letters alphabetically
   public static final byte A = 0;   
   public static final byte C = 1;
   public static final byte G = 2;
   public static final byte T = 3;
   public static final byte GAP = 5; // in pairwise and mutiple alignment this represent a GAP in the alignment
   public static final byte N = 6; // any base
   public static final byte R = 7;
   public static final byte Y = 8;
   public static final byte S = 9;
   public static final byte W = 10;
   public static final byte K = 11;
   public static final byte M = 12;
   public static final byte B = 13;
   public static final byte d = 14;
   public static final byte H = 15;
   public static final byte V = 16;
   public static final byte UNKNOWN = -1;
   public static final byte D = -2; // base is deleted
   
   public static final double PVALUE_NORMAL = 0.05;
   public static final double PVALUE_STRICT = 0.15;
   
   public static final char[] TRIPLET_VAR = {'!','"','#','$','%','&','\'', '(',')','*','+',',','-','.','/','0','1','2','3','4','5','6','7','8','9', ':',';', '<','=','>','?','@','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','[','\\' ,']', '^','_', '`', 'a','b'}; // 64+1 characters   
   public static final char UPPER_BOUND = '~'; // ascii of this character is greator than the ones in TRIPLET_VAR
   
   public static final HashMap<String, String> codonMap = new HashMap();
   static{
       codonMap.put("TTT", "Phe");
   }
   
   public static int getMaximumIndex(int[] arr){
       // returns the index of the element in the array that has maximum number
       int idx = 0;
       int max = Integer.MIN_VALUE;
       for (int i=0; i<arr.length; i++){
           if (arr[i] > max){
               max=arr[i];
               idx = i;
           }
       }
       
       return idx;
   }
   
   public static byte nt_to_i(char nt){
           switch (nt){
                case 'A':
                case 'a':
                    return A;                       
                case 'T':                        
                case 't':
                case 'U':                        
                case 'u':                    
                    return T;                        
                case 'C':
                case 'c':
                    return C;                        
                case 'G':                        
                case 'g':
                    return G;   
                case '-':
                    return GAP;
                case 'n':
                case 'N':
                    return N;
                case 'd':
                case 'D':
                case '*':
                    return D; 
                case 'R':
                case 'r':
                    return A;
                case 'W':
                case 'w':
                    return A;
                case 'M':
                case 'm':
                    return A;
                case 'Y':
                case 'y':
                    return C;
                case 'S':
                case 's':
                    return C;
                case 'K':
                case 'k':
                    return G;
                case 'B':
                case 'b':
                    return C;
                case 'H':
                case 'h':
                    return A;
                case 'V':
                case 'v':
                    return A;
                    
                    
            }
           return UNKNOWN;       
    }  
   public static String relative_pos_string(byte rpos){
       switch(rpos){
           case POS_Q_Extended_S:
               return "Q extended by S";
           case POS_S_Extended_Q:
               return "S extended by Q";
           case  POS_Q_Inside_S:
               return "Q inside S";
           case POS_S_Inside_Q:
               return "S inside Q";
           case POS_Distinct:
               return "Distinct";
       }
       
       return "UNKNOWN";
   }
   public static int getTripletIndex(char triplet){
       for(int i=0; i<TRIPLET_VAR.length; i++)
           if (TRIPLET_VAR[i]== triplet)
               return i;
       
       return -1;
   }
   public static boolean isGood(String seq){
       // a seq is considered good if it is not homopolymere or contains any n or N
       if (seq.indexOf('n')>=0 || seq.indexOf('N')>=0)
           return false;
       //if (seq.matches("[Aa]+|[Tt]+|[Cc]+|[Gg]+"))
              // return false;
       return true;
   }
    public static int seqNNum(String seq){
       // accepts a seq and counts number of N and return 
       int cnt = 0;
       for (int i=0; i<seq.length(); i++){
           if(seq.charAt(i)=='N' ||  seq.charAt(i)=='n')
               cnt++;
       }
       
       return (cnt);
   }

   public static String getReverse(String seq){
       StringBuffer buf = new StringBuffer();
       for(int i=seq.length()-1; i>=0; i--)
           buf.append(seq.charAt(i));
       
       return buf.toString();
   }
   public static String getRC(String seq){
       // return reverse complement of a sequence (this method guarantees that the length of input parameter is the same as output) 
       StringBuffer buff = new StringBuffer();
       for(int i=seq.length()-1; i>=0; i--){
           char ch = seq.charAt(i);
           byte b = nt_to_i(ch);
           if (b==A)
               buff.append('T');
           else if (b==T)
               buff.append('A');
           else if (b==C)
               buff.append('G');
           else if (b==G)
               buff.append('C');
           else if (b==N)
               buff.append('N');
           else if (b==GAP)
               buff.append('-');           
           else if (b==UNKNOWN)
               buff.append('N');
       }
       return buff.toString();
   }
   public static String standardize(String seq){
       // return lower case and IUPAC code with 
       StringBuffer buff = new StringBuffer();
       for(int i=0; i<seq.length(); i++){
           char ch = seq.charAt(i);
           byte b = nt_to_i(ch);
           if (b==A)
               buff.append('A');
           else if (b==T)
               buff.append('T');
           else if (b==C)
               buff.append('C');
           else if (b==G)
               buff.append('G');
           else if (b==N)
               buff.append('N');
           else if (b==GAP)
               buff.append('-');           
           else if (b==UNKNOWN)
               buff.append('N');
       }
       return buff.toString();
   }
   
   
   public static String simplify(String seq){
       // remove gaps , replace any non-actgn characters with n
        StringBuffer buff = new StringBuffer();
        for(int i=0; i<seq.length(); i++){
           char ch = seq.charAt(i);
           byte b = nt_to_i(ch);
           if (b==A)
               buff.append('A');
           else if (b==T)
               buff.append('T');
           else if (b==C)
               buff.append('C');
           else if (b==G)
               buff.append('G');
           else if (b==N)
               buff.append('N');
           else if (b==UNKNOWN)
               buff.append('N');
            // if it is * D or - do nothing
        }
       return buff.toString();
   }
   public static String extractGeneID(String column9){
       // from column nine of a gff file
       String[] arr = column9.split(";");
       for (String el : arr){
           //if (el.startsWith("gene_id="))
           if (el.startsWith("ID="))
               return el.substring(3);
       }
       return null;
   }
   // another implementation using biojava
   public static String getReverseComplement(String seq){
        
            DNASequence dna = new DNASequence(seq);
            
            return dna.getReverseComplement().getSequenceAsString(); 
   }
   public static double getSimilarity(String seq1, String seq2){
        
        int cc = 0;
        for(int i=0; i<seq1.length(); i++)
            if (seq1.charAt(i)==seq2.charAt(i) || seq1.charAt(i)=='N' || seq2.charAt(i)=='N')
                cc++;
            
            return (double)cc/seq1.length();
   }
   
   
   public static <T> T[] subarray(T[] arr, Class<T> cl, int start_idx, int end_idx){
      //start_idx is inclusive and end_idx is exclusive
       //Class<T> clazz;
       T[] out = (T[]) Array.newInstance(cl, end_idx - start_idx);
       int j = 0;
       for(int i=start_idx; i<end_idx; i++)
           out[j++] = arr[i];
       
       return out;
  }
   
   public static String fixLength(String in, char ch, int len) throws Exception{
       // this methods , adds ch to in if in length is less than len
       if (in.length()>len)
           throw new Exception("Exception inside fixLength");
       if (in.length()== len)
           return in;
       StringBuffer out = new StringBuffer(in);
       for(int i=0; i<len-in.length(); i++)
           out.append(ch);
       
       return out.toString();
   }
  
   public static int hammingDistance(String str1, String str2){
       //find the hamming distance and return it, we assume str1 and str2 have the same length
        if (str1.length() != str2.length())
           return -1;
       
        int dist = 0;
        for(int i=0; i<str1.length(); i++){
           char ch1 = str1.charAt(i);
           char ch2 = str2.charAt(i);
           if(ch1!=ch2 && ch1!='N' && ch2!='N')
               dist++; //so for example if at least one of chars is not N then distance will not increase
        }
       
        return dist;
       
   }
   
   
  
    
 
   
   public static int getGaps(String alignment){
       // it returns number of - exist
       int gaps = 0;
       for (int i=0; i<alignment.length(); i++){
           if (alignment.charAt(i)=='-')
               gaps++;
       }
       return gaps;
   }
   public static TreeMap getGapsMap(String seq){
       // this will return position of the next characer to one or multiple - , for example if sequence is AA--TC-AAT-G then will return 4->2, 7->1, 11->1
       TreeMap<Integer, Integer> out = new TreeMap();
       Pattern p = Pattern.compile("\\-+"); // possessive(+)
        Matcher   m = p.matcher(seq); // remember that qseq contains Gaps or -, but we only look for continuous base without gap in the middle
        while(m.find()){
            // extract the equal part from reference (subject)
            out.put(m.end(), m.end()-m.start());
        } 
        
        return out;
       
   }
   
   
   public static String getVCFAllelesStr(VariantContext vc, int sample_index) throws Exception{
       if (vc.getGenotypes().get(sample_index).getAlleles().size() == 1)// its haploid sucha ax chrX and chrY
           return vc.getGenotypes().get(sample_index).getAllele(0).getBaseString().toUpperCase();
           //throw new Exception("Genotype/Specie is not diploid at contig " + vc.getContig() + ":" + vc.getStart());
       if (vc.getGenotypes().get(sample_index).getAllele(0).getBaseString().equals(vc.getGenotypes().get(sample_index).getAllele(1).getBaseString()))
           return vc.getGenotypes().get(sample_index).getAllele(0).getBaseString().toUpperCase();// homozygous
       else
           return (vc.getGenotypes().get(sample_index).getAllele(0).getBaseString().toUpperCase()+ "/" +vc.getGenotypes().get(sample_index).getAllele(1).getBaseString().toUpperCase());// heterozygous
   }
   public static String getVCFAlternateStr(VariantContext vc) throws Exception{
       StringBuffer buff = new StringBuffer();
       for(htsjdk.variant.variantcontext.Allele al : vc.getAlternateAlleles()){
           buff.append(al.getBaseString().toUpperCase()+"/");
       }
       if (buff.charAt(buff.length()-1)=='/')
           buff.deleteCharAt(buff.length()-1);
       return buff.toString(); 
   }
   
   public static byte ZG_String_Byte(String zg){
       if(zg.toLowerCase().startsWith("hom"))
           return ZG_HOMOZYGOUS;
       if (zg.toLowerCase().startsWith("het"))
           return ZG_HETROZYGOUS;
       return ZG_UNKNOWN;
   }
   public static String Base_String(byte b){
       switch(b){
           case A:
               return "A";
           case T:
               return "T";
           case C:
               return "C";
           case G:
               return "G";
            case N:
               return "N";
            case GAP:
                return "-";
            case D:
                return "D";
            case UNKNOWN:
                return "";
            case R:
                return "A";
            case Y:
                return "C";
            case S:
                return "C";
            case W:
                return "A";
            case K:
                return "G";
            case M:
                return "A";
            case B:
                return "C";
            case H:
                return "AH";
            case V:
                return "A";
                
       }
       
       return " ";
   }
  
   
   public static int scale_to(int scale_to, int cov, int max_cov){
      if (cov==0)
           return 0;

       int new_cov = 0; 
       if(scale_to>1){
           if (cov>=max_cov)
               new_cov = scale_to;
           else{
             int interval = max_cov/scale_to; // it return an integer 
             new_cov = cov/interval + 1;
           }
       }else if (scale_to==1){
           // in this mode we have binary mode, if cov is above max_cov it returns 1 otherwise it is 0
           if(cov>max_cov)
               new_cov = 1;
           else
               new_cov = 0;
       }
       
       
            
       return new_cov;
   }
   public static List<Byte> getProfileColumn(List<String> alignment, int col){
      List<Byte> out = new ArrayList() ;
      for(String alignedseq : alignment)
          out.add(nt_to_i(alignedseq.charAt(col)) );
      
      return out;
   }
   
   public static double phredToProb(int phred){
       return Math.pow(10, (double)phred/-10);
   }
   
   public static double getMAPQ(char ch){
       int ascii = ch; // ch is a character that represent quality of the map
       int phred = ascii - 33;
       return phredToProb(phred);
   }
   public static double getPhred(char ch){
       int ascii = ch; // ch is a character that represent quality of the map
       int phred = ascii - 33;
       return phred;
   }
   public static double getMAPQScale(char ch, int mean_mapq){
       int ascii = ch; // ch is a character that represent quality of the map
       int phred = ascii - 33;
       
       if (phred >= mean_mapq)
           return 1;
       
       double d = (1-phredToProb(phred))/(1-phredToProb(mean_mapq));
       return d;
   }
   
   public static int goToNextInteger(int idx, String[] arr){
       //from index position in the arr sacns the array to get to the first Integer, returns -1 if nothing has been found
       int i = 0;
       for (i=idx; i< arr.length; i++)
           if (arr[i].trim().matches("\\d+"))
               break;
       
       
       if (i==arr.length)
           return -1;
       
       return i;
   }
   public static int getReferencePosition(String polymorphicSequence, int idx, int start_pos){
       //returns reference position of the base at idx (zero based) position in the sequence: AACCDDNNNNNNNNC[+A/AA/AAA]C:ACTT[CD][AC]AADGCNNNNNNNNN
        boolean openBracket = false;
        int pos = start_pos;
        for (int i=0; i<polymorphicSequence.length();i++){
           if (i == idx)
               return pos; 
           char ch = polymorphicSequence.charAt(i);
           switch (ch){
               case '[':
                   openBracket = true;
                   if (polymorphicSequence.charAt(i+1) != '+')
                       pos++;
                   break;
               case ']':
                   openBracket = false;
                   break;
               case 'N':
               case 'A':
               case 'T':
               case 'C':
               case 'G':   
               case 'D':
                   if (!openBracket)
                       pos++;
           }//switch
       } //for
       return -1;
   }
   
   public static int getCharacterCount(char ch, String str){
       // counts number of characters in the string
       int cc = 0;
       for (int i=0; i<str.length(); i++)
           if(str.charAt(i)==ch)
               cc++;
       
       return cc;
   }
   public static int getStartRepeat(String str,int min_len ){
       // if the same character repeats at least min_len time in the begining of the string it will return length of repeat otherwise retunrs 0
       if (str.isEmpty())
           return 0;
       
       int cc = 1;
       char start_ch = str.charAt(0);
       for (int i=1; i<str.length(); i++)
           if(str.charAt(i)==start_ch)
               cc++;
           else 
               break;
       if(cc>= min_len)
         return cc;
       else
           return 0;
   }
   public static int getEndRepeat(String str,int min_len ){
       // if the same character repeats at least min_len time at the end of the string it will return length of repeat otherwise retunrs 0
       if (str.isEmpty())
           return 0;
       
       int cc = 1;
       char end_ch = str.charAt(str.length()-1);
       for (int i=str.length()-2; i>=0; i--)
           if(str.charAt(i)==end_ch)
               cc++;
           else 
               break;
       if(cc>= min_len)
         return cc;
       else
           return 0;
   }
   
   public static String removeLeadingTrailingN(String str){
       // if a sequence is like NNACCCTTNNNNN it will return ACCCTT
       
       int i = 0; //left index 
       for ( ; i<str.length(); i++){
           if (str.charAt(i)!='N')
               break;
       }
       //so when we get here i points to first non-N character
       int j = str.length()-1; //right index 
       for ( ; j>=0; j--){
           if (str.charAt(j)!='N')
               break;
       }
       //j points to first non N character from right hand side
       if (i<=j)
           return str.substring(i, j+1);       
       else
           return "";
   }
   public static Integer sumarr(Collection<Integer> arr){
       int sum = 0;
       for(int n : arr)
           sum += n; 
       return sum;
   }
   
   public static int sequenceDistance(String seq1, String seq2){       
       //this method returns sequence distance between 2 sequences of the same length
       int dist = 0;
       for (int i=0; i<seq1.length(); i++)
           if (nt_to_i(seq1.charAt(i))!= nt_to_i(seq2.charAt(i)) )
               dist++;
       
       return dist; // distance zero means 2 sequences are the same
   }
   
   
   public static List<SNP> getSNPSubList(List<SNP> list, int pos_start, int pos_end){
       // process the list and return part of that (sublist) located in the range from start_pos until pos_end
       //This is very inefficient and takes a long time, use next method
       List<SNP> subList = new ArrayList();
       int sz = list.size();
       if (sz<2)
           return null;
       int pos1 = list.get(0).position;
       int posn = list.get(sz-1).position;
       if (pos_end<pos1)
           return null;
       if (pos_start>posn)
           return null;    
       int start_search = 0;
       //for(int i=start_search)
       for(SNP snp : list){
           if (snp.position>= pos_start && snp.position<=pos_end)                
               subList.add(snp.clone()); // Note: we create a new SNP object for every element of the subList, and we didnt not use the element of main list, why? because our dag changes alleles of SNPs 
       }
       return subList;           
   }
   
      public static List<SNP> getSNPSubList(int start_search_idx, List<SNP> list, int pos_start, int pos_end){
       // process the list and return part of that (sublist) located in the range from start_pos until pos_end
       List<SNP> subList = new ArrayList();
       int sz = list.size();
       if (sz<2)
           return null;
       int pos1 = list.get(0).position;
       int posn = list.get(sz-1).position;
       if (pos_end<pos1)
           return null;
       if (pos_start>posn)
           return null;    
       for(int i=start_search_idx; i<sz; i++){
           SNP snp = list.get(i);
           if (snp.position>= pos_start && snp.position<=pos_end)                
               subList.add(snp.clone()); // Note: we create a new SNP object for every element of the subList, and we didnt not use the element of main list, why? because our dag changes alleles of SNPs 
           else if (snp.position>pos_end)
               break; // we passed by end of range
       }
       return subList;           
   }

      public static Set<Integer> getSNPPositions(List<SNP> list){
       // process the list and return part of that (sublist) located in the range from start_pos until pos_end
       Set<Integer> positions = new HashSet();
       for(SNP snp:list){
           for (int i=snp.position; i<snp.position+snp.ref.length();i++)
               positions.add(i);
       }
       return positions;           
   }

   
   public static void detectBroken(List<SNP> snpList,int max_insert_size){
       int last_snp_pos = 0; // this variable is used to determine breaking point in SNP list
       for(SNP snp: snpList){
            boolean broken = false;
            if (max_insert_size>0 && last_snp_pos>0 && (snp.position - last_snp_pos) > max_insert_size)
                snp.broken = true;
           
           last_snp_pos = snp.position;
       }
   }
  
   
  public static String getSNPAlleleFromRead(String read, Cigar ncg, int read_start, int snp_start, int snp_end){

       StringBuffer buff = new StringBuffer();
       int ref_pos = read_start; 
       int read_idx = 0; 
       int cigar_idx = 0; 
      
       int cigar_size = ncg.numCigarElements();
       try{
       boolean exit_loop = false;
       while (cigar_idx<cigar_size){
            CigarElement ce = ncg.getCigarElement(cigar_idx);
            for (int j=0; j<ce.getLength(); j++){
                if (ref_pos >= snp_start){
                   // here depending on the type of SNP, we set the value of SNP and return
                   if (!ce.getOperator().equals(CigarOperator.D))
                       buff.append(read.charAt(read_idx));
                   // see if we are at the end of current cigar element and also at the end of reference postion
                    if (ref_pos== snp_end && j==ce.getLength()-1)
                       if (cigar_idx < cigar_size-1){
                           CigarElement ce_next = ncg.getCigarElement(cigar_idx + 1);
                           if (ce_next.getOperator().equals(CigarOperator.I)){
                               int insertion_length = ce_next.getLength();
                               //we read next insertion_length bases from the read
                               //added logic on 5/3/2020 
                               if(ce.getOperator().equals(CigarOperator.D))
                                  buff.append(read.substring(read_idx , read_idx  + insertion_length));// then it is like 1D1I , does it work for 2D3I  ? It should
                               else 
                                   buff.append(read.substring(read_idx + 1, read_idx + 1 + insertion_length));// like 50M2I  
                           }
                    }
                   
                } 
                
              

            }// j for each base in a Cigar Element
            if(exit_loop)
                break;
            cigar_idx++;
       } //while Cigar element
       
   
       }catch(Exception ex){
                   System.out.println("exception inside getSNPAlleleFromRead" + ex.getMessage() );
       }
       
       return buff.toString();
   }
  
   
  
   
  
   
   public static boolean isVariation(String[] genotypes){
       //It checks if there is variation in genotypes of samples
       Set<String> gset = new HashSet();
       for (String g :genotypes)
           if(g!=null)
              gset.add(g); // in this way we get rid of duplicate bases
       
       if (gset.size()>1)
           return true;
       
        return false;  
   }
   
   public static boolean isVariationFromReference(String[] genotypes, String ref){
       //It checks if there is any variation from the reference then return true
       //current logic, if a genotype is empty string it is treated as difference from reference (or deletion); however it may be incporrect;we may change this logic in future.
       for (String g :genotypes)
           if(g!=null && !g.isEmpty() &&  !g.toLowerCase().equals(ref.toLowerCase()))
              return true; 
       
       
        return false;  
   }
   
  
   
  
   
     public static int readIdxToArrayIdx(int read_idx, int read_num){
        // we need to know if we get read idx (1-based) where in the array it should be stored
        if (read_idx>0)
            return read_idx-1;
        //so we get here when read_idx is negative
        return (read_num + Math.abs(read_idx)-1);
    }
     public static int arrayIdxToReadIdx(int arr_idx, int read_num){
        //we have array_idx, we want to know which read idx (1-based) is in the fastq file
        if (arr_idx<read_num)
            return arr_idx+1;
        //we get here if arr_idx>=read_num
        return (-(arr_idx-read_num+1));
    }
    public static int readIdxToProcessedArrayIdx(int read_idx){
        // we need to know if we get read idx (1-based) whether it has been processed or not
        
        return (Math.abs(read_idx)-1);
    }
    
    
}