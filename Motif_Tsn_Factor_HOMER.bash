#-------------------------------------------------------Motif Discovery and Transcription Factor Binding Site Analysis (HOMER)---------------------------------------------------

#Chapters:
  #Intro
  
  #What are motifs?
   #A short, recurring pattern of nucleotides that has biological meaning
   #Motifs can be found in DNA, RNA, and proteins.
   #They are often involved in the regulation of gene expression, protein folding, and other biological processes.
   # Eg:- TATA Box motif

  #Search motif  
   #Search sequences which are similar to known motif
   #Statstical method to detect the sequences which are most likely to be motif than just random sequences
  
  #Function
   #Protein folding (HTH) 
   #Gene Expression Regulation 
   #Protein Protein Interaction  

  #What are Transcription Factors?
   #Transcription factors are proteins that control the rate of transcription of genetic information from DNA to messenger RNA (mRNA).
   #Can be either DNA-binding proteins or RNA-binding proteins.
   #Bind to specific sequences of DNA called promoter sequences.
   #The binding to promoter sequences can either stimulate or repress transcription.
   #Essential for regulating gene expression, which is the process by which genes are turned on and off.
   #Can be activated by a variety of signals, including hormones, growth factors, and environmental cues.
   #Can be regulated by post-translational modifications, such as phosphorylation and acetylation
   #Based on their structure or function, some of the most common types of transcription factors include:
    #Activators:Stimulates transcription by binding to promoter sequences and recruiting RNA polymerase
    #Repressors:Represses transcription by binding to promoter sequences and preventing RNA polymerase from binding.
    #Coactivators:Work with activators to stimulate transcription.
    #Corepressors:Work with repressors to repress transcription

  #Structure of a gene and regulatory regions
   #Cis regulatory region-----------------Transcription Unit 
       #Module(enhancer)             Core Promoter--exon--intron--UTR
                                      #Promoters
                                      #Enhancers
                                      #Silencers
                                      #Insulators

  #Why is there a need for gene regulation?
   #Positive regulation:Activators bind to DNA sequences near the promoter of a gene, which is the region of DNA that RNA polymerase binds to initiate transcription.
                      #:Once bound, activators recruit RNA polymerase to the promoter, which allows transcription to begin.
   #Negative regulation:Repressors bind to DNA sequences near the promoter of a gene, which prevents RNA polymerase from binding.
                       #:This prevents transcription from occurring.
   #Long-range regulatory elements (LCR--Enhancer--Silencer-Insulator)
   #Promoter Elements(Proximal Promoter Elements--Core Promoter ELements)
   #Gene


  #Different ways to represent Motifs or Transcription Factor Binding Sites(TFBS)
   #Consensus Sequences: Most frequent occuring nucleotide at given position (A,G,C,T)(Y=C/T,H=A/C/T)
   #Position Frequency Matrix(PFM):Graph of nucleotide position vs counts
    #PFM is used to create: Position Weight Matrix(PWM)/Position Specific Weight Matrix(PSWM)/Position Specific scoring Matrix(PSSM)
     #PWM uses log odd score (Represent the Probability of occuring of each nucleotide at a speceific position)
      #Higher log odd score = High probability of occurence of that base at that position
    #Sequence Logos:Graphical representation of the PWM (Nucleotide position vs Bits(0 to 2))
 
  #When to perform a Transcription Factor Binding Site Analysis?
   ##Analysis
    #Motif discovery
     ##Research questions
      #What are the binding preferences of transcription factors?
      #How do transcription factors regulate gene expression?
      #How does DNA sequence variation affect transcription factor binding?
   ##Analysis
    #Motif enrichment
     ##Research questions
      #What genes are regulated by a particular transcription factor? 
      #How do motifs contribute to the regulation of gene expression in different cell types or under different conditions?
      #How do motifs contribute to the development of diseases?
      #How can transcription factors be targeted for drug development?

  #Tools available to perform motif analysis
   #HOMER(Hypergeometric Optimization of Motif EnRichment)
   #The MEME Suite
   #JASPAR  
   #geneXplain
   #Genome Enhancer(Identify master Regulator )


  #Motif databases
   #TRANSFAC 2.O
   #JASPAR
  
  #Case study for demonstration using Homer
   #Background:RUNX1 is a transcription factor that plays a critical role in the development of breast cancer.
              #:It is involved in the regulation of genes that control cell growth, differentiation, and apoptosis.
              #:Overexpression of RUNX1 can lead to the activation of genes that promote cancer growth and survival.
   #Study design:This study used ChIP-Seq to identify genomic binding of RUNX1. Study generated two samples from MCF10A cells: one input control and one RUNX1 ChIP-Seq.
   #Goal of the Analysis:To identify DNA motifs (de novo and canonical) that are enriched in the genomic regions where RUNX1 binds 

  #Fetch data from NCBI GEO
   #https://www.ncbi.nlm.nih.gov/geo/quer...https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129314 
   #gunzip GSE129314_RUNX1ChIP_MCF10AOEq0.05_peaks.txt.gz
   #cat GSE129314_RUNX1ChIP_MCF10AOEq0.05_peaks.txt | less  

  #About HOMER and notes on installation
   #http://homer.ucsd.edu/homer/
   #Configure conda channels (do this once)
     conda config --add channels defaults
     conda config --add channels bioconda
     conda config --add channels conda-forge
   #Create a dedicated environment (recommended)
     conda create -n homer_env homer -y
   #Activate it:
     conda activate homer_env
   #Verify installation
     findMotifsGenome.pl -h
   #Download hg38 FASTA (from UCSC mirror)
    perl /home/aman/miniconda3/envs/homer_env/share/homer/.//configureHomer.pl -install hg38
    curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip hg38.fa.gz
   #Install Genome Data (Required for Motif/Peak Analysis)
     perl $(which configureHomer.pl) -install hg38
   #Verify installation
     perl $(which configureHomer.pl) -list

  #Following Homer instructions to find motifs in genomic regions
   #http://homer.ucsd.edu/homer/ngs/peakMotifs.html
    findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
    i.e. findMotifsGenome.pl ERpeaks.txt hg18 ER_MotifOutput/ -size 200 -mask
   #Acceptable Input files
    #findMotifsGenome.pl accepts HOMER peak files or BED files:
     #HOMER peak files should have at minimum 5 columns (separated by TABs, additional columns will be ignored):
     #Column1: Unique Peak ID
     #Column2: chromosome
     #Column3: starting position
     #Column4: ending position
     #Column5: Strand (+/- or 0/1, where 0="+", 1="-")
    #BED files should have at minimum 6 columns (separated by TABs, additional columns will be ignored)
     #Column1: chromosome
     #Column2: starting position
     #Column3: ending position
     #Column4: Unique Peak ID
     #Column5: not used
     #Column6: Strand (+/- or 0/1, where 0="+", 1="-")
    #In theory, HOMER will accept BED files with only 4 columns (+/- in the 4th column), and files without unique IDs, 
    #but this is NOT recommended.  For one, if you don't have unique IDs for your regions, it's hard to go back and figure out which region contains which peak.


  


  #SCRIPT

#!/bin/bash


pwd

  # adjust paths to files based on your current working directory

  #Manipulate peak file to get it into the acceptable format
  # Run findMotifsGenome.pl on genomic regions
  # To find line number for header line in the peak file
gunzip GSE129314_RUNX1ChIP_MCF10AOEq0.05_peaks.txt.gz
cat GSE129314_RUNX1ChIP_MCF10AOEq0.05_peaks.txt | less
cat GSE129314_RUNX1ChIP_MCF10AOEq0.05_peaks.txt | grep -n ^chr | less
cat GSE129314_RUNX1ChIP_MCF10AOEq0.05_peaks.txt | tail -n +30 | less
cat GSE129314_RUNX1ChIP_MCF10AOEq0.05_peaks.txt | tail -n +30 | cut -f 1-3,10 |less
cat GSE129314_RUNX1ChIP_MCF10AOEq0.05_peaks.txt | tail -n +30 | cut -f 1-3,10 > GSE129314_RUNX1ChIP_peaks_homer.txt

  #Run Homer script to find motifs
  # Run Homer
findMotifsGenome.pl GSE129314_RUNX1ChIP_peaks_homer.txt hg38 Result/ -size 200
  # size = size of the region for motif finding
  # For recommendations to set this parameter, refer: http://homer.ucsd.edu/homer/ngs/peakMotifs.html

  #Method 1 to find motif instances
  # find motif instances - method 1
findMotifsGenome.pl GSE129314_RUNX1ChIP_peaks_homer.txt hg38 Results/ -find  Results/knownResults/known41.motif > Results/motif_instances.txt

  #Method 2 to find motif instances
  # find motif instances - method 2
annotatePeaks.pl GSE129314_RUNX1ChIP_peaks_homer.txt hg38 -m Results/knownResults/known41.motif > Results/motif_instances_annotatePeak_method.txt


  #Looking at resulting files
  #Understanding Homer’s de novo motif finding results
   #homerMotifs.all.motifs:All the motifs compiled in one file
   #homerResults:Hold motif files that were used to find instances of each motif
   #knownResults.txt: Stats associated with canonical motif that we found          
   #seq.autonorm.tsv: Autonormalization stats, steps that were performed by homer to account for the biased institute by shorter oligonucleotide.
   #homerMotifs.motifs10; homerMotifs.motifs12; homerMotifs.motifs8; Hold the result of the denovo motif searched for (10,12,8-K mer size)    
   #homerResults.html: Result of the canonical motif (IMPORTANT)
   #motifFindingParameters.txt:command that were used to execute the script 
   #knownResults:Hold motif files that were used in finding specific instances of each motif       
   #motif_instances.txt:A list of genomic locations where a specific motif is found seq peaks
   #knownResults.html:Result of the canonical result findings 
   #nonRedundant.motifs:A curated set of unique motifs after removing:Duplicate motifs; Highly similar motifs; Reverse complements

  #Understanding Homer’s known motif (canonical) finding results
  #Find which peaks do a specific motif bind to?
   #To find the location where p53 actually bind in the genomic region go on HOMER page (http://homer.ucsd.edu/homer/ngs/peakMotifs.html)
    #Go to: Finding Instance of Specific Motifs
     #By default, HOMER does not return the locations of each motif found in the motif discovery process.  To recover the motif locations, you must first select the motifs you're interested in by getting the "motif file" output by HOMER.  You can combine multiple motifs in single file if you like to form a "motif library".  To identify motif locations, you have two options:
    #1. Run findMotifsGenome.pl with the "-find <motif file>" option.  This will output a tab-delimited text file with each line containing an instance of the motif in the target peaks.  The output is sent to stdout.
     #For example: findMotifsGenome.pl ERalpha.peaks hg18 MotifOutputDirectory/ -find motif1.motif > outputfile.txt
findMotifsGenome.pl GSE129314_RUNX1ChIP_peaks_homer.txt hg38 ./ -find Results/knownResults/known44.motif > motif_instances.txt
     #motif_instances: Contain results
     #The output file will contain the columns:
     #Peak/Region ID
     #Offset from the center of the region
     #Sequence of the site
     #Name of the Motif
     #Strand
     #Motif Score (log odds score of the motif matrix, higher scores are better matches)
    #2. Run annotatePeaks.pl with the "-m <motif file>" option (see the annotation section for more info).  Chuck prefers doing it this way.  This will output a tab-delimited text file with each line containing a peak/region and a column containing instance of each motif separated by commas to stdout
     #For example: annotatePeaks.pl ERalpha.peaks hg18 -m motif1.motif > outputfile.txt
annotatePeaks.pl GSE129314_RUNX1ChIP_peaks_homer.txt hg38 -m Results/knownResults/known44.motif > motif_instances2.txt
     #The output file will contain columns:
     #Peak/Region ID
     #Chromosome
     #Start
     #End
     #Strand of Peaks
     #6-18: annotation information
     #19. CpG%
     #20. GC%
     #21. Motif Instances
     #Motif Instances have the following format:
     #<distance from center of region>(<sequence>,<strand>,<conservation>)
     #i.e -29(TAAATCAACA,+,0.00)
     #You can also find histogram of motif density this way by adding "-hist <#>" to the command.  For example:
annotatePeaks.pl ERalpha.peaks hg18 -m ere.motif foxa1.motif -size 1000 -hist 10 > outputfile.txt
 



#DATA SOURCES
#▸ Homer Vignette:
  #http://homer.ucsd.edu/homer/ngs/peakMotifs.html

#▸ Data used for demonstration:
  #https://www.ncbi.nlm.nih.gov/geo/quer...https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129314

#▸ Link to Code:
  #https://github.com/kpatel427/YouTubeTutorials/blob/main/runHomer.sh





#-------------------------------------------------------Motif Discovery and Transcription Factor Binding Site Analysis (HOMER)---------------------------------------------------
