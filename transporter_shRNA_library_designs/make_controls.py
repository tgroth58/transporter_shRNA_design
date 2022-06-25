import os
import random
import argparse

def read_seqs(seqLoc):
    '''
    Reads in a list of transcript files given a file location:
    input - seqLoc: location of transcript files to parse
    output : dictionary of transcript sequences
    '''
    seq_dict={}
    for f in os.listdir(seqLoc):
        txName=f.split('_')[0]
        with open(os.path.join(*[seqLoc,f]),'r') as fl:
            content=fl.readlines()
        seq_dict[txName]=''.join([sq.strip() for sq in content[1:]])
    return(seq_dict)

def permute_seq(seq):
    '''
    Permutes a given sequence, scrambling the order of the bases
    input - seq: A nucleotide sequence to be scrambled
    output : the given sequence scrambled
    '''
    scrambled=''.join(random.sample(list(seq),len(seq)))
    return(scrambled)

def create_design_template(seq,temp):
    '''
    Creates temporary scramble sequence file
    input - seq: DNA sequence to be written
    output - a tempate sequence
    '''
    with open(temp,'w') as f:
        f.write('>TmpSeq\n')
        f.write(seq)


def run_shRNA_design(inFile,outFile):
    '''
    Runs the shRNA designer at the OS level
    Location of script is hard-coded into this function
    input - seq: Transcript sequence to encode
    output : List of shRNAs 
    '''
    #Function Variables:
    shRNA_Tool='/home/ted/Documents/Programs/si_shRNA_design/si_shRNA_selector-x64'
    command=' '.join(['-dGmin -32','-dGmax -28','-i',inFile,'-o',outFile])
    os.system(command)
    print("Created designs for "+inFile)

def parse_designs(inFile):
    '''
    Parses the design file for shRNAs that are efficient
    input - seq: Design File created by si_shRNA_selector-x64
    output : List of shRNAs 
    '''
    with open(inFile,'r') as f:
        designs=f.readlines()
    parse_designs=[x for x in designs if x[-1]=='efficient']
    return(parse_designs)

def create_blast_query(designs,tmp_query):
    '''
    Creates temporary query file for blast to use to search for transcripts
    input - designs: list of lists of shRNA designs
    input - tmp_query: name of temporary file containing shRNA sequence query info
    '''
    u_trdict={'U':'T'}
    with open(tmp_query,'w') as f:
        for i,d in enumerate(designs):
            f.write('>'+i+'\n')
            f.write(d[1].upper().translate(u_trdict)+'\n')

def blast_designs(designs,tmp_query,db,tmp_out):
    '''
    Runs a BLAST at the system level to check for off target effects on any shRNAs
    input - designs: list of lists of shRNA designs
    output - blast results
    '''
    #Function Variables:
    blastn='blastn'
    command=' '.join([blastn,'-query',tmp_query,'-db',db,'-task blastn-short','-outfmt 6','-out',tmp_out])
    print('BLASTing designs against '+db+' database')
    os.system(command)

def get_offTargets(tmp_blastOut):
    '''
    Returns a list of sequence entries that have 
    high degree of matching to any shRNA designs
    '''
    with open(tmp_blastOut,'r') as f:
        blastRes=f.readlines()

    # Start and End of the query match must be 1 and 17, respectively
    # Number of mismatches must be less than or equal to 2
    offTargets=[i for i,x in enumerate(blastRes) if x[]==1 and x[]==19 and x[]<=2]
    return(offTargets)

def design_pipeline(seqs,
        tmp_design_in='tmp_design.fa',tmp_design_out='tmp_design_out.txt',
        tmp_blast_in='tmp_blast_query.fa',tmp_blast_out='tmp_blast_res.txt',
        db)
    ''' 
    Runs the entire design pipeline until "target" number of guides designed
    '''
    #1. Select a sequence randomly and permute it:
    seq=permute_seq(random.sample(seqs,1))
    create_design_template(seq,tmp_design_in)
    #2. Make shRNAs and parse efficient ones:
    run_shRNA_design(tmp_design_in,tmp_design_out)
    designs=parse_designs(tmp_design_out)
    #3. Run BLAST on designs:
    create_blast_query(designs,tmp_blast_in)
    blast_designs(designs,tmp_blast_in,db,tmp_blast_out)
    #4. Parse BLAST results with low matching rates to other genes
    offTargets=get_offTargets(tmp_blast_out)
    specific_shRNA=[designs[i] for i in range(len(designs)) if i not in offTargets]
    # Remove temporary files:
    os.remove([tmp_design_in,tmp_design_out,tmp_blast_in,tmp_blast_out])    
    #Return the shRNA designs:
    return(specific_shRNA)

def parse_arguments():
    parser=argparse.ArgumentParser(description='Random Control shRNA designer')

    parser.add_argument('--seqs',
            help="Location of sequence files to create transcripts.")
    parser.add_argument('--db',
            help="BLAST Database used for off-target search.")
    parser.add_argument('--target',type=int,
            help="Number of total control shRNAs desired.")
    parser.add_argument('--output',
            help="Name of the output file for created shRNAs.")

def design_main(seqs,db,target,output_designs):
    '''
    Runs the design pipeline until the "target" number of designs have been created
    '''
    #Parse arguments:
    seqs=args.seqs
    db=args.db
    target=args.target
    output_designs=args.output
    #Temporary file names
    tmp_design_in='tmp_design.fa',
    tmp_design_out='tmp_design_out.txt',
    tmp_blast_in='tmp_blast_query.fa',
    tmp_blast_out='tmp_blast_res.txt',

    final_designs=[]
    rounds=1
    while len(final_designs)<target:
        uniqueSeqs=[y[1] for y in final_designs]
        print('Round '+rounds+' of design')
        new_designs=design_pipeline(seqs,
            tmp_design_in,
            tmp_design_out,
            tmp_blast_in,
            tmp_blast_out,
            db
        )
        new_designs=[x for x in new designs if x[1] not in uniqueSeqs]
        print('{0} specific shRNAs found'.format(str(len(new_designs))))
        new_designs.append(final_designs)
        rounds+=1

    #Parse the requested number of targets
    final_designs=[0:target]
    return(final_designs)

if name=='__main__':

