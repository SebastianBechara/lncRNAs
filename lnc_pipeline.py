##########lnc pipeline script
import os, sys
from Bio import SeqIO
from Bio.Seq import Seq

F_exit_cnt = 0
S_exit_cnt = 0

def check_folders():#deepseq):
    if os.path.isdir('temp') != True:
        os.mkdir('temp')
    if os.path.isdir('result') != True:
        os.mkdir('result')
#        os.mkdir('result/')
#    os.mkdir('result/lnc_candidates'+str(deepseq).split('/')[-1])

def check_file(file): ### only if theres one pipeline; with SPAdes, will have to use 2 cnters
    if len( open(str(file)).readlines()) <=2:
        sys.exit()

def F_exit(file, F_exit_cnt):
    num_lines = len(open(str(file)).readlines())
    if num_lines > 0:
        print('\n####\n\nThe FEELnc output file contains '+ str(num_lines) +' lines!\n\n####\n')
    else:
        F_exit_cnt += 1
        print('\n####\n\nThe FEELnc filter found no candidates\n\n####\n')
    return F_exit_cnt

def S_exit(file, S_exit_cnt):
    num_lines = len(open(str(file)).readlines())
    if num_lines > 0:
        print('\n####\n\nThe SPAdes assembler pipeline output file contains '+ str(num_lines) +' lines!\n\n####\n')
    else:
        S_exit_cnt += 1
        print('\n####\n\nThe SPAdes assembler pipeline found no candidates\n\n####\n')
    return S_exit_cnt

def sys_exit(F_exit_cnt, S_exit_cnt):
    if F_exit_cnt == 1 and S_exit_cnt == 1:
        print('\n\n####\n\nBoth FEELnc and SPAdes pipeline found no candidates!\n\n####\n\nExiting now . . .\n\n####')
        sys.exit()

def trim_seq_data(file, minlen):
    print('\n####\nTrimming sequencing data . . .\n####\n')
    output= 'temp/'+str(file).split('/')[-1]+'_trimmed_15_'+str(minlen)+'.fq'
    os.system('bbduk.sh -Xmx2g in=' + str(file) + ' out=' + output + ' ref=/home/debian9/host/bbmap/resources/adapters.fa ktrim=l mink=11 qtrim=rl trimq=15 minlen=' + str(minlen))
    return output

def map_Hisat(file, ref_path):
    print('\n####\nMapping sequencing data . . .\n####\n')
    output= str(file) + '_mapped.sam'
    os.system('hisat2 --dta -x ' + str(ref_path) + ' -U ' + str(file) + ' -S ' + str(file) + '_mapped.sam')
    return output

def filter_sort(file):
    print('\n####\nFiltering and sorting mapped reads . . .\n')
    output = str(file) + '_trimmed_mapped_sorted'
    os.system('samtools view -Sb -F 4 '+ str(file) + ' | samtools sort -o '  + output) # > ' + str(file)+'trimmed_mapped_sorted.bam')
    return output

def bamToGtf(file, guide_ref):
    print('Done!\n\n####\nConverting mapped reads from BAM to GTF . . .\n####\n')
    output=str(file)+'.gtf'
    os.system('stringtie ' + str(file)+' -G ' + str(guide_ref) + ' -o ' + output)
    return output

def gff_to_gtf(ref_anno):
    print('Done!\n\n####\nChecking refference annotation file for FEELnc . . .\n####\n')
    output= str(ref_anno) + '.gtf'
    if ref_anno.endswith('.gtf') == False:
        os.system('gffread ' + str(ref_anno) + ' -T -o '+ output)
        return output
    else:
        return ref_anno

def FEELnc(file, ref):
    print('\n####\nRunning the FELLnc filter module . . .\n####\n')
    output= str(file) +'_candidates.gtf'
    os.system('FEELnc_filter.pl -i ' + str(file) + ' -a ' + str(ref) + ' > '+ output)
    return output

def gtf_to_fasta(file, genome):
    print('\n####\nConverting FEELnc output (GTF) to fasta . . .\n####\n')
    output=str(file) + '.fasta'
    os.system('gffread '+str(file)+' -g '+str(genome) +' -w '+str(file)+'.fasta')
    fasta = {str(i.description):str(i.seq) for i in SeqIO.parse(str(output), 'fasta')}
    with open(str(output),'w') as outfile:
        for k,v in fasta.items():
            outfile.write('>FEELnc_'+k+':\n'+v+'\n')
    return output
###################### not incorporated in main()
def SPAdes(S_input):
    print('\n####\nBuilding de novo transcriptom assembly . . .\n####\n')
    output = ('temp/SPAdes_'+str(S_input).split('/')[-1])
    os.system('\nspades.py --rna -s '+ str(S_input) + ' -o ' + output) #### if on kubuntu just spades but spades.yp also works
    return output
###################### not incorporated in main()
def sort_S_output(S_output):
    print('\n####\nSorting SPAdes output . . .\n####\n')
    trans_file = str(S_output)+'/transcripts.fasta'
    transcripts = {str(i.description):str(i.seq) for i in SeqIO.parse(trans_file, 'fasta')}
    plus200bp ={}
    output = 'temp/SPAdes_+200transcripts.fasta'
    for k,v in transcripts.items():
        if len(v) > 200:
            plus200bp[k] = v
    with open(str(output), 'w') as outfile:
        for k,v in plus200bp.items():
            outfile.write('>'+k+':\n'+v+'\n')

    return output, plus200bp ##########(supposed to be a fasta file not a dict)!!! is now both, should help later
###################### not incorporated in main()
def make_blast_database(fasta):
    print('Done!\n####\nBuilding BLAST database . . .\n####\n')
    output = 'temp/' + str(fasta).split('/')[-1].split('.')[0] + '_blastDB'
    print(fasta, str(fasta).split('/')[-1].split('.')[0])
    print(output)
    os.system('makeblastdb -in '+str(fasta)+' -dbtype nucl -parse_seqids -out ' + str(output))
    return output

def blastn(input, blast_database):
    print('\n####\nBlasting transcriptom against build database . . .\n####\n')
    print(blast_database+'\n')
    output = str(input)+'_BLAST_trans_vs_'+str(blast_database).split('/')[-1]
    os.system('blastn -query '+str(input)+' -db '+str(blast_database)+' -out '+str(output)+" -outfmt '6'") ################
    return output

def filter_blastn(input, transdict):
    print('Done!\n####\nFiltering BLAST output . . .\n####\n')
    good_cov_lines=[]
    good=[]
    bad=[]
    plus200fa = transdict
    print('plus200fa:', len(list(plus200fa.keys())))
    #plus200fa = {str(i.description):str(i.seq) for i in SeqIO.parse(str(plus200bp), 'fasta')}
    candidates = {}
    output = 'temp/'+str(input).split('/')[-1]+'SPAdes_lnc_candidates.fasta'
    print('\n\t## Building CDS dictionary . . . \t')
    CDSs = []
    for line in open('/home/debian9/host/annotation_files/ptet_51_mac_CDSsonly_annotation.gff3').readlines()[1:]:
        scaff = line.split('\t')[0]
        start = int(line.split('\t')[3])
        stop = int(line.split('\t')[4])
        CDSs.append((scaff, start, stop))

    print('\n\t## checking BLAST coverage . . . \t')
    lines = open(str(input)).readlines()
    for line in lines:
        print('\t##\tProcessing line: '+ str(lines.index(line))+' of '+ str(len(lines)), end='\t\t\r')
        if (100 / int(line.split('\t')[0].split('_')[3])) * int(line.split('\t')[3]) >= 50:
            good_cov_lines.append(line)

    print('\n\n\t## Filtering out candidates . . . \t')
    ############################################################################
    #    start and stop is not ordered, might mess up the actuall sorting      #
    #    idea, put start stop into lst, sort it and the call the list when     #
    #    sorting the input/output                                              #
    ############################################################################
    for line in good_cov_lines:
        scaf = line.split('\t')[1] ## is cds[0]
        start = int(line.split('\t')[8]) ## is cds[1]
        stop = int(line.split('\t')[9]) ## is cds[2]
        pos=[start, stop]
        pos.sort()
        print('\t##\tProcessing line: ', good_cov_lines.index(line),' of ', len(good_cov_lines),' @ ',scaf, end='\t\t\r')
        for cds in CDSs:
            ### new code that should make it better, now it catches transcripts it missed earlier
            if scaf == cds[0]:
                if pos[0] <= int(cds[2]):  ##  start <= nc[2]-end
                    ### Within CDS
                    if int(cds[1]) <= pos[1]:  ##  end >= nc[1]-start
                        bad.append(line)
                        break
                    if int(cds[1]) <= pos[0] <= int(cds[2]):  ##  nc[1]-start <= start <= nc[2]-end
                        bad.append(line)
                        break
                else:
                    continue
        ##### my own slightly messed up code
        # for cds in CDSs:
        #     if scaf == cds[0]:
        #         if pos[0] <= int(cds[1]):
        #             if int(cds[1]) <= pos[1] <= int(cds[2]):
        #                 bad.append(line)
        #                 break
        #             elif int(cds[2]) <= pos[1]:
        #                 bad.append(line)
        #                 break
        #         elif int(cds[1]) <= pos[0]:
        #             if pos[1] <= int(cds[2]):
        #                 bad.append(line)
        #                 break
        #             elif int(cds[2]) <= pos[1]:
        #                 bad.append(line)
        #                 break
    print('\t##\tProcessing line: ', good_cov_lines.index(line),' of ', len(good_cov_lines),' @ ',scaf)
    for line in good_cov_lines:
        if line not in bad:
            good.append(line)
    for line in good:
        candidates.setdefault(str(line.split('\t')[0]), []).append(line)
    #do_they_fit = []
    with open(output, 'w') as outfile:
        print('\n\n\t## Writing to file . . . \t')
        #outfile.write('## this file contains the lncRNA candidates obtained by the SPAdes pipeline\n')
        for key in list(candidates.keys()):
            for k,v in plus200fa.items():
                #do_they_fit.append((key, k))
                if key.strip(':') == k:
    #                do_they_fit.append((key, k))
                    outfile.write('>SPAdes_'+str(key)+'\n'+str(v)+'\n')
    #print(len(do_they_fit), do_they_fit[0])
    return output

def cat_fastas(FEELnc_fasta, SPAdes_fasta):
    print('\n####\nConcatenating fasta file obtained from FEELnc and SPAdes pipeline . . .\n###\n')
    output = 'temp/combined_lnc_candidates.fasta'
    os.system('cat ' + str(FEELnc_fasta) +' '+ str(SPAdes_fasta) + ' > '+output)
    return output

def all6RF(file):
    print('Done!\n\n####\nTranslating all reading frames . . .\n####\n')
    candidates = { str(i.description):str(i.seq) for i in SeqIO.parse(str(file), 'fasta')}
    cand6RF = []
    output='temp/lnc_candidates_all6RF.fasta'
    for k,v in candidates.items():
        if k.startswith('FEEL'):
            RF1 = str(Seq(v).translate(table=6)).replace('*','X')
            RF2 = str(Seq(v)[1:].translate(table=6)).replace('*', 'X')
            RF3 = str(Seq(v)[2:].translate(table=6)).replace('*', 'X')
            RF_1 = str(Seq(v).reverse_complement().translate(table=6)).replace('*','X')
            RF_2 = str(Seq(v).reverse_complement()[1:].translate(table=6)).replace('*','X')
            RF_3 = str(Seq(v).reverse_complement()[2:].translate(table=6)).replace('*','X')
            cand6RF.append('>'+str(k)+'_RF1:\n'+RF1+'\n')
            cand6RF.append('>'+str(k)+'_RF2:\n'+RF2+'\n')
            cand6RF.append('>'+str(k)+'_RF3:\n'+RF3+'\n')
            cand6RF.append('>'+str(k)+'_RF_1:\n'+RF_1+'\n')
            cand6RF.append('>'+str(k)+'_RF_2:\n'+RF_2+'\n')
            cand6RF.append('>'+str(k)+'_RF_3:\n'+RF_3+'\n')
        elif k.startswith('SPA'):
            RF1 = str(Seq(v).translate(table=6)).replace('*','X')
            RF1 = str(Seq(v).translate(table=6)).replace('*','X')
            RF2 = str(Seq(v)[1:].translate(table=6)).replace('*', 'X')
            RF3 = str(Seq(v)[2:].translate(table=6)).replace('*', 'X')
            RF_1 = str(Seq(v).reverse_complement().translate(table=6)).replace('*','X')
            RF_2 = str(Seq(v).reverse_complement()[1:].translate(table=6)).replace('*','X')
            RF_3 = str(Seq(v).reverse_complement()[2:].translate(table=6)).replace('*','X')
            cand6RF.append('>'+str(k)+'_RF1:\n'+RF1+'\n')
            cand6RF.append('>'+str(k)+'_RF2:\n'+RF2+'\n')
            cand6RF.append('>'+str(k)+'_RF3:\n'+RF3+'\n')
            cand6RF.append('>'+str(k)+'_RF_1:\n'+RF_1+'\n')
            cand6RF.append('>'+str(k)+'_RF_2:\n'+RF_2+'\n')
            cand6RF.append('>'+str(k)+'_RF_3:\n'+RF_3+'\n')
#    for k,v in candidates.items():
#        RF1 = str(Seq(v).translate(table=6)).replace('*','X')
#        RF2 = str(Seq(v)[1:].translate(table=6)).replace('*','X')
#        RF3 = str(Seq(v)[2:].translate(table=6)).replace('*','X')
#        RF_1 = str(Seq(v).reverse_complement().translate(table=6)).replace('*','X')
#        RF_2 = str(Seq(v).reverse_complement()[1:].translate(table=6)).replace('*','X')
#        RF_3 = str(Seq(v).reverse_complement()[2:].translate(table=6)).replace('*','X')
#        cand_6RF.append('>FEELnc_' + str(k) + '_RF1:\n' + RF1 + '\n')
#        cand_6RF.append('>FEELnc_' + str(k) + '_RF2:\n' + RF2 + '\n')
#        cand_6RF.append('>FEELnc_' + str(k) + '_RF3:\n' + RF3 + '\n')
#        cand_6RF.append('>FEELnc_' + str(k) + '_RF_1:\n' + RF_1 + '\n')
#        cand_6RF.append('>FEELnc_' + str(k) + '_RF_2:\n' + RF_2 + '\n')
#        cand_6RF.append('>FEELnc_' + str(k) + '_RF_3:\n' + RF_3 + '\n')

    with open(output, 'w') as outfile:
        for i in cand6RF:
            outfile.write(i)

    return output

def hmmsearch(database, file, comb_fasta, title):
    print('Done!\n\n####\nChecking candidates for significant protein domains . . .\n####\n')
    midput='temp/putative_lnc_candidates.txt'
    output = 'temp/FEELnc-SPAdes_putative_lncRNAs_'+str(title.split('/')[-1].split('_')[0])+ '.txt'
    output_fasta = 'temp/FEELnc-SPAdes_putative_lncRNAs_'+str(title.split('/')[-1].split('_')[0])+ '.fasta'
####output2 = 'result/FEELnc_putatie_lncRNAs.fasta'########################
    os.system('hmmsearch --tblout '+ midput +' '+ str(database)+ ' '+str(file))
#############
    print('Done!\n\n####\nCleaning up lncRNA candidates . . .\n####\n')
    lnc={ str(i.description):[] for i in SeqIO.parse(str(comb_fasta), 'fasta')}#'/home/debian9/host/analysis/lncRNA/lncRNA_candidates_mod_anno.fasta' ,'fasta')}
    lnc_lst=[]
    candidates={ str(i.description):str(i.seq) for i in SeqIO.parse(str(comb_fasta), 'fasta')}

    with open(midput, 'r') as infile:
        lines=infile.readlines()

        for line in lines:
            if line.startswith('#') == False:
                if float(line.split()[4]) < 1e-3 and float(line.split()[7]) < 1e-3:
                    lnc.setdefault(str(line.split()[0].split('RF')[0].strip('_')), []).append(line)

        for k,v in lnc.items():
            if len(v) == 0:
                lnc_lst.append(k)

    with open(output,'w') as outfile:
        outfile.write('## filtered putative lncRNAs from: ' + str(title.split('/')[-1].split('_')[0])+'\n')
        for i in lnc_lst:
            outfile.write(str(i)+'\n')
            #for k,v in lnc.items():
                #if i == k:
                    #outfile.write('>'+str(i)+'\n'+str(v)+'\n')
    with open(output_fasta, 'w') as outfile:
        for i in lnc_lst:
            for k,v in candidates.items():
                if i == k:
                    if i.startswith('FEEL'):
                        outfile.write('>'+str(k).strip(':')+'_length='+str(len(v))+'\n'+str(v)+'\n')
                    elif i.startswith('SPA'):
                        outfile.write('>'+'_'.join(str(k).split('_')[:5]+str(k).strip(':').split('_')[-2:])+'\n'+str(v)+'\n')
                    #outfile.write('>'+str(k)+'\n'+str(v)+'\n')

#                    if i.stratswith('FEEL'):
#                        outfile.write('>'+str(k).strip(':')+'_length='+str(len(v))+'\n'+str(v)+'\n')
    return output, output_fasta

def split_cand(fasta_file):
###############################################
    fasta={str(i.description):str(i.seq) for i in SeqIO.parse(str(fasta_file), 'fasta')}
    feel={}
    spades={}
    outfileFeel='result/FEELnc_lnc_cand.fasta'
    outfileSpades='result/SPAdes_lnc_cand.fasta'
################################################################################
#                   repair weird ass shit here
################################################################################
    for k,v in fasta.items():
        if k.startswith('FEEL'):
        #feel.setdefault(k, []).append(v)
            feel[k]=v
        if k.startswith('SPA'):
        #spades.setdefault(k, []).append(v)
            spades[k]=v

    with open(outfileFeel, 'w') as outfile:
        for k,v in feel.items():
            outfile.write('>'+str(k)+'\n'+str(v)+'\n')
        #outfile.write('>'+str(k).strip(':')+'_length='+str(len(v))+'\n'+str(v)+'\n')

    with open(outfileSpades, 'w') as outfile:
        for k,v in spades.items():
            outfile.write('>'+str(k)+'\n'+str(v)+'\n')
        #outfile.write('>'+'_'.join(str(k).split('_')[:5]+str(k).strip(':')\
        #.split('_')[-2:])+'\n'+str(v)+'\n') ###### you have to shorten the id otherwise you cant build a blast DB

    return outfileFeel, outfileSpades
    #############################
    # blast them against each other
    #############################
def blast_ungapped(fasta_file, database):
    print('\n\n#### Blasting '+str(fasta_file).split('/')[-1]+' against '+str(database).split('/')[-1])
    output = 'temp/'+str(fasta_file).split('/')[-1].split('_')[0]+'_vs_'+str(database).split('/')[-1]+'_ungapped_rep'
    os.system("blastn -query "+str(fasta_file)+" -db "+str(database)+" -ungapped -perc_identity 97 -out "+str(output)+" -outfmt '6'")
    return output

def check_same_fromOne(blast_rep):
    print('\n\n\tDone!\n\n#### Checking overlapping candidates\n')
    print('\t### checking for shared lncRNA candidates from:\n'+str(blast_rep)+'\n')
    keep = []
    final_same = []
    same_same= []
    lines = open(str(blast_rep)).readlines()
    for line in lines:
        print('\t\t##\tprocessing line '+str(lines.index(line))+' of '+str(len(lines)), end='\t\t\r')
        query_name=line.split('\t')[0]
        #print(query_name)
        ref_name=line.split('\t')[1]
        #print(ref_name)
        query_pos=[int(line.split('\t')[6]), int(line.split('\t')[7])]
        ref_pos=[int(line.split('\t')[8]), int(line.split('\t')[9])]
        query_pos.sort() ##  [0] == start, [1] == end
        ref_pos.sort() ##  [0] == start, [1] == end
        #print(query_pos, ref_pos)
        if line.startswith('FEEL'):
            query_length=int(line.split('\t')[0].split('=')[-1])
        #    print(query_length)
            ref_length=int(line.split('\t')[1].split('_')[4])
        #    print(ref_length)
        elif line.startswith('SPA'):
            query_length=int(line.split('\t')[0].split('_')[4])
        #    print(query_length)
            ref_length=int(line.split('\t')[1].split('=')[-1])
        #    print(ref_length)
        if ref_length == ref_pos[1]:
            if ref_length > query_length:
                keep.append(ref_name)
            elif query_length > ref_length:
                keep.append(query_name)
            elif ref_length == query_length:
                keep.append(query_name)
                if query_name.startswith('FEEL'):
                    same_same.append((query_name, ref_name))
                elif query_name.startswith('SPA'):
                    same_same.append((ref_name, query_name))
        elif query_length == query_pos[1]:
            if query_length > ref_length:
                keep.append(query_name)
            elif ref_length > query_length:
                keep.append(ref_name)
            elif ref_length == query_length:
                keep.append(query_name)
                if query_name.startswith('FEEL'):
                    same_same.append((query_name, ref_name))
                elif query_name.startswith('SPA'):
                    same_same.append((ref_name, query_name))
    for i in keep:
        if i not in final_same:
            final_same.append(i)
#    print('Same_same: ',same_same,'length of final_same:', len(final_same))
    return final_same, same_same

def check_same_final(final_FEELnc, final_SPAdes, deepseqFile, fasta_file, same_list_FEEL, same_list_SPA):
    final =[]
    fasta = {str(i.description):str(i.seq) for i in SeqIO.parse(str(fasta_file), 'fasta')}
    output = 'result/FEELnc-SPAdes_overlapping_lncRNA_candidates_from_'+str(deepseqFile).split('/')[-1].split('_')[0]+'.txt'
    output_fasta = 'result/FEELnc-SPAdes_overlapping_lncRNA_candidates_from_'+str(deepseqFile).split('/')[-1].split('_')[0]+'.fasta'
    for i in final_FEELnc:
        ###
        if same_list_FEEL == []:
            if i not in final:
                final.append(i)
        else:
        ###
            for j in same_list_FEEL:
                if i in j and i not in final and j[0] not in final:
                    final.append(j[0])
                    break
                elif i not in j and i not in final:
                    final.append(i)
                    break
            ####################################################################
            ########## CHECK if it works
            ##########################################################
            ####################################################################
    for i in final_SPAdes:
        ###
        if same_list_SPA == []:
            if i not in final:
                final.append(i)
        else:
        ###
            for j in same_list_SPA:
                if i in j and i not in final and j[0] not in final:
                    final.append(j[0])################build in i==j[0] or something
                    break
                elif i not in j and i not in final:
                    final.append(i)
                    break
#    print(len(final))
    with open(output,'w') as outfile:
        outfile.write('### putative overlapping lncRNA candidates, that were predicted by both pipelines\n')
        for i in final:
            outfile.write(str(i)+'\n')
#    print('used fastafile:',fasta_file,'\nfasta output:',output_fasta)
    with open(output_fasta, 'w') as outfile:
        for i in final:
            for k,v in fasta.items():
                if i == k:
                    outfile.write('>'+str(i)+'\n'+str(v)+'\n')
    return output_fasta

#check_same_final(final_feel, final_spa,'/home/debian9/host/Deepseq_data/lncRNA/lnc_late/rRNdepLTP_L4_R1_001_56r7DjFcMGSl.fastq'  , '/home/debian9/host/scripts/temp/combined_lnc_candidates.fasta' , same_feel, same_spa)

################################################################################
#                       under construction
################################################################################
def make_log(filename, minlen, dictFEEL, dictSPAdes, dictboth):
    print('\n\n\n\tDone!!\n\n####\n>>>generated items:\n\n   -overlapping candidates(id-list & fasta)\n   -all candidates by both pipelines(id-list & fasta)\
\n   -log-file containing some interesting specs\n>>> for this seq-data '+str(filename))
    filepath = 'result/log_lnc_pipeline_for_'+str(filename).split('/')[-1].split('_')[0]+'.txt'
    FEEL = {str(i.description):str(i.seq) for i in SeqIO.parse(str(dictFEEL), 'fasta')}
    SPAdes = {str(i.description):str(i.seq) for i in SeqIO.parse(str(dictSPAdes), 'fasta')}
    shared = {str(i.description):str(i.seq) for i in SeqIO.parse(str(dictboth), 'fasta')}
    with open(filepath,'w') as outfile:
        outfile.write('### this file serves as a log in order to check some parameters afterwards\n\n')
        outfile.write(' >>> Sequencing data path: '+str(filename)+'\n >>> Trimming: trimming quality 15 with a minlength of '+str(minlen)+'\n\
 >>> FEELnc only candidates: '+str(len(list(FEEL.keys())))+'\n >>> SPAdes only candidates: '+str(len(list(SPAdes.keys())))+'\n \
>>> shared candidates: '+ str(len(list(shared.keys()))))
    ## trimm parameters and reads that are left after timming
    ## initial FEELnc_candidates
    ## initial SPA_candidates
    ## filtered FEELnc_candidates
    ## SPA_blast
    ## filtered SPA_candidates
    ## candidates that overlap

def main():
    if len(sys.argv[1:]) != 6:
        print('>>> Too many or missing arguments!!\n    Usage:\n\n    lnc_pipeline.py\
 [sequencing_file] [minimum_length_for_read_trimming] [HiSat2_refference] [GFF_annotation_file]\
 [Fasta_refference_genome] [HMMer_database]')
        sys.exit()

    file = sys.argv[1]  ### sequencing file
    minlen = sys.argv[2]    ### minimum length of the reads
    ref_path = sys.argv[3]  ### hisat2 reference
    guide_ref = sys.argv[4] ### gff annotation file
    ref_anno = guide_ref
#    ref_anno = sys.argv[5]  ### theoreticaly same as [4] ;gff (or if available gtf) annotation file################(##############)###################
    genome = sys.argv[5]    ### genome used for mapping
    hmm_database = sys.argv[6]  ### hmmer database / pfam database
    blast_database_fasta = genome
#    blast_database_fasta = sys.argv[8]    ### theoreticaly same as [6]; fasta of sequences you want to build the blast database with

#    F_exit_cnt = 0
#    S_exit_cnt = 0

    check_folders()
    trim = trim_seq_data(file, minlen)
    bam = map_Hisat(trim, ref_path)
        #bam = map_Hisat(trim_seq_data(file, minlen), ref_path)
    sorted = filter_sort(bam)
    FEELnc_input = bamToGtf(sorted, guide_ref)
        #FEELnc_input = bamToGtf(filter_sort(sam), guide_ref)
    FEELnc_ref = gff_to_gtf(ref_anno)
    FEELnc_output = FEELnc(FEELnc_input, FEELnc_ref)
        #gtf_to_fasta(FEELnc(FEELnc_input, FEELnc_ref), genome))
    FEELnc_fasta = gtf_to_fasta(FEELnc_output, genome)
    f_exit = F_exit(FEELnc_fasta, F_exit_cnt)
    #print('\n'+str(f_exit)+'\n')
    sp = SPAdes(trim)
    sp_sorted, plus200bp = sort_S_output(sp)

    blastdb = make_blast_database(blast_database_fasta)
    blast_rep = blastn(sp_sorted, blastdb)
    #print('plus200pb:',len(list(plus200bp.keys())))
    SPAdes_fasta = filter_blastn(blast_rep, plus200bp)
    s_exit = S_exit(SPAdes_fasta, S_exit_cnt)
    #print('\n'+str(s_exit)+'\n')

    #sys_exit(F_exit_cnt, S_exit_cnt)
    if f_exit == 1:
        if s_exit == 0:
            fasta = SPAdes_fasta
        elif s_exit == 1:
            sys_exit(f_exit, s_exit)
    elif s_exit == 1:
        if f_exit == 0:
            fasta = FEELnc_fasta
        elif f_exit == 1:
            sys_exit(f_exit, s_exit)
    else:
        fasta = cat_fastas(FEELnc_fasta, SPAdes_fasta)

    cand6RF = all6RF(fasta)
    out, out_fasta = hmmsearch(hmm_database, cand6RF, fasta, file)
    feel_split, spa_split = split_cand(out_fasta)
    FEELnc_candidates_BLAST_DB = make_blast_database(feel_split)
    SPAdes_candidates_BLAST_DB = make_blast_database(spa_split)
    spa_vs_feelDB = blast_ungapped(spa_split, FEELnc_candidates_BLAST_DB)
    feel_vs_spaDB = blast_ungapped(feel_split, SPAdes_candidates_BLAST_DB)
    SPAvsFEEL_same_dict, same_SPAvsFEEL = check_same_fromOne(spa_vs_feelDB)
    FEELvsSPA_same_dict, same_FEELvsSPA = check_same_fromOne(feel_vs_spaDB)
    shared_candidates = check_same_final(FEELvsSPA_same_dict, SPAvsFEEL_same_dict, file, out_fasta, same_FEELvsSPA, same_SPAvsFEEL)
    make_log(file, minlen, feel_split, spa_split, shared_candidates)
    os.system('rm -r temp')


main()
