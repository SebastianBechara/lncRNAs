rep = open('/home/debian9/host/scripts/result/lnc_all_TP_new_20190811/final_lnc_cand_vs_MAC_rep').readlines()
lnc_names = open('/home/debian9/host/scripts/result/lnc_all_TP_new_20190811/lncRNA_Renamed.txt').readlines()
gff =[]

def get_lnc_rep_lines(rep, lnc_names):
    rep_dict={}
    for i in rep:
        aln_length = int(i.split('\t')[3])
        for j in lnc_names:
            if i.split('\t')[0] == j.split('\t')[1].strip():
                #if j.startswith('FEE'):
                #    key = i.split('\t')[0]+'_'+i.split('\t')[1]+'_'+inr(j.split('\t')[0].split(''))
                #else:
                key = i.split('\t')[0]+'_'+i.split('\t')[1]+'_'+j.split('\t')[0].split('_length_')[-1].split('_')[0]
                rep_dict.setdefault(key, []).append(aln_length)
                ################################################################
                # figure out why this didn't work !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ################################################################
    rep_lines={}
    for i in rep:
        for j in lnc_names:
            rep_lines.setdefault(i.split('\t')[0]+'_'+i.split('\t')[1]+'_'+j.split('\t')[0].split('_length_')[-1].split('_')[0], []).append(i)

    good_rep_lines = {}
    for k,v in rep_dict.items():
        for k1,v1 in rep_lines.items():
            if k == k1:
                if sum(v) / int(k.split('_')[-1]) >= 0.90:
                    good_rep_lines.setdefault(k, []).append(v1)
    return rep_dict

def build_gtf(good_rep_lines):
    gtf_lines =[]
    for k,v in good_rep_lines.items():
        for i in v:
            for j in i:
                scaff = j.split('\t')[1]
                source_prog = 'FEELnc_SPAdes_pipeline'
                feature = 'lncRNA'
                _start = int(j.split('\t')[8])
                _end = int(j.split('\t')[9])
                start = str(min(_start, _end))
                end = str(max(_start, _end))
                score = '.'
                if _start > _end:
                    orientation = '-'
                else:
                    orientation = '+'
                frame = '.'
                end_thing = 'id '+j.split('\t')[0]+'; coordinates '+'..'.join(j.split('\t')[6:8])+'; '+'biotype: putative lncRNA;'
                gtf_lines.append('\t'.join([scaff, source_prog, feature, start, end, score, orientation, frame, end_thing]))
    return gtf_lines

def write_gtf(gtf_lines):
    with open('lnc_candidates_anno.gtf', 'w') as outfile:
        for i in gtf_lines:
            outfile.write(i+'\n')
