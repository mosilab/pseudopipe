#!/usr/bin/env python3

import sys
import csv
from Bio import SeqIO


DFAST_gbk_dir = sys.argv[1]
prokka_gbk_dir = sys.argv[2]
PF_gbff = sys.argv[3]
prokka_pseudos = sys.argv[4]
output_dir = sys.argv[5]
# output_dir2 = sys.argv[6] # future additional input for different pipelines
predictor_type = sys.argv[6]


cd = output_dir
sys.stdout = open(f"{cd}/log.txt", "w")

########## Processing gbk files into csv ################
def gbk2csv(seq):
    
    print("Generating CSV file with all genome features from Genebank file \n")

    annots = SeqIO.index(seq, "genbank")
    keys_=set(())
    for annot in annots:
        for feature in annots[f'{annot}'].features:
            for dic in feature.qualifiers.keys():
                keys_.add(dic)
    keys_list = ['type','location']+list(keys_)
    fieldnames = keys_list
    
    with open(f"{seq}_genome.csv", "w") as output:
        writer = csv.DictWriter(output, fieldnames, restval='NA')
        writer.writeheader() 

        for annot in annots:
            for feature in annots[f'{annot}'].features:
    #             if feature.type=='CDS':
                row_data = {'type':f'{feature.type}', 'location':f'{feature.location}'}
                for key,value in feature.qualifiers.items():
                    if len(value)>1:
                        row_data[key] = value[0:len(value)]
                    else:
                        row_data[key] = value[0]
    #                 print(len(value))
                writer.writerow(row_data)  
    output.close()
    
    data = pd.read_csv(f'{seq}_genome.csv')
    return data

def gbk2csv_cds(seq):
    
    print("Generating CSV file with CDS from Genebank file \n")

    annots = SeqIO.index(seq, "genbank")
    keys_=set(())
    for annot in annots:
        for feature in annots[f'{annot}'].features:
            for dic in feature.qualifiers.keys():
                keys_.add(dic)
    keys_list = ['type','location']+list(keys_)
    fieldnames = keys_list
    
    with open(f"{seq}_cds.csv", "w") as output:
        writer = csv.DictWriter(output, fieldnames, restval='NA')
        writer.writeheader() 

        for annot in annots:
            for feature in annots[f'{annot}'].features:
                if feature.type=='CDS':
                    row_data = {'type':f'{feature.type}', 'location':f'{feature.location}'}
                    for key,value in feature.qualifiers.items():
                        if len(value)>1:
                            row_data[key] = value[0:len(value)]
                        else:
                            row_data[key] = value[0]
        #                 print(len(value))
                    writer.writerow(row_data)  
    output.close()
    
    data = pd.read_csv(f'{seq}_cds.csv')
    return data



####################### DFAST ########################
############# Getting DFAST pseudogenes ##############
import pandas as pd

#use gbk2csv before this and provide cds.csv path

def get_DFAST_pseudos(file_path):
    
    print("Getting DFAST CDS and pseudos in to fasta \n")
    # cds = pd.read_csv({file_path})
    # cds
    cds = file_path
    
    nonsense = cds.loc[cds['note'].str.contains("stop", case=False, na=False)]
    frameshift = cds.loc[cds['note'].str.contains("frameshift", case=False, na=False)]

    cds["mutation"]='none'
    cds.loc[nonsense.index, "mutation"] = "nonsense"
    cds.loc[frameshift.index, "mutation"] = "frameshift"
    # print(cds)
    # cds
    if "transl_except" in cds.keys():
        print("transl_except present in column \n")
        DFAST_pseuds = cds[(cds["mutation"] != 'none') | (cds["transl_except"].notnull())]
    # DFAST_pseuds
    else:
        print("transl_except absent in column \n")
        DFAST_pseuds = cds.loc[cds["mutation"] != 'none']
    # print(DFAST_pseuds)
    print (cds.mutation.value_counts())
    DFAST_pseuds.to_csv(f'{cd}/DFAST_pseudos.csv')
    
    print("Writing all DFAST CDS to fasta file \n")
    counter = 0
    with open(f"{cd}/DFAST_cdsAA.fasta", "w") as output:
        for _,series in cds.iterrows():
    #         print(key,value[0])
            if not pd.isna(series.translation):
                output.write(f'>{series.locus_tag} {series.location}\n{series.translation}\n')
                counter += 1
    output.close()
    print(f'{counter} DFAST CDS written to fasta \n')

    print("Writing all DFAST pseudos CDS to fasta file \n")
    counterDFAST = 0
    with open(f"{cd}/DFAST_pseudos_cdsAA.fasta", "w") as output:
        for _,series in DFAST_pseuds.iterrows():
    #         print(key,value[0])
            if not pd.isna(series.translation):
                output.write(f'>{series.locus_tag} {series.location}\n{series.translation}\n')
                counterDFAST += 1
    output.close()
    print(f'{counterDFAST} DFAST pseudos written to fasta \n')

    return counterDFAST


###################################### PSEUDOFINDER ############################################
################# Get Pseudofinder old_locus_tags in pseudos gbff file ##########################
import os
#import gffutils

#use PF gbff file
def get_PF_old_loc_tags(file_path):
    print(file_path)
    print("Getting pseudofinder locus tags \n")
    # cd = os.getcwd()
    # cd

    gff = file_path
    # print(gff)
    #fn = gffutils.example_filename(f'{gff}')
    # print(fn)
    # print(open(fn).read().split('\n'))
    fn = gff

    pseuds=[]
    for item in open(fn).read().split('\n'):
        if not "#" in item:
            pseuds.append(item.split("\t"))
    # print(pseuds)

    columns = ['seqname', 'source', 'feature', 'start', 'end', 'mark1', 'strand', 'mark2', 'attribute']
    df = pd.DataFrame(pseuds,columns=columns)

    df[['note','locus_tag','old_locus_tag']] = df.attribute.str.split(';',expand=True)
    # df
    
    #TO-DO split reason in note by fullstop and remove trailing empties in list
    
    df['old_locus_tag'] = df['old_locus_tag'].str.replace(r'old_locus_tag=', '')
    df.drop(df.tail(1).index,inplace=True)
    # df.head(10)

    old_loc_list = []

    for item in df.old_locus_tag.tolist():
        old_loc_list.append(item.split(','))
    # print(old_loc_list)

    old_loc_len = []
    counter = 0
    for i in old_loc_list:
        old_loc_len.append(len(i))
        counter += 1
    # print(old_loc_len)
    print(f'{counter} lines of PF old location tags found')

    df = pd.DataFrame(list(zip(old_loc_list, old_loc_len)),
                   columns =['old_loc_tags', 'number_of_loc_tags'])
    # df

    df.to_csv(f'{cd}/PF_old_loc_tags.csv')
    
    return old_loc_list

################# Get pseudofinder pseudos ###########################
# import csv
# from Bio import SeqIO

def get_PF_pseudos(old_loc_tags,file_path):
    
    print("Getting pseudofinder pseudos CDS \n")
    PF_pseudos = old_loc_tags
    len(PF_pseudos)

    # data = pd.read_csv({file_path})
    data = file_path

    seen = {}
    not_found = {}
    for item in PF_pseudos:
    #     for loc_tag in item:
    #     seen[f'{item}']=f'{data.loc[data.locus_tag==item[0]]}'
    #     print(item,'\n',data.loc[data.locus_tag==item[0]].index)
    #     data.loc[data.locus_tag==item[0]][translation]
        no_hit=[]
        if len(item)>1:
            for i in range(len(item)):
                found = data[data.locus_tag==item[i].strip()].translation.tolist()
                if len(found)!=0:
                    seen[f'{item[i].strip()} {data[data.locus_tag==item[i].strip()].location.tolist()[0]}']=found
                    break
    #             print(found)
                else:
                    no_hit.append(item[i].strip())
    #                 not_found[f'{item[i]}']=found ## append to a list and finally write it to the dict if nothing is found
                not_found[f'{no_hit}']=[]
        else:
            found = data[data.locus_tag==item[0].strip()].translation.tolist()
            if len(found)!=0:
                seen[f'{item[0].strip()} {data[data.locus_tag==item[0].strip()].location.tolist()[0]}']=found
    #         seen[f'{item}']= f'{found}'
    #         print(found)
            else:
                not_found[f'{item[0].strip()}']=found
    print("Pseudogenes found ", len(seen.keys()))
    print("Pseudogenes not found ", len(not_found.keys()))

    counterPF = 0
    with open(f"{cd}/PF_pseudos_cdsAA.fasta", "w") as output:
        for key,value in seen.items():
            print(key,value[0])
            key = key.split("_")
            key = key[0]+"_PF_"+key[1]
            output.write(f'>{key}\n{value[1]}\n')
            counterPF+=1
    output.close()
    print(counterPF, " total pseudogenes written into fasta \n")
    
    return counterPF


##################### PROKKA #############################
############ Get prokka pseudogenes ######################
#use tsv
output_path = cd

def get_prokka_cds(file_path):

    print("Getting Prokka CDS \n")

    cds = file_path
    
    with open(f"{cd}/prokka_cdsAA.fasta", "w") as output:
        for _,series in cds.iterrows():
    #         print(key,value[0])
            if not pd.isna(series.translation):
                output.write(f'>{series.locus_tag} {series.location}\n{series.translation}\n')
    output.close()


def get_prokka_pseudos(file_path,reference_gbkcsv):
    
    print("Getting Prokka pseudos \n")
    df = pd.read_csv(file_path, sep='\t', skiprows=1, skipfooter=1, engine='python')
#     df

    prokka_pseudos = df['Start'].tolist()
#     prokka_pseudos

    # data = pd.read_csv(reference_gbkcsv)
    data = reference_gbkcsv
#     data

    seen = {}

    for item in prokka_pseudos:
        print(item)
        found = [x for x in data[data.locus_tag==item.strip()].translation.tolist() if str(x) != "nan"]
        print(f'{found} found')
        if found:
            seen[f'{item.strip()} {data[data.locus_tag==item.strip()].location.tolist()[0]}']=found

#     print(seen)

    counterPK = 0
    with open(f"{output_path}/prokka_pseudos_cdsAA.fasta", "w") as output:
        for key,value in seen.items():
            # print(key,value[0])
            output.write(f'>{key}\n{value[0]}\n')
            counterPK+=1
    output.close()
    print(f'{counterPK} Prokka pseudos found \n')

    return counterPK
    

################### Combining Pseudogene AA sequences ######################

def combined_pseudos_writer(pseudo_filename,outfile):
    # Open each file in read mode
    with open(f'{output_path}/{pseudo_filename}') as infile:

        # read the data from file1 and
        # file2 and write it in file3
        outfile.write(infile.read())

    # Add '\n' to enter data of file2
    # from next line
    outfile.write("\n")

def combining_all_pseudos(predictor_type='all'):
    
    print("Combining pseudogenes into single fasta file \n")
    all_filenames = ["PF_pseudos_cdsAA.fasta","prokka_pseudos_cdsAA.fasta","DFAST_pseudos_cdsAA.fasta"]

    # Open file3 in write mode
    with open(f'{output_path}/combined_pseudos.fasta', 'w') as outfile:

        if predictor_type == 'all':
            # Iterate through list
            for names in all_filenames:
                combined_pseudos_writer(names,outfile)
        elif predictor_type == 'PF':
            combined_pseudos_writer(all_filenames[0],outfile)
        elif predictor_type == 'prokka':
            combined_pseudos_writer(all_filenames[1],outfile)
        elif predictor_type == 'DFAST':
            combined_pseudos_writer(all_filenames[2],outfile)
            
if predictor_type == 'all':
    DFAST_gbkcsv = gbk2csv(DFAST_gbk_dir)
    prokka_gbkcsv = gbk2csv(prokka_gbk_dir)
    prokka_gbkcsv_cds = gbk2csv_cds(prokka_gbk_dir)
    DFAST_gbkcsv_cds = gbk2csv_cds(DFAST_gbk_dir)

    counterDFAST = get_DFAST_pseudos(DFAST_gbkcsv_cds)
    # # get_DFAST_pseudos(f"{DFAST_gbk_dir}_genome.csv")
    get_prokka_cds(prokka_gbkcsv_cds)


    old_loc_tags = get_PF_old_loc_tags(PF_gbff) #working with a few tweaks
    counterPF = get_PF_pseudos(old_loc_tags,prokka_gbkcsv)

    counterPK = get_prokka_pseudos(prokka_pseudos,prokka_gbkcsv) #working with a few tweaks

    combining_all_pseudos()

    counter = counterDFAST + counterPF + counterPK        
    print(f'{counter} pseudogenes written to combined fasta file \n')

elif predictor_type == 'PF':
    prokka_gbkcsv = gbk2csv(prokka_gbk_dir)
    prokka_gbkcsv_cds = gbk2csv_cds(prokka_gbk_dir)

    get_prokka_cds(prokka_gbkcsv_cds)
    old_loc_tags = get_PF_old_loc_tags(PF_gbff) #working with a few tweaks
    counterPF = get_PF_pseudos(old_loc_tags,prokka_gbkcsv)
    print(f'{counterPF} pseudogenes written to combined fasta file \n')
    combining_all_pseudos('PF')

elif predictor_type == 'prokka':
    prokka_gbkcsv = gbk2csv(prokka_gbk_dir)
    prokka_gbkcsv_cds = gbk2csv_cds(prokka_gbk_dir)
    get_prokka_cds(prokka_gbkcsv_cds)
    counterPK = get_prokka_pseudos(prokka_pseudos,prokka_gbkcsv) #working with a few tweaks    
    print(f'{counterPK} pseudogenes written to combined fasta file \n')
    combining_all_pseudos('prokka')

elif predictor_type == 'DFAST':
    DFAST_gbkcsv = gbk2csv(DFAST_gbk_dir)
    DFAST_gbkcsv_cds = gbk2csv_cds(DFAST_gbk_dir)
    counterDFAST = get_DFAST_pseudos(DFAST_gbkcsv_cds)
    print(f'{counterDFAST} pseudogenes written to combined fasta file \n')
    combining_all_pseudos('DFAST')