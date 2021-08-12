#!/usr/bin/python3
import os
import argparse

def parse_hits(lines,skip_first):
    if len(lines)<=4:
        return None
    ##Query Id
    record_result = list()
    query_id = lines[0].split(":")[1].strip()
    record_result.append(query_id)
    
    ## skip to Hit table
    line_p = 1
    while line_p < len(lines):
        if lines[line_p].endswith("hits found"):
            break
        else:
            line_p+=1
    line_p += ( 2 if skip_first else 1)

    identity_sum = 0
    bitscore_sum = 0
    positives_sum=0
    score_sum = 0
    if line_p == len(lines):
        return None
    hit_nums=0
    while line_p < len(lines):
        line = lines[line_p]
        pid,bitscore,ppos,score = line.rstrip().split("\t")
        hit_nums+=1
        line_p+=1
        identity_sum += float(pid)
        bitscore_sum += float(bitscore)
        positives_sum += float(ppos)
        score_sum += float(score)
    record_result.append(f"{round(identity_sum/hit_nums,4)},{round(bitscore_sum/hit_nums,4)},{round(positives_sum/hit_nums,4)},{round(score_sum/hit_nums,4)}")
    return ",".join(record_result)

def parse_blastp(filename,skip_first,filename_out):
    if not os.path.exists(filename):
        raise FileNotFoundError(filename)
        
    handler = open(filename,'r')
    dirname = os.path.dirname(filename)
    basename = os.path.basename(filename)
    

    output_name = filename_out if filename_out!= None else os.path.join(dirname,"blastp_processed_"+basename.split(".")[0]+".csv")
    
    out = open(output_name,'w+')
    
    out.write("query_id,")
    out.write("% Identity,bitscore,% Positives,score\n")
    counter = 0
    lines = list()
    while True:
        try:
            line = next(handler)
            if line.startswith("# BLASTP"):
                parse_result = parse_hits(lines,skip_first)
                if parse_result != None:
                    out.write(parse_result+"\n")
                    counter += 1
                    if (counter%5000)==0:
                        print(f'{counter} processed')
                lines = list()
            else:
                lines.append(line.rstrip())
        except StopIteration:
            parse_result = parse_hits(lines[:-1],skip_first)
            if parse_result != None:
                out.write(parse_result+"\n")
                counter += 1
                print(f'{counter} processed')
            break
    handler.close()
    out.close()

def main():
        parser = argparse.ArgumentParser(description='blastp result processor',add_help=True)
        parser.add_argument('--input','-i',help='result file from blastp(outfmt "7 pident bitscore ppos score")')
        parser.add_argument('-skip',help='1 if skip the first hit, in default 0',default=0)
        parser.add_argument('--out',help='The output file path')
        args = parser.parse_args()
        fin = args.input
        fout = args.out
        skip_first=args.skip
        if fin == None:
            parser.print_help()
            print("Input file not found")
            exit(1)
        parse_blastp(fin,skip_first,fout)
    
main()
