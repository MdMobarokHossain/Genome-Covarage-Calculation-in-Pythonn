from Bio import SeqIO
fastaFile = "EpiCoV_BulkUpload_Batch_27&26_30thSep2021_study.fa"
handle = open(fastaFile,"r")
for record in SeqIO.parse(handle,"fasta"):
    if record.id=="hCoV-19/Bangladesh/CHRF-0631/2021":
        # print('>'+record.id+'\n'+record.seq)
        x=len(record.seq)
        y=record.seq.count("N")
        m=record.seq.count("G")+record.seq.count("C")+record.seq.count("A")+record.seq.count("T")
        gc= ((record.seq.count("G")+record.seq.count("C"))/(record.seq.count("G")+record.seq.count("C")+record.seq.count("A")+record.seq.count("T")))*100
        print ('The GC content of provided sequence:',round(gc,2),"%")
        print ('The number of N:' ,y)
        print ('The genome coverage:', round((m-y)/m*100,2),"%")