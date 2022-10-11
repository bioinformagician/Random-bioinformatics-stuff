def fasta_chopper_gc_extractor(fasta_txt_file):
    with open(fasta_txt_file, 'r') as f: #open the fastafile and put it in a list where each index corresponds to a row
        myNames = [line.strip() for line in f]
    ids=[]
    sequences=[]
    arrow_counter=0
    index_saver=[]
    for i in range(len(myNames)): #gets fasta IDs and their index to isolate sequences
        if myNames[i][0]=='>':
            ids.append(myNames[i])
            index_saver.append(i)

    for i in range(len(index_saver)-1): #isolates sequences
        sequences.append(myNames[index_saver[i]+1:index_saver[i+1]])
    sequences.append(myNames[index_saver[-1]+1:len(myNames)])
    for i in range(len(sequences)):
        sequences[i]=''.join(sequences[i])
    currentwinner = 0
    highest_gc = 0
    g_counter = 0
    c_counter = 0
    for row in range(len(sequences)): #find highest gc content
        for pos in range(len(sequences[row])):
            if sequences[row][pos] == 'G':
                g_counter += 1
            elif sequences[row][pos] == 'C':
                c_counter += 1
        current_gc = (g_counter + c_counter) / len(sequences[row])
        if current_gc > highest_gc:
            highest_gc = current_gc
            currentwinner = row
        g_counter = 0
        c_counter = 0
    thelist = [currentwinner, highest_gc]
    print(ids[thelist[0]], 'with a GC content of',thelist[1])
