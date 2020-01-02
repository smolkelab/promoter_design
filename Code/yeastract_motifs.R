
find.yeastract = function(motifs.pgpd) {

  motifs.pgpd$Mcm1p_F = grepl('TTACC.AATT.GGTAA',motifs.pgpd$Seq) | grepl('TT[AT]CC[CT]AA[AT]..GG[AT]AA[AT][AT]',motifs.pgpd$Seq) | grepl('[AGT]CC[CT][AT][AT][AT]..[AG]G',motifs.pgpd$Seq) | grepl('CC[CT][AT][AT][AT]..[AG]G',motifs.pgpd$Seq) | grepl('CCTAATTAGG',motifs.pgpd$Seq)
  motifs.pgpd$Mcm1p_R = grepl('TTACC.AATT.GGTAA',motifs.pgpd$Seq) | grepl('[AT][AT]TT[AT]CC..[AT]TT[AG]GG[AT]AA',motifs.pgpd$Seq) | grepl('C[CT]..[AT][AT][AT][AG]GG[ACT]',motifs.pgpd$Seq) | grepl('C[CT]..[AT][AT][AT][AG]GG',motifs.pgpd$Seq) | grepl('CCTAATTAGG',motifs.pgpd$Seq)
  motifs.pgpd$Usv1p_F = grepl('AGGGG',motifs.pgpd$Seq) | grepl('TCAGGGGT',motifs.pgpd$Seq) | grepl('TCAGGTAC',motifs.pgpd$Seq)
  motifs.pgpd$Usv1p_R = grepl('CCCCT',motifs.pgpd$Seq) | grepl('ACCCCTGA',motifs.pgpd$Seq) | grepl('GTACCTGA',motifs.pgpd$Seq)
  motifs.pgpd$Mig1p_F = grepl('AATT[AG]TCCGGGG',motifs.pgpd$Seq) | grepl('TATAGTGCGGGG',motifs.pgpd$Seq) | grepl('ATTTTTGTGGGG',motifs.pgpd$Seq) | grepl('GAATA[CT]CTGGGG',motifs.pgpd$Seq) | grepl('TTATTTCTGGGG',motifs.pgpd$Seq) | grepl('TAAAAGCCGGGG',motifs.pgpd$Seq) | grepl('TTAAAAGCGGGG',motifs.pgpd$Seq) | grepl('TTAAACGTGGGG',motifs.pgpd$Seq) | grepl('ATAAAACTGGGG',motifs.pgpd$Seq) | grepl('ATAAACGCGGGG',motifs.pgpd$Seq) | grepl('ATTAATGTGGGG',motifs.pgpd$Seq) | grepl('AAAAATGCGGGG',motifs.pgpd$Seq) | grepl('ATTTTGCGGGG',motifs.pgpd$Seq) | grepl('CCCC[AG]C',motifs.pgpd$Seq) | grepl('[GC][CT]GG[AG]G',motifs.pgpd$Seq)
  motifs.pgpd$Mig1p_R = grepl('CCCCGGA[CT]AATT',motifs.pgpd$Seq) | grepl('CCCCGCACTATA',motifs.pgpd$Seq) | grepl('CCCCACAAAAAT',motifs.pgpd$Seq) | grepl('CCCCAG[AG]TATTC',motifs.pgpd$Seq) | grepl('CCCCAGAAATAA',motifs.pgpd$Seq) | grepl('CCCCGGCTTTTA',motifs.pgpd$Seq) | grepl('CCCCGCTTTTAA',motifs.pgpd$Seq) | grepl('CCCCACGTTTAA',motifs.pgpd$Seq) | grepl('CCCCAGTTTTAT',motifs.pgpd$Seq) | grepl('CCCCGCGTTTAT',motifs.pgpd$Seq) | grepl('CCCCACATTAAT',motifs.pgpd$Seq) | grepl('CCCCGCATTTTT',motifs.pgpd$Seq) | grepl('CCCCGCAAAAT',motifs.pgpd$Seq) | grepl('G[CT]GGGG',motifs.pgpd$Seq) | grepl('C[CT]CC[AG][GC]',motifs.pgpd$Seq)
  motifs.pgpd$Nrg2p_F = grepl('[AC]AGGGTCC',motifs.pgpd$Seq)
  motifs.pgpd$Nrg2p_R = grepl('GGACCCT[GT]',motifs.pgpd$Seq)
  motifs.pgpd$Dal81p_F = grepl('[GC]AAA[AT].TGCG[CGT]T',motifs.pgpd$Seq) | grepl('GAAAATTGCGTT',motifs.pgpd$Seq) | grepl('CGGC......GCCG',motifs.pgpd$Seq)
  motifs.pgpd$Dal81p_R = grepl('A[ACG]CGCA.[AT]TTT[GC]',motifs.pgpd$Seq) | grepl('AACGCAATTTTC',motifs.pgpd$Seq) | grepl('CGGC......GCCG',motifs.pgpd$Seq)
  motifs.pgpd$Msn4p_F = grepl('CCCCT',motifs.pgpd$Seq) | grepl('AGGGG',motifs.pgpd$Seq)
  motifs.pgpd$Msn4p_R = grepl('AGGGG',motifs.pgpd$Seq) | grepl('CCCCT',motifs.pgpd$Seq)
  motifs.pgpd$Rox1p_F = grepl('TCTATTGTTTCCC',motifs.pgpd$Seq) | grepl('[CT][CT][CT]ATTGTTCTC',motifs.pgpd$Seq) | grepl('ACAAT',motifs.pgpd$Seq)
  motifs.pgpd$Rox1p_R = grepl('GGGAAACAATAGA',motifs.pgpd$Seq) | grepl('GAGAACAAT[AG][AG][AG]',motifs.pgpd$Seq) | grepl('ATTGT',motifs.pgpd$Seq)
  motifs.pgpd$Cup2p_F = grepl('[ACT]T[ACT]..GCTG[AGT]',motifs.pgpd$Seq) | grepl('GCGTCTTTTCCGCTGA',motifs.pgpd$Seq) | grepl('TCTTTTGCTG',motifs.pgpd$Seq) | grepl('TCTTTTTGCTG',motifs.pgpd$Seq) | grepl('TCTTTTTTGCTG',motifs.pgpd$Seq)
  motifs.pgpd$Cup2p_R = grepl('[ACT]CAGC..[AGT]A[AGT]',motifs.pgpd$Seq) | grepl('TCAGCGGAAAAGACGC',motifs.pgpd$Seq) | grepl('CAGCAAAAGA',motifs.pgpd$Seq) | grepl('CAGCAAAAAGA',motifs.pgpd$Seq) | grepl('CAGCAAAAAAGA',motifs.pgpd$Seq)
  motifs.pgpd$Gzf3p_F = grepl('GATAAG',motifs.pgpd$Seq)
  motifs.pgpd$Gzf3p_R = grepl('CTTATC',motifs.pgpd$Seq)
  motifs.pgpd$Sip4p_F = grepl('TCCATT[GC][AG]TCCG[AG]',motifs.pgpd$Seq) | grepl('.CC[AGT]T[CT].[ACG].CCG.',motifs.pgpd$Seq) | grepl('[CT]CC.[CT]T.[AG][AG]CCG.',motifs.pgpd$Seq)
  motifs.pgpd$Sip4p_R = grepl('[CT]CGGA[CT][GC]AATGGA',motifs.pgpd$Seq) | grepl('.CGG.[CGT].[AG]A[ACT]GG.',motifs.pgpd$Seq) | grepl('.CGG[CT][CT].A[AG].GG[AG]',motifs.pgpd$Seq)
  motifs.pgpd$Hac1p_F = grepl('CAGCGTG',motifs.pgpd$Seq) | grepl('ATGGTATCAT',motifs.pgpd$Seq) | grepl('TGACGTCA',motifs.pgpd$Seq) | grepl('CCAGC',motifs.pgpd$Seq)
  motifs.pgpd$Hac1p_R = grepl('CACGCTG',motifs.pgpd$Seq) | grepl('ATGATACCAT',motifs.pgpd$Seq) | grepl('TGACGTCA',motifs.pgpd$Seq) | grepl('GCTGG',motifs.pgpd$Seq)
  motifs.pgpd$Aca1p_F = grepl('TGACGTCA',motifs.pgpd$Seq) | grepl('TTACGTAA',motifs.pgpd$Seq)
  motifs.pgpd$Aca1p_R = grepl('TGACGTCA',motifs.pgpd$Seq) | grepl('TTACGTAA',motifs.pgpd$Seq)
  motifs.pgpd$Leu3p_F = grepl('CCG....CGG',motifs.pgpd$Seq) | grepl('CCGG..CCGG',motifs.pgpd$Seq) | grepl('CCGTTAACGG',motifs.pgpd$Seq)
  motifs.pgpd$Leu3p_R = grepl('CCG....CGG',motifs.pgpd$Seq) | grepl('CCGG..CCGG',motifs.pgpd$Seq) | grepl('CCGTTAACGG',motifs.pgpd$Seq)
  motifs.pgpd$Rfx1p_F = grepl('TCGCCATGGCAAC',motifs.pgpd$Seq) | grepl('TTGCCATGGCAAC',motifs.pgpd$Seq) | grepl('TCGCCATGACAAC',motifs.pgpd$Seq) | grepl('CTATTGCTGCAAC',motifs.pgpd$Seq) | grepl('TCTCTGTGGCAAC',motifs.pgpd$Seq) | grepl('TTGTCACAGCAAC',motifs.pgpd$Seq) | grepl('TTGTCGCAGCAAC',motifs.pgpd$Seq) | grepl('TTGTTGTGGCAAC',motifs.pgpd$Seq) | grepl('TTTCCACAGCAAC',motifs.pgpd$Seq) | grepl('T[CT]GCCATGCCAAC',motifs.pgpd$Seq)
  motifs.pgpd$Rfx1p_R = grepl('GTTGCCATGGCGA',motifs.pgpd$Seq) | grepl('GTTGCCATGGCAA',motifs.pgpd$Seq) | grepl('GTTGTCATGGCGA',motifs.pgpd$Seq) | grepl('GTTGCAGCAATAG',motifs.pgpd$Seq) | grepl('GTTGCCACAGAGA',motifs.pgpd$Seq) | grepl('GTTGCTGTGACAA',motifs.pgpd$Seq) | grepl('GTTGCTGCGACAA',motifs.pgpd$Seq) | grepl('GTTGCCACAACAA',motifs.pgpd$Seq) | grepl('GTTGCTGTGGAAA',motifs.pgpd$Seq) | grepl('GTTGGCATGGC[AG]A',motifs.pgpd$Seq)
  motifs.pgpd$Stp1p_F = grepl('CGGCTC',motifs.pgpd$Seq) | grepl('CGGC......CGGC',motifs.pgpd$Seq) | grepl('[AG][CT][AG]CGGC[AG]C',motifs.pgpd$Seq)
  motifs.pgpd$Stp1p_R = grepl('GAGCCG',motifs.pgpd$Seq) | grepl('GCCG......GCCG',motifs.pgpd$Seq) | grepl('G[CT]GCCG[CT][AG][CT]',motifs.pgpd$Seq)
  motifs.pgpd$Sko1p_F = grepl('TGACGTTT',motifs.pgpd$Seq) | grepl('TGACGTCA',motifs.pgpd$Seq) | grepl('TTACGTAA',motifs.pgpd$Seq) | grepl('ATGACGTA',motifs.pgpd$Seq)
  motifs.pgpd$Sko1p_R = grepl('AAACGTCA',motifs.pgpd$Seq) | grepl('TGACGTCA',motifs.pgpd$Seq) | grepl('TTACGTAA',motifs.pgpd$Seq) | grepl('TACGTCAT',motifs.pgpd$Seq)
  motifs.pgpd$Pdr8p_F = grepl('TCCG[ACT]GGA',motifs.pgpd$Seq)
  motifs.pgpd$Pdr8p_R = grepl('TCC[AGT]CGGA',motifs.pgpd$Seq)
  motifs.pgpd$Lys14p_F = grepl('TCC[AG].[CT]GGA',motifs.pgpd$Seq)
  motifs.pgpd$Lys14p_R = grepl('TCC[AG].[CT]GGA',motifs.pgpd$Seq)
  motifs.pgpd$Crz1p_F = grepl('CACCAGTCGGTGGCTGTGCGCTTG',motifs.pgpd$Seq) | grepl('GAATGGCTG',motifs.pgpd$Seq) | grepl('GGGTGGCTG',motifs.pgpd$Seq) | grepl('G.GGC[GT]CA',motifs.pgpd$Seq) | grepl('GGGGCTG',motifs.pgpd$Seq) | grepl('GAGCCC',motifs.pgpd$Seq) | grepl('CAGCCAC',motifs.pgpd$Seq)
  motifs.pgpd$Crz1p_R = grepl('CAAGCGCACAGCCACCGACTGGTG',motifs.pgpd$Seq) | grepl('CAGCCATTC',motifs.pgpd$Seq) | grepl('CAGCCACCC',motifs.pgpd$Seq) | grepl('TG[AC]GCC.C',motifs.pgpd$Seq) | grepl('CAGCCCC',motifs.pgpd$Seq) | grepl('GGGCTC',motifs.pgpd$Seq) | grepl('GTGGCTG',motifs.pgpd$Seq)
  motifs.pgpd$Ndt80p_F = grepl('G.C[AG]CAAA[AT]',motifs.pgpd$Seq) | grepl('CACAAA',motifs.pgpd$Seq)
  motifs.pgpd$Ndt80p_R = grepl('[AT]TTTG[CT]G.C',motifs.pgpd$Seq) | grepl('TTTGTG',motifs.pgpd$Seq)
  motifs.pgpd$Hap5p_F = grepl('T.ATTGGT',motifs.pgpd$Seq)
  motifs.pgpd$Hap5p_R = grepl('ACCAAT.A',motifs.pgpd$Seq)
  motifs.pgpd$Dal82p_F = grepl('[GC]AAA[AT].TGCG[CGT]T',motifs.pgpd$Seq) | grepl('GAAAATTGCGTT',motifs.pgpd$Seq)
  motifs.pgpd$Dal82p_R = grepl('A[ACG]CGCA.[AT]TTT[GC]',motifs.pgpd$Seq) | grepl('AACGCAATTTTC',motifs.pgpd$Seq)
  motifs.pgpd$Rtg3p_F = grepl('GTCAC',motifs.pgpd$Seq) | grepl('GGTAC',motifs.pgpd$Seq)
  motifs.pgpd$Rtg3p_R = grepl('GTGAC',motifs.pgpd$Seq) | grepl('GTACC',motifs.pgpd$Seq)
  motifs.pgpd$Gis1p_F = grepl('T[AT]AGGGAT',motifs.pgpd$Seq) | grepl('AGGGG',motifs.pgpd$Seq)
  motifs.pgpd$Gis1p_R = grepl('ATCCCT[AT]A',motifs.pgpd$Seq) | grepl('CCCCT',motifs.pgpd$Seq)
  motifs.pgpd$Oaf1p_F = grepl('CGG...TAA',motifs.pgpd$Seq)
  motifs.pgpd$Oaf1p_R = grepl('TTA...CCG',motifs.pgpd$Seq)
  motifs.pgpd$Kar4p_F = grepl('[CT]AA..CAAA..C.G.[CT]T',motifs.pgpd$Seq) | grepl('CTAA..CAAA..C.G.[CT]T',motifs.pgpd$Seq) | grepl('TTAA..CAAA..C.G.[CT]T',motifs.pgpd$Seq)
  motifs.pgpd$Kar4p_R = grepl('A[AG].C.G..TTTG..TT[AG]',motifs.pgpd$Seq) | grepl('A[AG].C.G..TTTG..TTAG',motifs.pgpd$Seq) | grepl('A[AG].C.G..TTTG..TTAA',motifs.pgpd$Seq)
  motifs.pgpd$Hap4p_F = grepl('T.ATTGGT',motifs.pgpd$Seq)
  motifs.pgpd$Hap4p_R = grepl('ACCAAT.A',motifs.pgpd$Seq)
  motifs.pgpd$Mac1p_F = grepl('TTTGC[GT]C[AG]',motifs.pgpd$Seq) | grepl('[AT][AT][AT]TTTGCTC[AG]',motifs.pgpd$Seq)
  motifs.pgpd$Mac1p_R = grepl('[CT]G[AC]GCAAA',motifs.pgpd$Seq) | grepl('[CT]GAGCAAA[AT][AT][AT]',motifs.pgpd$Seq)
  motifs.pgpd$Upc2p_F = grepl('TCGTATA',motifs.pgpd$Seq) | grepl('TCGTT[CT]AG',motifs.pgpd$Seq) | grepl('TCGT[AT]CC',motifs.pgpd$Seq) | grepl('TATACGA',motifs.pgpd$Seq) | grepl('TAAACGA',motifs.pgpd$Seq)
  motifs.pgpd$Upc2p_R = grepl('TATACGA',motifs.pgpd$Seq) | grepl('CT[AG]AACGA',motifs.pgpd$Seq) | grepl('GG[AT]ACGA',motifs.pgpd$Seq) | grepl('TCGTATA',motifs.pgpd$Seq) | grepl('TCGTTTA',motifs.pgpd$Seq)
  motifs.pgpd$Mbp1p_F = grepl('ACGCGT',motifs.pgpd$Seq)
  motifs.pgpd$Mbp1p_R = grepl('ACGCGT',motifs.pgpd$Seq)
  motifs.pgpd$Rim101p_F = grepl('TGCCAAG',motifs.pgpd$Seq) | grepl('CTTGGC[AG]',motifs.pgpd$Seq)
  motifs.pgpd$Rim101p_R = grepl('CTTGGCA',motifs.pgpd$Seq) | grepl('[CT]GCCAAG',motifs.pgpd$Seq)
  motifs.pgpd$Arg80p_F = grepl('CCTCTAAAGG',motifs.pgpd$Seq)
  motifs.pgpd$Arg80p_R = grepl('CCTTTAGAGG',motifs.pgpd$Seq)
  motifs.pgpd$Pho2p_F = grepl('[AT]TA[AT]T[AT]',motifs.pgpd$Seq)
  motifs.pgpd$Pho2p_R = grepl('[AT]A[AT]TA[AT]',motifs.pgpd$Seq)
  motifs.pgpd$Rds1p_F = grepl('[CGT]CGGCCG',motifs.pgpd$Seq)
  motifs.pgpd$Rds1p_R = grepl('CGGCCG[ACG]',motifs.pgpd$Seq)
  motifs.pgpd$Sum1p_F = grepl('G.C[AG]CAAA[AT]',motifs.pgpd$Seq)
  motifs.pgpd$Sum1p_R = grepl('[AT]TTTG[CT]G.C',motifs.pgpd$Seq)
  motifs.pgpd$Cbf1p_F = grepl('[AG]TCACGTG',motifs.pgpd$Seq) | grepl('[AG]TCAC[AG]TG',motifs.pgpd$Seq) | grepl('TCACGTT',motifs.pgpd$Seq) | grepl('CACGTG',motifs.pgpd$Seq)
  motifs.pgpd$Cbf1p_R = grepl('CACGTGA[CT]',motifs.pgpd$Seq) | grepl('CA[CT]GTGA[CT]',motifs.pgpd$Seq) | grepl('AACGTGA',motifs.pgpd$Seq) | grepl('CACGTG',motifs.pgpd$Seq)
  motifs.pgpd$Yap3p_F = grepl('TTACTAA',motifs.pgpd$Seq) | grepl('TGACTCA',motifs.pgpd$Seq) | grepl('ATTACGTAAT',motifs.pgpd$Seq)
  motifs.pgpd$Yap3p_R = grepl('TTAGTAA',motifs.pgpd$Seq) | grepl('TGAGTCA',motifs.pgpd$Seq) | grepl('ATTACGTAAT',motifs.pgpd$Seq)
  motifs.pgpd$Sfl1p_F = grepl('AGAA.TTCC',motifs.pgpd$Seq)
  motifs.pgpd$Sfl1p_R = grepl('GGAA.TTCT',motifs.pgpd$Seq)
  motifs.pgpd$Azf1p_F = grepl('AAGAAAAA',motifs.pgpd$Seq) | grepl('AAAAGAAA',motifs.pgpd$Seq)
  motifs.pgpd$Azf1p_R = grepl('TTTTTCTT',motifs.pgpd$Seq) | grepl('TTTCTTTT',motifs.pgpd$Seq)
  motifs.pgpd$Hap1p_F = grepl('CGG...TA.CGG',motifs.pgpd$Seq) | grepl('CGG......CGG',motifs.pgpd$Seq) | grepl('CGG...TA.CGG...TA',motifs.pgpd$Seq) | grepl('TCGCGTTCATATATATCGGATATACTCTTGAGGCACAATT',motifs.pgpd$Seq) | grepl('CGCATTTATCGC',motifs.pgpd$Seq) | grepl('CCG.TAT.TCC',motifs.pgpd$Seq) | grepl('CCGATA',motifs.pgpd$Seq)
  motifs.pgpd$Hap1p_R = grepl('CCG.TA...CCG',motifs.pgpd$Seq) | grepl('CCG......CCG',motifs.pgpd$Seq) | grepl('TA...CCG.TA...CCG',motifs.pgpd$Seq) | grepl('AATTGTGCCTCAAGAGTATATCCGATATATATGAACGCGA',motifs.pgpd$Seq) | grepl('GCGATAAATGCG',motifs.pgpd$Seq) | grepl('GGA.ATA.CGG',motifs.pgpd$Seq) | grepl('TATCGG',motifs.pgpd$Seq)
  motifs.pgpd$Tea1p_F = grepl('CGG..........CCG',motifs.pgpd$Seq)
  motifs.pgpd$Tea1p_R = grepl('CGG..........CCG',motifs.pgpd$Seq)
  motifs.pgpd$Cin5p_F = grepl('TTACTAA',motifs.pgpd$Seq)
  motifs.pgpd$Cin5p_R = grepl('TTAGTAA',motifs.pgpd$Seq)
  motifs.pgpd$Ino4p_F = grepl('ATGTGAAAT',motifs.pgpd$Seq) | grepl('CATGTGAAAT',motifs.pgpd$Seq) | grepl('[AT][CT]TTCA[CT][AG]TG[GC]',motifs.pgpd$Seq) | grepl('T[CT]TTCACATG[CT]',motifs.pgpd$Seq) | grepl('CACATGC',motifs.pgpd$Seq)
  motifs.pgpd$Ino4p_R = grepl('ATTTCACAT',motifs.pgpd$Seq) | grepl('ATTTCACATG',motifs.pgpd$Seq) | grepl('[GC]CA[CT][AG]TGAA[AG][AT]',motifs.pgpd$Seq) | grepl('[AG]CATGTGAA[AG]A',motifs.pgpd$Seq) | grepl('GCATGTG',motifs.pgpd$Seq)
  motifs.pgpd$Yap1p_F = grepl('TTACTAA',motifs.pgpd$Seq) | grepl('TGACTCA',motifs.pgpd$Seq) | grepl('TGACTAA',motifs.pgpd$Seq) | grepl('T[GT]ACAAA',motifs.pgpd$Seq) | grepl('TTACTCA',motifs.pgpd$Seq) | grepl('TGACAA',motifs.pgpd$Seq)
  motifs.pgpd$Yap1p_R = grepl('TTAGTAA',motifs.pgpd$Seq) | grepl('TGAGTCA',motifs.pgpd$Seq) | grepl('TTAGTCA',motifs.pgpd$Seq) | grepl('TTTGT[AC]A',motifs.pgpd$Seq) | grepl('TGAGTAA',motifs.pgpd$Seq) | grepl('TTGTCA',motifs.pgpd$Seq)
  motifs.pgpd$Put3p_F = grepl('CGG..........CCG',motifs.pgpd$Seq)
  motifs.pgpd$Put3p_R = grepl('CGG..........CCG',motifs.pgpd$Seq)
  motifs.pgpd$Gcr1p_F = grepl('CTTCC',motifs.pgpd$Seq) | grepl('C[AT]TCC',motifs.pgpd$Seq)
  motifs.pgpd$Gcr1p_R = grepl('GGAAG',motifs.pgpd$Seq) | grepl('GGA[AT]G',motifs.pgpd$Seq)
  motifs.pgpd$Cat8p_F = grepl('.CC[AGT]T[CT].[ACG].CCG.',motifs.pgpd$Seq) | grepl('[CT]CC.[CT]T.[AG][AG]CCG.',motifs.pgpd$Seq)
  motifs.pgpd$Cat8p_R = grepl('.CGG.[CGT].[AG]A[ACT]GG.',motifs.pgpd$Seq) | grepl('.CGG[CT][CT].A[AG].GG[AG]',motifs.pgpd$Seq)
  motifs.pgpd$Rap1p_F = grepl('ACACCC[AG][CT]ACA[CT]',motifs.pgpd$Seq) | grepl('ACACCCA[CT]ACA[CT][CT][CT]',motifs.pgpd$Seq)
  motifs.pgpd$Rap1p_R = grepl('[AG]TGT[AG][CT]GGGTGT',motifs.pgpd$Seq) | grepl('[AG][AG][AG]TGT[AG]TGGGTGT',motifs.pgpd$Seq)
  motifs.pgpd$Fkh2p_F = grepl('GT[AC]AACAA',motifs.pgpd$Seq) | grepl('[AG]TAAA[CT]AA',motifs.pgpd$Seq) | grepl('[AG][CT][AC]AA[CT]A',motifs.pgpd$Seq) | grepl('[AG][CT]AAACA[AT][AT]',motifs.pgpd$Seq) | grepl('AAACA',motifs.pgpd$Seq)
  motifs.pgpd$Fkh2p_R = grepl('TTGTT[GT]AC',motifs.pgpd$Seq) | grepl('TT[AG]TTTA[CT]',motifs.pgpd$Seq) | grepl('T[AG]TT[GT][AG][CT]',motifs.pgpd$Seq) | grepl('[AT][AT]TGTTT[AG][CT]',motifs.pgpd$Seq) | grepl('TGTTT',motifs.pgpd$Seq)
  motifs.pgpd$Ime1p_F = grepl('TTTTC[ACT][ACT]CG',motifs.pgpd$Seq)
  motifs.pgpd$Ime1p_R = grepl('CG[AGT][AGT]GAAAA',motifs.pgpd$Seq)
  motifs.pgpd$Tec1p_F = grepl('TTCTCACATTCTTC',motifs.pgpd$Seq) | grepl('CATTCT',motifs.pgpd$Seq) | grepl('CATTCC',motifs.pgpd$Seq) | grepl('CATTCTT',motifs.pgpd$Seq) | grepl('[AG][AC]ATTC[CT][CT]',motifs.pgpd$Seq) | grepl('ACATTCTTACATTCTT',motifs.pgpd$Seq) | grepl('GAATGT',motifs.pgpd$Seq) | grepl('[AG]GAATGT',motifs.pgpd$Seq)
  motifs.pgpd$Tec1p_R = grepl('GAAGAATGTGAGAA',motifs.pgpd$Seq) | grepl('AGAATG',motifs.pgpd$Seq) | grepl('GGAATG',motifs.pgpd$Seq) | grepl('AAGAATG',motifs.pgpd$Seq) | grepl('[AG][AG]GAAT[GT][CT]',motifs.pgpd$Seq) | grepl('AAGAATGTAAGAATGT',motifs.pgpd$Seq) | grepl('ACATTC',motifs.pgpd$Seq) | grepl('ACATTC[CT]',motifs.pgpd$Seq)
  motifs.pgpd$Met32p_F = grepl('AAACTGTGG',motifs.pgpd$Seq) | grepl('AAACTGTG',motifs.pgpd$Seq) | grepl('CTGTGGC',motifs.pgpd$Seq)
  motifs.pgpd$Met32p_R = grepl('CCACAGTTT',motifs.pgpd$Seq) | grepl('CACAGTTT',motifs.pgpd$Seq) | grepl('GCCACAG',motifs.pgpd$Seq)
  motifs.pgpd$Msn2p_F = grepl('CCCCT',motifs.pgpd$Seq) | grepl('[AG]GGGG',motifs.pgpd$Seq) | grepl('AGGGG',motifs.pgpd$Seq) | grepl('TCAGGGGT',motifs.pgpd$Seq) | grepl('AACGGGGT',motifs.pgpd$Seq)
  motifs.pgpd$Msn2p_R = grepl('AGGGG',motifs.pgpd$Seq) | grepl('CCCC[CT]',motifs.pgpd$Seq) | grepl('CCCCT',motifs.pgpd$Seq) | grepl('ACCCCTGA',motifs.pgpd$Seq) | grepl('ACCCCGTT',motifs.pgpd$Seq)
  motifs.pgpd$Rof1p_F = grepl('TTAAGGTT',motifs.pgpd$Seq)
  motifs.pgpd$Rof1p_R = grepl('AACCTTAA',motifs.pgpd$Seq)
  motifs.pgpd$Ppr1p_F = grepl('CGG......CCG',motifs.pgpd$Seq)
  motifs.pgpd$Ppr1p_R = grepl('CGG......CCG',motifs.pgpd$Seq)
  motifs.pgpd$Pdr3p_F = grepl('TCCGCGGA',motifs.pgpd$Seq) | grepl('TCCGTGGA',motifs.pgpd$Seq) | grepl('TCCACGGA',motifs.pgpd$Seq) | grepl('TCCGCGCA',motifs.pgpd$Seq) | grepl('TCCGCGGG',motifs.pgpd$Seq) | grepl('TTCCGCGGA[AT]',motifs.pgpd$Seq)
  motifs.pgpd$Pdr3p_R = grepl('TCCGCGGA',motifs.pgpd$Seq) | grepl('TCCACGGA',motifs.pgpd$Seq) | grepl('TCCGTGGA',motifs.pgpd$Seq) | grepl('TGCGCGGA',motifs.pgpd$Seq) | grepl('CCCGCGGA',motifs.pgpd$Seq) | grepl('[AT]TCCGCGGAA',motifs.pgpd$Seq)
  motifs.pgpd$Rlm1p_F = grepl('CTA[AT][AT][AT][AT]TAG',motifs.pgpd$Seq) | grepl('TA[AT][AT][AT][AT]TAG[AC]',motifs.pgpd$Seq)
  motifs.pgpd$Rlm1p_R = grepl('CTA[AT][AT][AT][AT]TAG',motifs.pgpd$Seq) | grepl('[GT]CTA[AT][AT][AT][AT]TA',motifs.pgpd$Seq)
  motifs.pgpd$Arr1p_F = grepl('TGATTACTAATCA',motifs.pgpd$Seq) | grepl('TGATTAATAATCA',motifs.pgpd$Seq)
  motifs.pgpd$Arr1p_R = grepl('TGATTAGTAATCA',motifs.pgpd$Seq) | grepl('TGATTATTAATCA',motifs.pgpd$Seq)
  motifs.pgpd$Stp2p_F = grepl('CGGCTC',motifs.pgpd$Seq) | grepl('CGGGGTG.......CGCACCG',motifs.pgpd$Seq) | grepl('C.CACC.G',motifs.pgpd$Seq)
  motifs.pgpd$Stp2p_R = grepl('GAGCCG',motifs.pgpd$Seq) | grepl('CGGTGCG.......CACCCCG',motifs.pgpd$Seq) | grepl('C.GGTG.G',motifs.pgpd$Seq)
  motifs.pgpd$Haa1p_F = grepl('GGCGAGGGG',motifs.pgpd$Seq) | grepl('GGCGAGAGG',motifs.pgpd$Seq) | grepl('GGCGCGGGG',motifs.pgpd$Seq) | grepl('GGCCAGGGG',motifs.pgpd$Seq) | grepl('AGCGAGGGG',motifs.pgpd$Seq) | grepl('[GC][AC]GG[GC]G',motifs.pgpd$Seq) | grepl('GAGGCG',motifs.pgpd$Seq)
  motifs.pgpd$Haa1p_R = grepl('CCCCTCGCC',motifs.pgpd$Seq) | grepl('CCTCTCGCC',motifs.pgpd$Seq) | grepl('CCCCGCGCC',motifs.pgpd$Seq) | grepl('CCCCTGGCC',motifs.pgpd$Seq) | grepl('CCCCTCGCT',motifs.pgpd$Seq) | grepl('C[GC]CC[GT][GC]',motifs.pgpd$Seq) | grepl('CGCCTC',motifs.pgpd$Seq)
  motifs.pgpd$Mot3p_F = grepl('[AT]AGGTA',motifs.pgpd$Seq) | grepl('CAGG[CT]A',motifs.pgpd$Seq) | grepl('AAGAGG',motifs.pgpd$Seq) | grepl('ATGGAT',motifs.pgpd$Seq) | grepl('TAGGTA',motifs.pgpd$Seq) | grepl('TAGGAT',motifs.pgpd$Seq) | grepl('T[AC]GGAA',motifs.pgpd$Seq) | grepl('AAGG[GT]A',motifs.pgpd$Seq) | grepl('AAGG[AT]T',motifs.pgpd$Seq)
  motifs.pgpd$Mot3p_R = grepl('TACCT[AT]',motifs.pgpd$Seq) | grepl('T[AG]CCTG',motifs.pgpd$Seq) | grepl('CCTCTT',motifs.pgpd$Seq) | grepl('ATCCAT',motifs.pgpd$Seq) | grepl('TACCTA',motifs.pgpd$Seq) | grepl('ATCCTA',motifs.pgpd$Seq) | grepl('TTCC[GT]A',motifs.pgpd$Seq) | grepl('T[AC]CCTT',motifs.pgpd$Seq) | grepl('A[AT]CCTT',motifs.pgpd$Seq)
  motifs.pgpd$Rgt1p_F = grepl('CGGA..A',motifs.pgpd$Seq)
  motifs.pgpd$Rgt1p_R = grepl('T..TCCG',motifs.pgpd$Seq)
  motifs.pgpd$Com2p_F = grepl('AGGGG',motifs.pgpd$Seq) | grepl('TCAGGGGT',motifs.pgpd$Seq) | grepl('AATAGGAG',motifs.pgpd$Seq) | grepl('ATAGGGG',motifs.pgpd$Seq) | grepl('ATAGGAG',motifs.pgpd$Seq) | grepl('ATAGGGT',motifs.pgpd$Seq)
  motifs.pgpd$Com2p_R = grepl('CCCCT',motifs.pgpd$Seq) | grepl('ACCCCTGA',motifs.pgpd$Seq) | grepl('CTCCTATT',motifs.pgpd$Seq) | grepl('CCCCTAT',motifs.pgpd$Seq) | grepl('CTCCTAT',motifs.pgpd$Seq) | grepl('ACCCTAT',motifs.pgpd$Seq)
  motifs.pgpd$Skn7p_F = grepl('GGC[CT]GGC',motifs.pgpd$Seq) | grepl('AGAACGTTC',motifs.pgpd$Seq) | grepl('ATTTGG[CT]TGGGCC',motifs.pgpd$Seq) | grepl('GGCGAGATCT',motifs.pgpd$Seq) | grepl('GGCCCAGA',motifs.pgpd$Seq) | grepl('GGCCGGC',motifs.pgpd$Seq) | grepl('GGCCAGA',motifs.pgpd$Seq) | grepl('GGCTGGC',motifs.pgpd$Seq) | grepl('GGCC[AG]',motifs.pgpd$Seq)
  motifs.pgpd$Skn7p_R = grepl('GCC[AG]GCC',motifs.pgpd$Seq) | grepl('GAACGTTCT',motifs.pgpd$Seq) | grepl('GGCCCA[AG]CCAAAT',motifs.pgpd$Seq) | grepl('AGATCTCGCC',motifs.pgpd$Seq) | grepl('TCTGGGCC',motifs.pgpd$Seq) | grepl('GCCGGCC',motifs.pgpd$Seq) | grepl('TCTGGCC',motifs.pgpd$Seq) | grepl('GCCAGCC',motifs.pgpd$Seq) | grepl('[CT]GGCC',motifs.pgpd$Seq)
  motifs.pgpd$Gal4p_F = grepl('CGG...........CCG',motifs.pgpd$Seq)
  motifs.pgpd$Gal4p_R = grepl('CGG...........CCG',motifs.pgpd$Seq)
  motifs.pgpd$Ino2p_F = grepl('ATGTGAAAT',motifs.pgpd$Seq) | grepl('CATGTGAAAT',motifs.pgpd$Seq) | grepl('[AT][CT]TTCA[CT][AG]TG[GC]',motifs.pgpd$Seq) | grepl('CACATGC',motifs.pgpd$Seq)
  motifs.pgpd$Ino2p_R = grepl('ATTTCACAT',motifs.pgpd$Seq) | grepl('ATTTCACATG',motifs.pgpd$Seq) | grepl('[GC]CA[CT][AG]TGAA[AG][AT]',motifs.pgpd$Seq) | grepl('GCATGTG',motifs.pgpd$Seq)
  motifs.pgpd$Ash1p_F = grepl('[CT]TGAT',motifs.pgpd$Seq)
  motifs.pgpd$Ash1p_R = grepl('ATCA[AG]',motifs.pgpd$Seq)
  motifs.pgpd$Tda9p_F = grepl('GAGGGG',motifs.pgpd$Seq)
  motifs.pgpd$Tda9p_R = grepl('CCCCTC',motifs.pgpd$Seq)
  motifs.pgpd$Fzf1p_F = grepl('CGTATCGTATAAGGCAACAATAG',motifs.pgpd$Seq) | grepl('[CT]G[GC][AC].[AC]CTATCA[CT]TT[CT][CT]',motifs.pgpd$Seq)
  motifs.pgpd$Fzf1p_R = grepl('CTATTGTTGCCTTATACGATACG',motifs.pgpd$Seq) | grepl('[AG][AG]AA[AG]TGATAG[GT].[GT][GC]C[AG]',motifs.pgpd$Seq)
  motifs.pgpd$Gsm1p_F = grepl('CGG........CGG',motifs.pgpd$Seq) | grepl('CGG.........CGG',motifs.pgpd$Seq)
  motifs.pgpd$Gsm1p_R = grepl('CCG........CCG',motifs.pgpd$Seq) | grepl('CCG.........CCG',motifs.pgpd$Seq)
  motifs.pgpd$Aft1p_F = grepl('[CT][AG]CACCC[AG]',motifs.pgpd$Seq) | grepl('TGCACCC',motifs.pgpd$Seq) | grepl('GGCACCC',motifs.pgpd$Seq) | grepl('TGCACCCA',motifs.pgpd$Seq)
  motifs.pgpd$Aft1p_R = grepl('[CT]GGGTG[CT][AG]',motifs.pgpd$Seq) | grepl('GGGTGCA',motifs.pgpd$Seq) | grepl('GGGTGCC',motifs.pgpd$Seq) | grepl('TGGGTGCA',motifs.pgpd$Seq)
  motifs.pgpd$Ecm22p_F = grepl('TCGTATA',motifs.pgpd$Seq)
  motifs.pgpd$Ecm22p_R = grepl('TATACGA',motifs.pgpd$Seq)
  motifs.pgpd$Xbp1p_F = grepl('CTCGA',motifs.pgpd$Seq)
  motifs.pgpd$Xbp1p_R = grepl('TCGAG',motifs.pgpd$Seq)
  motifs.pgpd$Reb1p_F = grepl('CCGGGTAA',motifs.pgpd$Seq) | grepl('TGTTACCCGT',motifs.pgpd$Seq) | grepl('CCGGGTGGAT',motifs.pgpd$Seq)
  motifs.pgpd$Reb1p_R = grepl('TTACCCGG',motifs.pgpd$Seq) | grepl('ACGGGTAACA',motifs.pgpd$Seq) | grepl('ATCCACCCGG',motifs.pgpd$Seq)
  motifs.pgpd$Hcm1p_F = grepl('[AT]AA[CT]AAACAA[AT]',motifs.pgpd$Seq)
  motifs.pgpd$Hcm1p_R = grepl('[AT]TTGTTT[AG]TT[AT]',motifs.pgpd$Seq)
  motifs.pgpd$Rpn4p_F = grepl('GGTGGCAAA',motifs.pgpd$Seq)
  motifs.pgpd$Rpn4p_R = grepl('TTTGCCACC',motifs.pgpd$Seq)
  motifs.pgpd$Hap2p_F = grepl('T.ATTGGT',motifs.pgpd$Seq)
  motifs.pgpd$Hap2p_R = grepl('ACCAAT.A',motifs.pgpd$Seq)
  motifs.pgpd$Mal63p_F = grepl('[AC]GC.........[AC]G[GC]',motifs.pgpd$Seq) | grepl('CGG.........CGG',motifs.pgpd$Seq) | grepl('CGC.........CGC',motifs.pgpd$Seq) | grepl('CGG.........CGC',motifs.pgpd$Seq)
  motifs.pgpd$Mal63p_R = grepl('[GC]C[GT].........GC[GT]',motifs.pgpd$Seq) | grepl('CCG.........CCG',motifs.pgpd$Seq) | grepl('GCG.........GCG',motifs.pgpd$Seq) | grepl('GCG.........CCG',motifs.pgpd$Seq)
  motifs.pgpd$Flo8p_F = grepl('CGGGGTTTTCT',motifs.pgpd$Seq)
  motifs.pgpd$Flo8p_R = grepl('AGAAAACCCCG',motifs.pgpd$Seq)
  motifs.pgpd$Tye7p_F = grepl('TCACGTG',motifs.pgpd$Seq)
  motifs.pgpd$Tye7p_R = grepl('CACGTGA',motifs.pgpd$Seq)
  motifs.pgpd$Swi4p_F = grepl('CACGAAA',motifs.pgpd$Seq) | grepl('TTTTCGCT',motifs.pgpd$Seq)
  motifs.pgpd$Swi4p_R = grepl('TTTCGTG',motifs.pgpd$Seq) | grepl('AGCGAAAA',motifs.pgpd$Seq)
  motifs.pgpd$Gln3p_F = grepl('GATAAG',motifs.pgpd$Seq) | grepl('GATTAG',motifs.pgpd$Seq)
  motifs.pgpd$Gln3p_R = grepl('CTTATC',motifs.pgpd$Seq) | grepl('CTAATC',motifs.pgpd$Seq)
  motifs.pgpd$Adr1p_F = grepl('TTGG[AG]G',motifs.pgpd$Seq)
  motifs.pgpd$Adr1p_R = grepl('C[CT]CCAA',motifs.pgpd$Seq)
  motifs.pgpd$Fkh1p_F = grepl('GT[AC]AACAA',motifs.pgpd$Seq) | grepl('[AG]TAAA[CT]AA',motifs.pgpd$Seq) | grepl('[AG][CT][AC]AA[CT]A',motifs.pgpd$Seq) | grepl('[AG][CT]AAACA[AT][AT]',motifs.pgpd$Seq) | grepl('AAACA',motifs.pgpd$Seq)
  motifs.pgpd$Fkh1p_R = grepl('TTGTT[GT]AC',motifs.pgpd$Seq) | grepl('TT[AG]TTTA[CT]',motifs.pgpd$Seq) | grepl('T[AG]TT[GT][AG][CT]',motifs.pgpd$Seq) | grepl('[AT][AT]TGTTT[AG][CT]',motifs.pgpd$Seq) | grepl('TGTTT',motifs.pgpd$Seq)
  motifs.pgpd$Mig2p_F = grepl('AATT[AG]TCCGGGG',motifs.pgpd$Seq) | grepl('TATAGTGCGGGG',motifs.pgpd$Seq) | grepl('ATTTTTGTGGGG',motifs.pgpd$Seq) | grepl('GAATA[CT]CTGGGG',motifs.pgpd$Seq) | grepl('TTATTTCTGGGG',motifs.pgpd$Seq) | grepl('TAAAAGCCGGGG',motifs.pgpd$Seq) | grepl('TTAAAAGCGGGG',motifs.pgpd$Seq) | grepl('TTAAACGTGGGG',motifs.pgpd$Seq) | grepl('ATAAAACTGGGG',motifs.pgpd$Seq) | grepl('ATAAACGCGGGG',motifs.pgpd$Seq) | grepl('ATTAATGTGGGG',motifs.pgpd$Seq) | grepl('AAAAATGCGGGG',motifs.pgpd$Seq) | grepl('ATTTTGCGGGG',motifs.pgpd$Seq) | grepl('CCCC[AG]C',motifs.pgpd$Seq)
  motifs.pgpd$Mig2p_R = grepl('CCCCGGA[CT]AATT',motifs.pgpd$Seq) | grepl('CCCCGCACTATA',motifs.pgpd$Seq) | grepl('CCCCACAAAAAT',motifs.pgpd$Seq) | grepl('CCCCAG[AG]TATTC',motifs.pgpd$Seq) | grepl('CCCCAGAAATAA',motifs.pgpd$Seq) | grepl('CCCCGGCTTTTA',motifs.pgpd$Seq) | grepl('CCCCGCTTTTAA',motifs.pgpd$Seq) | grepl('CCCCACGTTTAA',motifs.pgpd$Seq) | grepl('CCCCAGTTTTAT',motifs.pgpd$Seq) | grepl('CCCCGCGTTTAT',motifs.pgpd$Seq) | grepl('CCCCACATTAAT',motifs.pgpd$Seq) | grepl('CCCCGCATTTTT',motifs.pgpd$Seq) | grepl('CCCCGCAAAAT',motifs.pgpd$Seq) | grepl('G[CT]GGGG',motifs.pgpd$Seq)
  motifs.pgpd$Ixr1p_F = grepl('[GT]TT[GC]AA[CT][GT]GTT[CT]A[GC]A',motifs.pgpd$Seq)
  motifs.pgpd$Ixr1p_R = grepl('T[GC]T[AG]AAC[AC][AG]TT[GC]AA[AC]',motifs.pgpd$Seq)
  motifs.pgpd$Stb5p_F = grepl('TCTCCGCGAAC',motifs.pgpd$Seq) | grepl('CGG.[GC]',motifs.pgpd$Seq)
  motifs.pgpd$Stb5p_R = grepl('GTTCGCGGAGA',motifs.pgpd$Seq) | grepl('[GC].CCG',motifs.pgpd$Seq)
  motifs.pgpd$Met4p_F = grepl('TCACGTG',motifs.pgpd$Seq)
  motifs.pgpd$Met4p_R = grepl('CACGTGA',motifs.pgpd$Seq)
  motifs.pgpd$Abf1p_F = grepl('AGCCGTAAATAGTTATCTTCCAAG',motifs.pgpd$Seq) | grepl('[AG]TC[AG][CT][CT][CT]...ACG',motifs.pgpd$Seq) | grepl('[AG]TC[AG][CT][CGT]....ACG',motifs.pgpd$Seq) | grepl('TC.......ACG',motifs.pgpd$Seq) | grepl('[AG]TC[AG]......ACG.[AG]',motifs.pgpd$Seq) | grepl('TC[AG]T.....A[CT]GA',motifs.pgpd$Seq) | grepl('T..CGT......TGAT',motifs.pgpd$Seq)
  motifs.pgpd$Abf1p_R = grepl('CTTGGAAGATAACTATTTACGGCT',motifs.pgpd$Seq) | grepl('CGT...[AG][AG][AG][CT]GA[CT]',motifs.pgpd$Seq) | grepl('CGT....[ACG][AG][CT]GA[CT]',motifs.pgpd$Seq) | grepl('CGT.......GA',motifs.pgpd$Seq) | grepl('[CT].CGT......[CT]GA[CT]',motifs.pgpd$Seq) | grepl('TC[AG]T.....A[CT]GA',motifs.pgpd$Seq) | grepl('ATCA......ACG..A',motifs.pgpd$Seq)
  motifs.pgpd$Cad1p_F = grepl('TTACTAA',motifs.pgpd$Seq)
  motifs.pgpd$Cad1p_R = grepl('TTAGTAA',motifs.pgpd$Seq)
  motifs.pgpd$Mig3p_F = grepl('AATT[AG]TCCGGGG',motifs.pgpd$Seq) | grepl('TATAGTGCGGGG',motifs.pgpd$Seq) | grepl('ATTTTTGTGGGG',motifs.pgpd$Seq) | grepl('GAATA[CT]CTGGGG',motifs.pgpd$Seq) | grepl('TTATTTCTGGGG',motifs.pgpd$Seq) | grepl('TAAAAGCCGGGG',motifs.pgpd$Seq) | grepl('TTAAAAGCGGGG',motifs.pgpd$Seq) | grepl('TTAAACGTGGGG',motifs.pgpd$Seq) | grepl('ATAAAACTGGGG',motifs.pgpd$Seq) | grepl('ATAAACGCGGGG',motifs.pgpd$Seq) | grepl('ATTAATGTGGGG',motifs.pgpd$Seq) | grepl('AAAAATGCGGGG',motifs.pgpd$Seq) | grepl('ATTTTGCGGGG',motifs.pgpd$Seq)
  motifs.pgpd$Mig3p_R = grepl('CCCCGGA[CT]AATT',motifs.pgpd$Seq) | grepl('CCCCGCACTATA',motifs.pgpd$Seq) | grepl('CCCCACAAAAAT',motifs.pgpd$Seq) | grepl('CCCCAG[AG]TATTC',motifs.pgpd$Seq) | grepl('CCCCAGAAATAA',motifs.pgpd$Seq) | grepl('CCCCGGCTTTTA',motifs.pgpd$Seq) | grepl('CCCCGCTTTTAA',motifs.pgpd$Seq) | grepl('CCCCACGTTTAA',motifs.pgpd$Seq) | grepl('CCCCAGTTTTAT',motifs.pgpd$Seq) | grepl('CCCCGCGTTTAT',motifs.pgpd$Seq) | grepl('CCCCACATTAAT',motifs.pgpd$Seq) | grepl('CCCCGCATTTTT',motifs.pgpd$Seq) | grepl('CCCCGCAAAAT',motifs.pgpd$Seq)
  motifs.pgpd$Rtg1p_F = grepl('GTCAC',motifs.pgpd$Seq) | grepl('GGTAC',motifs.pgpd$Seq)
  motifs.pgpd$Rtg1p_R = grepl('GTGAC',motifs.pgpd$Seq) | grepl('GTACC',motifs.pgpd$Seq)
  motifs.pgpd$Zap1p_F = grepl('ACC[CT][CT].AAGGT',motifs.pgpd$Seq) | grepl('ACCGTCACTGC',motifs.pgpd$Seq) | grepl('TTCTTTATGGT',motifs.pgpd$Seq) | grepl('ATAACC',motifs.pgpd$Seq) | grepl('ACCTT.AAGGT',motifs.pgpd$Seq) | grepl('ACCTTGAAGGT',motifs.pgpd$Seq)
  motifs.pgpd$Zap1p_R = grepl('ACCTT.[AG][AG]GGT',motifs.pgpd$Seq) | grepl('GCAGTGACGGT',motifs.pgpd$Seq) | grepl('ACCATAAAGAA',motifs.pgpd$Seq) | grepl('GGTTAT',motifs.pgpd$Seq) | grepl('ACCTT.AAGGT',motifs.pgpd$Seq) | grepl('ACCTTCAAGGT',motifs.pgpd$Seq)
  motifs.pgpd$Pho4p_F = grepl('CACGT[GT]',motifs.pgpd$Seq) | grepl('CACGTGGG',motifs.pgpd$Seq)
  motifs.pgpd$Pho4p_R = grepl('[AC]ACGTG',motifs.pgpd$Seq) | grepl('CCCACGTG',motifs.pgpd$Seq)
  motifs.pgpd$Rme1p_F = grepl('GTACCACAAAA',motifs.pgpd$Seq) | grepl('GTACCTCAAAA',motifs.pgpd$Seq) | grepl('GAACCTCAA[AG]A',motifs.pgpd$Seq)
  motifs.pgpd$Rme1p_R = grepl('TTTTGTGGTAC',motifs.pgpd$Seq) | grepl('TTTTGAGGTAC',motifs.pgpd$Seq) | grepl('T[CT]TTGAGGTTC',motifs.pgpd$Seq)
  motifs.pgpd$Uga3p_F = grepl('[GC]GCGG.[AT]TTT',motifs.pgpd$Seq)
  motifs.pgpd$Uga3p_R = grepl('AAA[AT].CCGC[GC]',motifs.pgpd$Seq)
  motifs.pgpd$Gat1p_F = grepl('GATAAG',motifs.pgpd$Seq)
  motifs.pgpd$Gat1p_R = grepl('CTTATC',motifs.pgpd$Seq)
  motifs.pgpd$Mot2p_F = grepl('ATATA',motifs.pgpd$Seq)
  motifs.pgpd$Mot2p_R = grepl('TATAT',motifs.pgpd$Seq)
  motifs.pgpd$Nrg1p_F = grepl('CCCCT',motifs.pgpd$Seq) | grepl('CCCTC',motifs.pgpd$Seq) | grepl('GGACCCT',motifs.pgpd$Seq)
  motifs.pgpd$Nrg1p_R = grepl('AGGGG',motifs.pgpd$Seq) | grepl('GAGGG',motifs.pgpd$Seq) | grepl('AGGGTCC',motifs.pgpd$Seq)
  motifs.pgpd$Smp1p_F = grepl('ACTACTA[AT][AT][AT][AT]TAG',motifs.pgpd$Seq)
  motifs.pgpd$Smp1p_R = grepl('CTA[AT][AT][AT][AT]TAGTAGT',motifs.pgpd$Seq)
  motifs.pgpd$Hsf1p_F = grepl('.GAA..TTC.',motifs.pgpd$Seq) | grepl('.TTC..GAA.',motifs.pgpd$Seq) | grepl('.GAA..TTC..GAA.',motifs.pgpd$Seq) | grepl('.GAA.......GAA.......GAA.',motifs.pgpd$Seq) | grepl('CTTCATGAAA',motifs.pgpd$Seq) | grepl('TTCTAGAA..TTCTAGAA',motifs.pgpd$Seq) | grepl('TTC..GAA',motifs.pgpd$Seq) | grepl('TTC..GAA..TTC',motifs.pgpd$Seq) | grepl('GAA..TTC..GAA',motifs.pgpd$Seq) | grepl('TTC..GAA..TTC..GAA',motifs.pgpd$Seq) | grepl('TTC..GAA.......GAA',motifs.pgpd$Seq) | grepl('TTC.......TTC..GAA',motifs.pgpd$Seq) | grepl('TTC.......TTC.......TTC',motifs.pgpd$Seq) | grepl('GAA.......GAA.......GAA',motifs.pgpd$Seq) | grepl('TTCTTC.......TTC',motifs.pgpd$Seq)
  motifs.pgpd$Hsf1p_R = grepl('.GAA..TTC.',motifs.pgpd$Seq) | grepl('.TTC..GAA.',motifs.pgpd$Seq) | grepl('.TTC..GAA..TTC.',motifs.pgpd$Seq) | grepl('.TTC.......TTC.......TTC.',motifs.pgpd$Seq) | grepl('TTTCATGAAG',motifs.pgpd$Seq) | grepl('TTCTAGAA..TTCTAGAA',motifs.pgpd$Seq) | grepl('TTC..GAA',motifs.pgpd$Seq) | grepl('GAA..TTC..GAA',motifs.pgpd$Seq) | grepl('TTC..GAA..TTC',motifs.pgpd$Seq) | grepl('TTC..GAA..TTC..GAA',motifs.pgpd$Seq) | grepl('TTC.......TTC..GAA',motifs.pgpd$Seq) | grepl('TTC..GAA.......GAA',motifs.pgpd$Seq) | grepl('GAA.......GAA.......GAA',motifs.pgpd$Seq) | grepl('TTC.......TTC.......TTC',motifs.pgpd$Seq) | grepl('GAA.......GAAGAA',motifs.pgpd$Seq)
  motifs.pgpd$Bas1p_F = grepl('TGACTC',motifs.pgpd$Seq)
  motifs.pgpd$Bas1p_R = grepl('GAGTCA',motifs.pgpd$Seq)
  motifs.pgpd$Gcn4p_F = grepl('TTACGTAA',motifs.pgpd$Seq) | grepl('TGATTCA',motifs.pgpd$Seq) | grepl('TGACTGA',motifs.pgpd$Seq) | grepl('TGACT[AC]T',motifs.pgpd$Seq) | grepl('[AG][AG]TGACTC',motifs.pgpd$Seq) | grepl('TGACTC',motifs.pgpd$Seq) | grepl('TGA[GC]TCA',motifs.pgpd$Seq) | grepl('TTGCGCAA',motifs.pgpd$Seq) | grepl('TTGCGTGA',motifs.pgpd$Seq) | grepl('GCACGTAG',motifs.pgpd$Seq) | grepl('CACGTG',motifs.pgpd$Seq)
  motifs.pgpd$Gcn4p_R = grepl('TTACGTAA',motifs.pgpd$Seq) | grepl('TGAATCA',motifs.pgpd$Seq) | grepl('TCAGTCA',motifs.pgpd$Seq) | grepl('A[GT]AGTCA',motifs.pgpd$Seq) | grepl('GAGTCA[CT][CT]',motifs.pgpd$Seq) | grepl('GAGTCA',motifs.pgpd$Seq) | grepl('TGA[GC]TCA',motifs.pgpd$Seq) | grepl('TTGCGCAA',motifs.pgpd$Seq) | grepl('TCACGCAA',motifs.pgpd$Seq) | grepl('CTACGTGC',motifs.pgpd$Seq) | grepl('CACGTG',motifs.pgpd$Seq)
  motifs.pgpd$Arg81p_F = grepl('[AT][AT]...[CT]A[AG].[ACT]...[ACG][ACG]CG[AG]',motifs.pgpd$Seq) | grepl('AAGTACAGTTAATAACGA',motifs.pgpd$Seq) | grepl('AAGTACAGTTAATAACGG',motifs.pgpd$Seq) | grepl('AAGTGCAACTGACTGCGA',motifs.pgpd$Seq) | grepl('AATGGAAATGGATAGCGA',motifs.pgpd$Seq)
  motifs.pgpd$Arg81p_R = grepl('[CT]CG[CGT][CGT]...[AGT].[CT]T[AG]...[AT][AT]',motifs.pgpd$Seq) | grepl('TCGTTATTAACTGTACTT',motifs.pgpd$Seq) | grepl('CCGTTATTAACTGTACTT',motifs.pgpd$Seq) | grepl('TCGCAGTCAGTTGCACTT',motifs.pgpd$Seq) | grepl('TCGCTATCCATTTCCATT',motifs.pgpd$Seq)
  motifs.pgpd$Stb4p_F = grepl('[CT]TCGGAA',motifs.pgpd$Seq)
  motifs.pgpd$Stb4p_R = grepl('TTCCGA[AG]',motifs.pgpd$Seq)
  motifs.pgpd$Hot1p_F = grepl('GGGACAAA',motifs.pgpd$Seq) | grepl('CATTTGGC',motifs.pgpd$Seq) | grepl('CACTTTGAC',motifs.pgpd$Seq)
  motifs.pgpd$Hot1p_R = grepl('TTTGTCCC',motifs.pgpd$Seq) | grepl('GCCAAATG',motifs.pgpd$Seq) | grepl('GTCAAAGTG',motifs.pgpd$Seq)
  motifs.pgpd$Pip2p_F = grepl('CCG...TA',motifs.pgpd$Seq)
  motifs.pgpd$Pip2p_R = grepl('TA...CGG',motifs.pgpd$Seq)
  motifs.pgpd$Swi5p_F = grepl('ACCAGC',motifs.pgpd$Seq)
  motifs.pgpd$Swi5p_R = grepl('GCTGGT',motifs.pgpd$Seq)
  motifs.pgpd$Cha4p_F = grepl('CGG..........CCG',motifs.pgpd$Seq) | grepl('CGG..........CCA',motifs.pgpd$Seq)
  motifs.pgpd$Cha4p_R = grepl('CGG..........CCG',motifs.pgpd$Seq) | grepl('TGG..........CCG',motifs.pgpd$Seq)
  motifs.pgpd$Yap5p_F = grepl('TTACTAA',motifs.pgpd$Seq)
  motifs.pgpd$Yap5p_R = grepl('TTAGTAA',motifs.pgpd$Seq)
  motifs.pgpd$Yrr1p_F = grepl('[AT]CCG[CT][GT][GT][AT][AT]',motifs.pgpd$Seq)
  motifs.pgpd$Yrr1p_R = grepl('[AT][AT][AC][AC][AG]CGG[AT]',motifs.pgpd$Seq)
  motifs.pgpd$Ace2p_F = grepl('ACCAGC',motifs.pgpd$Seq)
  motifs.pgpd$Ace2p_R = grepl('GCTGGT',motifs.pgpd$Seq)
  motifs.pgpd$Pdr1p_F = grepl('TCCGCGGA',motifs.pgpd$Seq) | grepl('TCCGTGGA',motifs.pgpd$Seq) | grepl('TCCACGGA',motifs.pgpd$Seq) | grepl('TCCGCGCA',motifs.pgpd$Seq) | grepl('TCCGCGGG',motifs.pgpd$Seq)
  motifs.pgpd$Pdr1p_R = grepl('TCCGCGGA',motifs.pgpd$Seq) | grepl('TCCACGGA',motifs.pgpd$Seq) | grepl('TCCGTGGA',motifs.pgpd$Seq) | grepl('TGCGCGGA',motifs.pgpd$Seq) | grepl('CCCGCGGA',motifs.pgpd$Seq)
  motifs.pgpd$Rph1p_F = grepl('CCCCT',motifs.pgpd$Seq) | grepl('AGGGG',motifs.pgpd$Seq)
  motifs.pgpd$Rph1p_R = grepl('AGGGG',motifs.pgpd$Seq) | grepl('CCCCT',motifs.pgpd$Seq)
  motifs.pgpd$Met31p_F = grepl('AAACTGTGG',motifs.pgpd$Seq) | grepl('AAACTGTG',motifs.pgpd$Seq) | grepl('CTGTGGC',motifs.pgpd$Seq)
  motifs.pgpd$Met31p_R = grepl('CCACAGTTT',motifs.pgpd$Seq) | grepl('CACAGTTT',motifs.pgpd$Seq) | grepl('GCCACAG',motifs.pgpd$Seq)
  motifs.pgpd$Znf1p_F = grepl('CGG.........CGG',motifs.pgpd$Seq)
  motifs.pgpd$Znf1p_R = grepl('CCG.........CCG',motifs.pgpd$Seq)
  motifs.pgpd$Aft2p_F = grepl('[CT][AG]CACCC[AG]',motifs.pgpd$Seq) | grepl('TGCACCC',motifs.pgpd$Seq) | grepl('GGCACCC',motifs.pgpd$Seq) | grepl('CGCACCC',motifs.pgpd$Seq)
  motifs.pgpd$Aft2p_R = grepl('[CT]GGGTG[CT][AG]',motifs.pgpd$Seq) | grepl('GGGTGCA',motifs.pgpd$Seq) | grepl('GGGTGCC',motifs.pgpd$Seq) | grepl('GGGTGCG',motifs.pgpd$Seq)
  motifs.pgpd$Ste12p_F = grepl('TGAAACA',motifs.pgpd$Seq) | grepl('[AC]TGAAACA',motifs.pgpd$Seq) | grepl('ATGAAACA',motifs.pgpd$Seq) | grepl('GGAAACA',motifs.pgpd$Seq) | grepl('ATAAAACA',motifs.pgpd$Seq) | grepl('ATGCAACA',motifs.pgpd$Seq) | grepl('ATGAAACAATGAAACA',motifs.pgpd$Seq) | grepl('TTGAAACAATGAAACA',motifs.pgpd$Seq) | grepl('ATGAAACAATGAGACA',motifs.pgpd$Seq) | grepl('ATGAAACAATGAAACG',motifs.pgpd$Seq) | grepl('TTGAAACAATGAAACG',motifs.pgpd$Seq) | grepl('CGTTTCAAAATGAAACA',motifs.pgpd$Seq) | grepl('ATGAACA...ATGAAACA',motifs.pgpd$Seq) | grepl('TGTTTCA',motifs.pgpd$Seq) | grepl('GTTTCA...TGAAAC',motifs.pgpd$Seq)
  motifs.pgpd$Ste12p_R = grepl('TGTTTCA',motifs.pgpd$Seq) | grepl('TGTTTCA[GT]',motifs.pgpd$Seq) | grepl('TGTTTCAT',motifs.pgpd$Seq) | grepl('TGTTTCC',motifs.pgpd$Seq) | grepl('TGTTTTAT',motifs.pgpd$Seq) | grepl('TGTTGCAT',motifs.pgpd$Seq) | grepl('TGTTTCATTGTTTCAT',motifs.pgpd$Seq) | grepl('TGTTTCATTGTTTCAA',motifs.pgpd$Seq) | grepl('TGTCTCATTGTTTCAT',motifs.pgpd$Seq) | grepl('CGTTTCATTGTTTCAT',motifs.pgpd$Seq) | grepl('CGTTTCATTGTTTCAA',motifs.pgpd$Seq) | grepl('TGTTTCATTTTGAAACG',motifs.pgpd$Seq) | grepl('TGTTTCAT...TGTTCAT',motifs.pgpd$Seq) | grepl('TGAAACA',motifs.pgpd$Seq) | grepl('GTTTCA...TGAAAC',motifs.pgpd$Seq)
  motifs.pgpd$Cst6p_F = grepl('TGACGTCA',motifs.pgpd$Seq) | grepl('TTACGTAA',motifs.pgpd$Seq) | grepl('GTGACGT',motifs.pgpd$Seq)
  motifs.pgpd$Cst6p_R = grepl('TGACGTCA',motifs.pgpd$Seq) | grepl('TTACGTAA',motifs.pgpd$Seq) | grepl('ACGTCAC',motifs.pgpd$Seq)
  motifs.pgpd$Ume6p_F = grepl('TAGCCGCCGA',motifs.pgpd$Seq) | grepl('TAGGCGGC',motifs.pgpd$Seq) | grepl('TCGGCGGCT',motifs.pgpd$Seq) | grepl('TGGGTGGCTA',motifs.pgpd$Seq)
  motifs.pgpd$Ume6p_R = grepl('TCGGCGGCTA',motifs.pgpd$Seq) | grepl('GCCGCCTA',motifs.pgpd$Seq) | grepl('AGCCGCCGA',motifs.pgpd$Seq) | grepl('TAGCCACCCA',motifs.pgpd$Seq)
  motifs.pgpd$Hap3p_F = grepl('T.ATTGGT',motifs.pgpd$Seq)
  motifs.pgpd$Hap3p_R = grepl('ACCAAT.A',motifs.pgpd$Seq)
  
return(motifs.pgpd)
}

merge.fr = function(motifs.pgpd) {
  motifs.pgpd$Gis1p = motifs.pgpd$Gis1p_F | motifs.pgpd$Gis1p_R; motifs.pgpd$Gis1p_F=NULL; motifs.pgpd$Gis1p_R=NULL
  motifs.pgpd$Msn4p = motifs.pgpd$Msn4p_F | motifs.pgpd$Msn4p_R; motifs.pgpd$Msn4p_F=NULL; motifs.pgpd$Msn4p_R=NULL
  motifs.pgpd$Rds1p = motifs.pgpd$Rds1p_F | motifs.pgpd$Rds1p_R; motifs.pgpd$Rds1p_F=NULL; motifs.pgpd$Rds1p_R=NULL
  motifs.pgpd$Reb1p = motifs.pgpd$Reb1p_F | motifs.pgpd$Reb1p_R; motifs.pgpd$Reb1p_F=NULL; motifs.pgpd$Reb1p_R=NULL
  motifs.pgpd$Abf1p = motifs.pgpd$Abf1p_F | motifs.pgpd$Abf1p_R; motifs.pgpd$Abf1p_F=NULL; motifs.pgpd$Abf1p_R=NULL
  motifs.pgpd$Crz1p = motifs.pgpd$Crz1p_F | motifs.pgpd$Crz1p_R; motifs.pgpd$Crz1p_F=NULL; motifs.pgpd$Crz1p_R=NULL
  motifs.pgpd$Dal82p = motifs.pgpd$Dal82p_F | motifs.pgpd$Dal82p_R; motifs.pgpd$Dal82p_F=NULL; motifs.pgpd$Dal82p_R=NULL
  motifs.pgpd$Pdr3p = motifs.pgpd$Pdr3p_F | motifs.pgpd$Pdr3p_R; motifs.pgpd$Pdr3p_F=NULL; motifs.pgpd$Pdr3p_R=NULL
  motifs.pgpd$Hap2p = motifs.pgpd$Hap2p_F | motifs.pgpd$Hap2p_R; motifs.pgpd$Hap2p_F=NULL; motifs.pgpd$Hap2p_R=NULL
  motifs.pgpd$Ino4p = motifs.pgpd$Ino4p_F | motifs.pgpd$Ino4p_R; motifs.pgpd$Ino4p_F=NULL; motifs.pgpd$Ino4p_R=NULL
  motifs.pgpd$Stb4p = motifs.pgpd$Stb4p_F | motifs.pgpd$Stb4p_R; motifs.pgpd$Stb4p_F=NULL; motifs.pgpd$Stb4p_R=NULL
  motifs.pgpd$Sfl1p = motifs.pgpd$Sfl1p_F | motifs.pgpd$Sfl1p_R; motifs.pgpd$Sfl1p_F=NULL; motifs.pgpd$Sfl1p_R=NULL
  motifs.pgpd$Met31p = motifs.pgpd$Met31p_F | motifs.pgpd$Met31p_R; motifs.pgpd$Met31p_F=NULL; motifs.pgpd$Met31p_R=NULL
  motifs.pgpd$Adr1p = motifs.pgpd$Adr1p_F | motifs.pgpd$Adr1p_R; motifs.pgpd$Adr1p_F=NULL; motifs.pgpd$Adr1p_R=NULL
  motifs.pgpd$Aft2p = motifs.pgpd$Aft2p_F | motifs.pgpd$Aft2p_R; motifs.pgpd$Aft2p_F=NULL; motifs.pgpd$Aft2p_R=NULL
  motifs.pgpd$Rof1p = motifs.pgpd$Rof1p_F | motifs.pgpd$Rof1p_R; motifs.pgpd$Rof1p_F=NULL; motifs.pgpd$Rof1p_R=NULL
  motifs.pgpd$Yap5p = motifs.pgpd$Yap5p_F | motifs.pgpd$Yap5p_R; motifs.pgpd$Yap5p_F=NULL; motifs.pgpd$Yap5p_R=NULL
  motifs.pgpd$Ndt80p = motifs.pgpd$Ndt80p_F | motifs.pgpd$Ndt80p_R; motifs.pgpd$Ndt80p_F=NULL; motifs.pgpd$Ndt80p_R=NULL
  motifs.pgpd$Uga3p = motifs.pgpd$Uga3p_F | motifs.pgpd$Uga3p_R; motifs.pgpd$Uga3p_F=NULL; motifs.pgpd$Uga3p_R=NULL
  motifs.pgpd$Tea1p = motifs.pgpd$Tea1p_F | motifs.pgpd$Tea1p_R; motifs.pgpd$Tea1p_F=NULL; motifs.pgpd$Tea1p_R=NULL
  motifs.pgpd$Nrg1p = motifs.pgpd$Nrg1p_F | motifs.pgpd$Nrg1p_R; motifs.pgpd$Nrg1p_F=NULL; motifs.pgpd$Nrg1p_R=NULL
  motifs.pgpd$Yrr1p = motifs.pgpd$Yrr1p_F | motifs.pgpd$Yrr1p_R; motifs.pgpd$Yrr1p_F=NULL; motifs.pgpd$Yrr1p_R=NULL
  motifs.pgpd$Sko1p = motifs.pgpd$Sko1p_F | motifs.pgpd$Sko1p_R; motifs.pgpd$Sko1p_F=NULL; motifs.pgpd$Sko1p_R=NULL
  motifs.pgpd$Pdr1p = motifs.pgpd$Pdr1p_F | motifs.pgpd$Pdr1p_R; motifs.pgpd$Pdr1p_F=NULL; motifs.pgpd$Pdr1p_R=NULL
  motifs.pgpd$Sip4p = motifs.pgpd$Sip4p_F | motifs.pgpd$Sip4p_R; motifs.pgpd$Sip4p_F=NULL; motifs.pgpd$Sip4p_R=NULL
  motifs.pgpd$Msn2p = motifs.pgpd$Msn2p_F | motifs.pgpd$Msn2p_R; motifs.pgpd$Msn2p_F=NULL; motifs.pgpd$Msn2p_R=NULL
  motifs.pgpd$Gsm1p = motifs.pgpd$Gsm1p_F | motifs.pgpd$Gsm1p_R; motifs.pgpd$Gsm1p_F=NULL; motifs.pgpd$Gsm1p_R=NULL
  motifs.pgpd$Met4p = motifs.pgpd$Met4p_F | motifs.pgpd$Met4p_R; motifs.pgpd$Met4p_F=NULL; motifs.pgpd$Met4p_R=NULL
  motifs.pgpd$Yap1p = motifs.pgpd$Yap1p_F | motifs.pgpd$Yap1p_R; motifs.pgpd$Yap1p_F=NULL; motifs.pgpd$Yap1p_R=NULL
  motifs.pgpd$Mot3p = motifs.pgpd$Mot3p_F | motifs.pgpd$Mot3p_R; motifs.pgpd$Mot3p_F=NULL; motifs.pgpd$Mot3p_R=NULL
  motifs.pgpd$Stp2p = motifs.pgpd$Stp2p_F | motifs.pgpd$Stp2p_R; motifs.pgpd$Stp2p_F=NULL; motifs.pgpd$Stp2p_R=NULL
  motifs.pgpd$Ime1p = motifs.pgpd$Ime1p_F | motifs.pgpd$Ime1p_R; motifs.pgpd$Ime1p_F=NULL; motifs.pgpd$Ime1p_R=NULL
  motifs.pgpd$Rpn4p = motifs.pgpd$Rpn4p_F | motifs.pgpd$Rpn4p_R; motifs.pgpd$Rpn4p_F=NULL; motifs.pgpd$Rpn4p_R=NULL
  motifs.pgpd$Gzf3p = motifs.pgpd$Gzf3p_F | motifs.pgpd$Gzf3p_R; motifs.pgpd$Gzf3p_F=NULL; motifs.pgpd$Gzf3p_R=NULL
  motifs.pgpd$Upc2p = motifs.pgpd$Upc2p_F | motifs.pgpd$Upc2p_R; motifs.pgpd$Upc2p_F=NULL; motifs.pgpd$Upc2p_R=NULL
  motifs.pgpd$Swi4p = motifs.pgpd$Swi4p_F | motifs.pgpd$Swi4p_R; motifs.pgpd$Swi4p_F=NULL; motifs.pgpd$Swi4p_R=NULL
  motifs.pgpd$Zap1p = motifs.pgpd$Zap1p_F | motifs.pgpd$Zap1p_R; motifs.pgpd$Zap1p_F=NULL; motifs.pgpd$Zap1p_R=NULL
  motifs.pgpd$Stp1p = motifs.pgpd$Stp1p_F | motifs.pgpd$Stp1p_R; motifs.pgpd$Stp1p_F=NULL; motifs.pgpd$Stp1p_R=NULL
  motifs.pgpd$Sum1p = motifs.pgpd$Sum1p_F | motifs.pgpd$Sum1p_R; motifs.pgpd$Sum1p_F=NULL; motifs.pgpd$Sum1p_R=NULL
  motifs.pgpd$Hap5p = motifs.pgpd$Hap5p_F | motifs.pgpd$Hap5p_R; motifs.pgpd$Hap5p_F=NULL; motifs.pgpd$Hap5p_R=NULL
  motifs.pgpd$Rgt1p = motifs.pgpd$Rgt1p_F | motifs.pgpd$Rgt1p_R; motifs.pgpd$Rgt1p_F=NULL; motifs.pgpd$Rgt1p_R=NULL
  motifs.pgpd$Tda9p = motifs.pgpd$Tda9p_F | motifs.pgpd$Tda9p_R; motifs.pgpd$Tda9p_F=NULL; motifs.pgpd$Tda9p_R=NULL
  motifs.pgpd$Cbf1p = motifs.pgpd$Cbf1p_F | motifs.pgpd$Cbf1p_R; motifs.pgpd$Cbf1p_F=NULL; motifs.pgpd$Cbf1p_R=NULL
  motifs.pgpd$Haa1p = motifs.pgpd$Haa1p_F | motifs.pgpd$Haa1p_R; motifs.pgpd$Haa1p_F=NULL; motifs.pgpd$Haa1p_R=NULL
  motifs.pgpd$Cin5p = motifs.pgpd$Cin5p_F | motifs.pgpd$Cin5p_R; motifs.pgpd$Cin5p_F=NULL; motifs.pgpd$Cin5p_R=NULL
  motifs.pgpd$Pip2p = motifs.pgpd$Pip2p_F | motifs.pgpd$Pip2p_R; motifs.pgpd$Pip2p_F=NULL; motifs.pgpd$Pip2p_R=NULL
  motifs.pgpd$Gln3p = motifs.pgpd$Gln3p_F | motifs.pgpd$Gln3p_R; motifs.pgpd$Gln3p_F=NULL; motifs.pgpd$Gln3p_R=NULL
  motifs.pgpd$Flo8p = motifs.pgpd$Flo8p_F | motifs.pgpd$Flo8p_R; motifs.pgpd$Flo8p_F=NULL; motifs.pgpd$Flo8p_R=NULL
  motifs.pgpd$Stb5p = motifs.pgpd$Stb5p_F | motifs.pgpd$Stb5p_R; motifs.pgpd$Stb5p_F=NULL; motifs.pgpd$Stb5p_R=NULL
  motifs.pgpd$Hac1p = motifs.pgpd$Hac1p_F | motifs.pgpd$Hac1p_R; motifs.pgpd$Hac1p_F=NULL; motifs.pgpd$Hac1p_R=NULL
  motifs.pgpd$Xbp1p = motifs.pgpd$Xbp1p_F | motifs.pgpd$Xbp1p_R; motifs.pgpd$Xbp1p_F=NULL; motifs.pgpd$Xbp1p_R=NULL
  motifs.pgpd$Ume6p = motifs.pgpd$Ume6p_F | motifs.pgpd$Ume6p_R; motifs.pgpd$Ume6p_F=NULL; motifs.pgpd$Ume6p_R=NULL
  motifs.pgpd$Azf1p = motifs.pgpd$Azf1p_F | motifs.pgpd$Azf1p_R; motifs.pgpd$Azf1p_F=NULL; motifs.pgpd$Azf1p_R=NULL
  motifs.pgpd$Mig1p = motifs.pgpd$Mig1p_F | motifs.pgpd$Mig1p_R; motifs.pgpd$Mig1p_F=NULL; motifs.pgpd$Mig1p_R=NULL
  motifs.pgpd$Com2p = motifs.pgpd$Com2p_F | motifs.pgpd$Com2p_R; motifs.pgpd$Com2p_F=NULL; motifs.pgpd$Com2p_R=NULL
  motifs.pgpd$Ppr1p = motifs.pgpd$Ppr1p_F | motifs.pgpd$Ppr1p_R; motifs.pgpd$Ppr1p_F=NULL; motifs.pgpd$Ppr1p_R=NULL
  motifs.pgpd$Ixr1p = motifs.pgpd$Ixr1p_F | motifs.pgpd$Ixr1p_R; motifs.pgpd$Ixr1p_F=NULL; motifs.pgpd$Ixr1p_R=NULL
  motifs.pgpd$Arg81p = motifs.pgpd$Arg81p_F | motifs.pgpd$Arg81p_R; motifs.pgpd$Arg81p_F=NULL; motifs.pgpd$Arg81p_R=NULL
  motifs.pgpd$Tye7p = motifs.pgpd$Tye7p_F | motifs.pgpd$Tye7p_R; motifs.pgpd$Tye7p_F=NULL; motifs.pgpd$Tye7p_R=NULL
  motifs.pgpd$Znf1p = motifs.pgpd$Znf1p_F | motifs.pgpd$Znf1p_R; motifs.pgpd$Znf1p_F=NULL; motifs.pgpd$Znf1p_R=NULL
  motifs.pgpd$Oaf1p = motifs.pgpd$Oaf1p_F | motifs.pgpd$Oaf1p_R; motifs.pgpd$Oaf1p_F=NULL; motifs.pgpd$Oaf1p_R=NULL
  motifs.pgpd$Rfx1p = motifs.pgpd$Rfx1p_F | motifs.pgpd$Rfx1p_R; motifs.pgpd$Rfx1p_F=NULL; motifs.pgpd$Rfx1p_R=NULL
  motifs.pgpd$Hap4p = motifs.pgpd$Hap4p_F | motifs.pgpd$Hap4p_R; motifs.pgpd$Hap4p_F=NULL; motifs.pgpd$Hap4p_R=NULL
  motifs.pgpd$Pdr8p = motifs.pgpd$Pdr8p_F | motifs.pgpd$Pdr8p_R; motifs.pgpd$Pdr8p_F=NULL; motifs.pgpd$Pdr8p_R=NULL
  motifs.pgpd$Skn7p = motifs.pgpd$Skn7p_F | motifs.pgpd$Skn7p_R; motifs.pgpd$Skn7p_F=NULL; motifs.pgpd$Skn7p_R=NULL
  motifs.pgpd$Ste12p = motifs.pgpd$Ste12p_F | motifs.pgpd$Ste12p_R; motifs.pgpd$Ste12p_F=NULL; motifs.pgpd$Ste12p_R=NULL
  motifs.pgpd$Rap1p = motifs.pgpd$Rap1p_F | motifs.pgpd$Rap1p_R; motifs.pgpd$Rap1p_F=NULL; motifs.pgpd$Rap1p_R=NULL
  motifs.pgpd$Mig2p = motifs.pgpd$Mig2p_F | motifs.pgpd$Mig2p_R; motifs.pgpd$Mig2p_F=NULL; motifs.pgpd$Mig2p_R=NULL
  motifs.pgpd$Hap1p = motifs.pgpd$Hap1p_F | motifs.pgpd$Hap1p_R; motifs.pgpd$Hap1p_F=NULL; motifs.pgpd$Hap1p_R=NULL
  motifs.pgpd$Mal63p = motifs.pgpd$Mal63p_F | motifs.pgpd$Mal63p_R; motifs.pgpd$Mal63p_F=NULL; motifs.pgpd$Mal63p_R=NULL
  motifs.pgpd$Swi5p = motifs.pgpd$Swi5p_F | motifs.pgpd$Swi5p_R; motifs.pgpd$Swi5p_F=NULL; motifs.pgpd$Swi5p_R=NULL
  motifs.pgpd$Ace2p = motifs.pgpd$Ace2p_F | motifs.pgpd$Ace2p_R; motifs.pgpd$Ace2p_F=NULL; motifs.pgpd$Ace2p_R=NULL
  motifs.pgpd$Mac1p = motifs.pgpd$Mac1p_F | motifs.pgpd$Mac1p_R; motifs.pgpd$Mac1p_F=NULL; motifs.pgpd$Mac1p_R=NULL
  motifs.pgpd$Rph1p = motifs.pgpd$Rph1p_F | motifs.pgpd$Rph1p_R; motifs.pgpd$Rph1p_F=NULL; motifs.pgpd$Rph1p_R=NULL
  motifs.pgpd$Yap3p = motifs.pgpd$Yap3p_F | motifs.pgpd$Yap3p_R; motifs.pgpd$Yap3p_F=NULL; motifs.pgpd$Yap3p_R=NULL
  motifs.pgpd$Aft1p = motifs.pgpd$Aft1p_F | motifs.pgpd$Aft1p_R; motifs.pgpd$Aft1p_F=NULL; motifs.pgpd$Aft1p_R=NULL
  motifs.pgpd$Pho4p = motifs.pgpd$Pho4p_F | motifs.pgpd$Pho4p_R; motifs.pgpd$Pho4p_F=NULL; motifs.pgpd$Pho4p_R=NULL
  motifs.pgpd$Mcm1p = motifs.pgpd$Mcm1p_F | motifs.pgpd$Mcm1p_R; motifs.pgpd$Mcm1p_F=NULL; motifs.pgpd$Mcm1p_R=NULL
  motifs.pgpd$Gal4p = motifs.pgpd$Gal4p_F | motifs.pgpd$Gal4p_R; motifs.pgpd$Gal4p_F=NULL; motifs.pgpd$Gal4p_R=NULL
  motifs.pgpd$Hap3p = motifs.pgpd$Hap3p_F | motifs.pgpd$Hap3p_R; motifs.pgpd$Hap3p_F=NULL; motifs.pgpd$Hap3p_R=NULL
  motifs.pgpd$Fkh1p = motifs.pgpd$Fkh1p_F | motifs.pgpd$Fkh1p_R; motifs.pgpd$Fkh1p_F=NULL; motifs.pgpd$Fkh1p_R=NULL
  motifs.pgpd$Cad1p = motifs.pgpd$Cad1p_F | motifs.pgpd$Cad1p_R; motifs.pgpd$Cad1p_F=NULL; motifs.pgpd$Cad1p_R=NULL
  motifs.pgpd$Mot2p = motifs.pgpd$Mot2p_F | motifs.pgpd$Mot2p_R; motifs.pgpd$Mot2p_F=NULL; motifs.pgpd$Mot2p_R=NULL
  motifs.pgpd$Ecm22p = motifs.pgpd$Ecm22p_F | motifs.pgpd$Ecm22p_R; motifs.pgpd$Ecm22p_F=NULL; motifs.pgpd$Ecm22p_R=NULL
  motifs.pgpd$Gcr1p = motifs.pgpd$Gcr1p_F | motifs.pgpd$Gcr1p_R; motifs.pgpd$Gcr1p_F=NULL; motifs.pgpd$Gcr1p_R=NULL
  motifs.pgpd$Cha4p = motifs.pgpd$Cha4p_F | motifs.pgpd$Cha4p_R; motifs.pgpd$Cha4p_F=NULL; motifs.pgpd$Cha4p_R=NULL
  motifs.pgpd$Ino2p = motifs.pgpd$Ino2p_F | motifs.pgpd$Ino2p_R; motifs.pgpd$Ino2p_F=NULL; motifs.pgpd$Ino2p_R=NULL
  motifs.pgpd$Cat8p = motifs.pgpd$Cat8p_F | motifs.pgpd$Cat8p_R; motifs.pgpd$Cat8p_F=NULL; motifs.pgpd$Cat8p_R=NULL
  motifs.pgpd$Rox1p = motifs.pgpd$Rox1p_F | motifs.pgpd$Rox1p_R; motifs.pgpd$Rox1p_F=NULL; motifs.pgpd$Rox1p_R=NULL
  motifs.pgpd$Pho2p = motifs.pgpd$Pho2p_F | motifs.pgpd$Pho2p_R; motifs.pgpd$Pho2p_F=NULL; motifs.pgpd$Pho2p_R=NULL
  motifs.pgpd$Kar4p = motifs.pgpd$Kar4p_F | motifs.pgpd$Kar4p_R; motifs.pgpd$Kar4p_F=NULL; motifs.pgpd$Kar4p_R=NULL
  motifs.pgpd$Gcn4p = motifs.pgpd$Gcn4p_F | motifs.pgpd$Gcn4p_R; motifs.pgpd$Gcn4p_F=NULL; motifs.pgpd$Gcn4p_R=NULL
  motifs.pgpd$Nrg2p = motifs.pgpd$Nrg2p_F | motifs.pgpd$Nrg2p_R; motifs.pgpd$Nrg2p_F=NULL; motifs.pgpd$Nrg2p_R=NULL
  motifs.pgpd$Hcm1p = motifs.pgpd$Hcm1p_F | motifs.pgpd$Hcm1p_R; motifs.pgpd$Hcm1p_F=NULL; motifs.pgpd$Hcm1p_R=NULL
  motifs.pgpd$Fzf1p = motifs.pgpd$Fzf1p_F | motifs.pgpd$Fzf1p_R; motifs.pgpd$Fzf1p_F=NULL; motifs.pgpd$Fzf1p_R=NULL
  motifs.pgpd$Gat1p = motifs.pgpd$Gat1p_F | motifs.pgpd$Gat1p_R; motifs.pgpd$Gat1p_F=NULL; motifs.pgpd$Gat1p_R=NULL
  motifs.pgpd$Fkh2p = motifs.pgpd$Fkh2p_F | motifs.pgpd$Fkh2p_R; motifs.pgpd$Fkh2p_F=NULL; motifs.pgpd$Fkh2p_R=NULL
  motifs.pgpd$Mbp1p = motifs.pgpd$Mbp1p_F | motifs.pgpd$Mbp1p_R; motifs.pgpd$Mbp1p_F=NULL; motifs.pgpd$Mbp1p_R=NULL
  motifs.pgpd$Rlm1p = motifs.pgpd$Rlm1p_F | motifs.pgpd$Rlm1p_R; motifs.pgpd$Rlm1p_F=NULL; motifs.pgpd$Rlm1p_R=NULL
  motifs.pgpd$Leu3p = motifs.pgpd$Leu3p_F | motifs.pgpd$Leu3p_R; motifs.pgpd$Leu3p_F=NULL; motifs.pgpd$Leu3p_R=NULL
  motifs.pgpd$Rtg1p = motifs.pgpd$Rtg1p_F | motifs.pgpd$Rtg1p_R; motifs.pgpd$Rtg1p_F=NULL; motifs.pgpd$Rtg1p_R=NULL
  motifs.pgpd$Bas1p = motifs.pgpd$Bas1p_F | motifs.pgpd$Bas1p_R; motifs.pgpd$Bas1p_F=NULL; motifs.pgpd$Bas1p_R=NULL
  motifs.pgpd$Dal81p = motifs.pgpd$Dal81p_F | motifs.pgpd$Dal81p_R; motifs.pgpd$Dal81p_F=NULL; motifs.pgpd$Dal81p_R=NULL
  motifs.pgpd$Rim101p = motifs.pgpd$Rim101p_F | motifs.pgpd$Rim101p_R; motifs.pgpd$Rim101p_F=NULL; motifs.pgpd$Rim101p_R=NULL
  motifs.pgpd$Cst6p = motifs.pgpd$Cst6p_F | motifs.pgpd$Cst6p_R; motifs.pgpd$Cst6p_F=NULL; motifs.pgpd$Cst6p_R=NULL
  motifs.pgpd$Hsf1p = motifs.pgpd$Hsf1p_F | motifs.pgpd$Hsf1p_R; motifs.pgpd$Hsf1p_F=NULL; motifs.pgpd$Hsf1p_R=NULL
  motifs.pgpd$Hot1p = motifs.pgpd$Hot1p_F | motifs.pgpd$Hot1p_R; motifs.pgpd$Hot1p_F=NULL; motifs.pgpd$Hot1p_R=NULL
  motifs.pgpd$Cup2p = motifs.pgpd$Cup2p_F | motifs.pgpd$Cup2p_R; motifs.pgpd$Cup2p_F=NULL; motifs.pgpd$Cup2p_R=NULL
  motifs.pgpd$Aca1p = motifs.pgpd$Aca1p_F | motifs.pgpd$Aca1p_R; motifs.pgpd$Aca1p_F=NULL; motifs.pgpd$Aca1p_R=NULL
  motifs.pgpd$Put3p = motifs.pgpd$Put3p_F | motifs.pgpd$Put3p_R; motifs.pgpd$Put3p_F=NULL; motifs.pgpd$Put3p_R=NULL
  motifs.pgpd$Ash1p = motifs.pgpd$Ash1p_F | motifs.pgpd$Ash1p_R; motifs.pgpd$Ash1p_F=NULL; motifs.pgpd$Ash1p_R=NULL
  motifs.pgpd$Rtg3p = motifs.pgpd$Rtg3p_F | motifs.pgpd$Rtg3p_R; motifs.pgpd$Rtg3p_F=NULL; motifs.pgpd$Rtg3p_R=NULL
  motifs.pgpd$Arr1p = motifs.pgpd$Arr1p_F | motifs.pgpd$Arr1p_R; motifs.pgpd$Arr1p_F=NULL; motifs.pgpd$Arr1p_R=NULL
  motifs.pgpd$Lys14p = motifs.pgpd$Lys14p_F | motifs.pgpd$Lys14p_R; motifs.pgpd$Lys14p_F=NULL; motifs.pgpd$Lys14p_R=NULL
  motifs.pgpd$Arg80p = motifs.pgpd$Arg80p_F | motifs.pgpd$Arg80p_R; motifs.pgpd$Arg80p_F=NULL; motifs.pgpd$Arg80p_R=NULL
  motifs.pgpd$Smp1p = motifs.pgpd$Smp1p_F | motifs.pgpd$Smp1p_R; motifs.pgpd$Smp1p_F=NULL; motifs.pgpd$Smp1p_R=NULL
  motifs.pgpd$Mig3p = motifs.pgpd$Mig3p_F | motifs.pgpd$Mig3p_R; motifs.pgpd$Mig3p_F=NULL; motifs.pgpd$Mig3p_R=NULL
  motifs.pgpd$Rme1p = motifs.pgpd$Rme1p_F | motifs.pgpd$Rme1p_R; motifs.pgpd$Rme1p_F=NULL; motifs.pgpd$Rme1p_R=NULL
  motifs.pgpd$Met32p = motifs.pgpd$Met32p_F | motifs.pgpd$Met32p_R; motifs.pgpd$Met32p_F=NULL; motifs.pgpd$Met32p_R=NULL
  motifs.pgpd$Tec1p = motifs.pgpd$Tec1p_F | motifs.pgpd$Tec1p_R; motifs.pgpd$Tec1p_F=NULL; motifs.pgpd$Tec1p_R=NULL
  motifs.pgpd$Usv1p = motifs.pgpd$Usv1p_F | motifs.pgpd$Usv1p_R; motifs.pgpd$Usv1p_F=NULL; motifs.pgpd$Usv1p_R=NULL
  return(motifs.pgpd)
}


