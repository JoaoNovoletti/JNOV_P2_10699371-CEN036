#!usr/bin/env python3

#Importando os módulos necessários 
import re
import sys

seqs={}
codons={}
arquivo = sys.argv[1]

#Abrinido o arquivo e correspondendo o gene e a sequencia em um dicionário
with open(arquivo, "r") as R_arquivo:
  for line in R_arquivo:
    line = line.rstrip()

    if line[0] == '>':

     achandogene = re.search(r">(\w+)\s", line)
     gene = (achandogene.group(1))

     if gene not in seqs.keys():
       seqs[gene]=''

    else:
      seqs[gene] += line

#Montando um dicionário com os genes e seus 6 frames
for gene in seqs.keys():
  if gene not in codons.keys():
    codons[gene]={}
for gene in codons:
  codons[gene]["frame1"] = re.findall(r"(.{3})",seqs[gene])
  codons[gene]["frame2"] = re.findall(r"(.{3})",seqs[gene][1:])
  codons[gene]["frame3"] = re.findall(r"(.{3})",seqs[gene][2:])
  codons[gene]["frame4"] = re.findall(r"(.{3})",seqs[gene][::-1])
  codons[gene]["frame5"] = re.findall(r"(.{3})",seqs[gene][-2::-1])
  codons[gene]["frame6"] = re.findall(r"(.{3})",seqs[gene][-3::-1])

#Achando os maiores ORFS de cada gene comparando todos os frames, imprimindo em ORF.fna e armazenando os dados para uso posterior
dadoscodonpepmaior = {}

with open ("ORF.fna", "w") as Tlongest:
  for gene in codons:
    maior = 0
    codonpepmaior = "Orf not found in any frame"
    tamanho = 0
    sequencia = ""
    for frame in codons[gene]:
      framejunto = 'x'.join(codons[gene][frame])
      
      framejuntotrocado = re.sub (r"TAG|TGA|TAA", "***", framejunto)
      
      for codonpep in re.finditer(r"(ATG\w*\*\*\*)", framejuntotrocado):
        
        tamanho = len(codonpep.group(0))
        if tamanho > maior:
          maior = tamanho
          sequencia = frame
          inicio = codonpep.start(0)
          fim = codonpep.end(0)
          codonpepmaior = framejunto[inicio:fim]
          
          codonpepmaiorwrite = codonpepmaior.replace("x", "")
          
    if sequencia != "":      
      framesequencia = ''.join(codons[gene][sequencia])
      achandoposicao = re.search(rf"{codonpepmaiorwrite}", framesequencia)
      inicioorf = achandoposicao.start(0)+1
      fimorf = achandoposicao.end(0)
      
      Tlongest.write(">" + gene + "_" + sequencia + "_" + str(inicioorf) + "_" + str(fimorf) + '\n' + codonpepmaiorwrite + '\n')
    else:
      Tlongest.write(">" + gene + "_" + '\n' + codonpepmaior + '\n')
    
    dadoscodonpepmaior[gene] = {}
    dadoscodonpepmaior[gene][sequencia] = []
    if sequencia != "": 
      dadoscodonpepmaior[gene][sequencia].append(codonpepmaior)
      dadoscodonpepmaior[gene][sequencia].append(inicioorf)
      dadoscodonpepmaior[gene][sequencia].append(fimorf)


#Montando um dicionário com as traduções dos códons
translation_table = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'AAT':'N', 'AAC':'N',
    'GAT':'D', 'GAC':'D',
    'TGT':'C', 'TGC':'C',
    'CAA':'Q', 'CAG':'Q',
    'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'CAT':'H', 'CAC':'H',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'AAA':'K', 'AAG':'K',
    'ATG':'M',
    'TTT':'F', 'TTC':'F',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'TGG':'W',
    'TAT':'Y', 'TAC':'Y',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TAA':'*', 'TGA':'*', 'TAG':'*'
}

#Abrindo o arquivo Python_08.translated.aa e escrevendo nele a tradução dos códons
with open ("ORF.faa", "w") as traducao:
  for gene in dadoscodonpepmaior:
    for sequencia in dadoscodonpepmaior[gene]:
      if sequencia != "":
        line = dadoscodonpepmaior[gene][sequencia][0]
        line = line.rstrip()
        codonslist = line.split('x')
        translatedcodons = []
        for codon in codonslist:
          codonfound = re.search(r"[ATCG][ATCG][ATCG]", codon)
          if codonfound:
            translatedcodons.append(translation_table[codon])
      translatedlist = ''
      for codons in translatedcodons:
        translatedlist += codons

    if sequencia != "":
      traducao.write(">" + gene + "_" + sequencia + "_" + str(dadoscodonpepmaior[gene][sequencia][1]) + "_" + str(dadoscodonpepmaior[gene][sequencia][2]) + '\n' + translatedlist + '\n')
    else:
      traducao.write(">" + gene + "_" + '\n' + "Orf not found in any frame" + '\n')
