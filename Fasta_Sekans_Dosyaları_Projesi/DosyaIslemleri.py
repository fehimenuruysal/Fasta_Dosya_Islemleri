import Needleman_Wunsch

class FASTA_DOSYASI():
    def __init__(self,dosya_adi,dosya_yolu):
        self.dosya_adi = dosya_adi
        self.dosya_yolu = dosya_yolu


    def genomDosyaDonusturucu(self):
        fasta_dosyasi = open(self.dosya_yolu, "r")
        data = fasta_dosyasi.read()
        sekans_baslangic = data.find("genome")
        dosya = open('Cikti_Dosyasi\\'+ self.dosya_adi +'.txt', "w")
        dosya.write(data[sekans_baslangic+6:])

    def proteinDosyaDonusturucu(self):
        fasta_dosyasi = open(self.dosya_yolu, "r")
        data = fasta_dosyasi.read()
        protein_dizilim_baslangic = data.find("]")
        dosya = open('Cikti_Dosyasi\\' + self.dosya_adi + '.txt', "w")
        dosya.write(data[protein_dizilim_baslangic+1:])

    def dosyaOkuma(self):
        dosya = open(self.dosya_yolu, "r")
        sekanslar = dosya.read()
        sekanslar = sekanslar.replace("\n", "")
        sekanslar = sekanslar.replace("\r", "")
        self.orfBulNukleotit(sekanslar)
        self.donusturucu(sekanslar)

    def donusturucu(self,sekanslar):
        amino_asit_tablosu = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '-', 'TAG': '-',
            'TGC': 'C', 'TGT': 'C', 'TGA': '-', 'TGG': 'W',
        }
        protein = ""
        for i in range(0, len(sekanslar), 3):
            if len(sekanslar) > i+3:
                kodon = sekanslar[i:i + 3]
                protein += amino_asit_tablosu[kodon]
        self.orfBulAminoAsit(protein)
        protein += "\n\n\n"
        for j in range(1, len(sekanslar) - 2, 3):
            if len(sekanslar) > (j + 3):
                kodon = sekanslar[j:j + 3]
                protein += amino_asit_tablosu[kodon]
        protein += "\n\n\n"
        for k in range(2, len(sekanslar) - 1, 3):
            if len(sekanslar) > (k + 3):
                kodon = sekanslar[k:k + 3]
                protein += amino_asit_tablosu[kodon]
        self.dosyaYazma(protein)

    def dosyaYazma(self,amino_asit_sekansi):
        dosya = open('Cikti_Dosyasi\\'+ self.dosya_adi +'_AA.txt', 'w')
        dosya.write(amino_asit_sekansi)

    def orfBulAminoAsit(self, protein):
        orf_protein = ""
        dosya = open('Cikti_Dosyasi\\'+ self.dosya_adi +'_AA_ORF.txt', 'w')
        for i in range(0, len(protein) - 1, 1):
            if (protein[i] == 'M'):
                i += 1
                while (protein[i] != '-' and i < len(protein) - 1):
                    orf_protein += protein[i]
                    i += 1
                if (protein[i] == '-'):
                    dosya.write(orf_protein + '\n\n')
                orf_protein = ""
        dosya.write(orf_protein)

    def orfBulNukleotit(self, sekanslar):
        orf_sekans= []
        orf_gen= ""
        dosya = open('Cikti_Dosyasi\\'+ self.dosya_adi +'_GEN_ORF.txt', 'w')

        for i in range(0, len(sekanslar)-1, 3):
            if len(sekanslar) > i+3:
                kodon = sekanslar[i:i + 3]
                orf_sekans.append(kodon)

        for i in range(0,len(orf_sekans)-1,1) :
            if orf_sekans[i] == "ATG" and (i+1) < len(orf_sekans):
                i += 1
                while (orf_sekans[i] != "TAA") and (orf_sekans[i] != "TAG") and (orf_sekans[i] !="TGA"):
                    orf_gen += orf_sekans[i]
                    if (i + 1) < len(orf_sekans):
                        i+=1
                    else:
                        break
                orf_gen += "\n\n\n"
        dosya.write(orf_gen)

def main():
    fasta_sars_nucleocapsid = FASTA_DOSYASI("SARS_Nucleocapsid", "SARS_Nucleocapsid.fasta")
    fasta_sars_nucleocapsid.proteinDosyaDonusturucu()

    fasta_sars_polyprotein = FASTA_DOSYASI("SARS_Polyprotein","SARS_Polyprotein.fasta")
    fasta_sars_polyprotein.proteinDosyaDonusturucu()

    fasta_covid_nucleocapsid = FASTA_DOSYASI("COVID_Nucleocapsid","COVID_Nucleocapsid.fasta")
    fasta_covid_nucleocapsid.proteinDosyaDonusturucu()

    fasta_covid_polyprotein = FASTA_DOSYASI("COVID_Polyprotein","COVID_Polyprotein.fasta")
    fasta_covid_polyprotein.proteinDosyaDonusturucu()

    fasta_sars = FASTA_DOSYASI("SARS", "SARS.fasta")
    fasta_sars.genomDosyaDonusturucu()

    fasta_covid = FASTA_DOSYASI("COVID", "COVID.fasta")
    fasta_covid.genomDosyaDonusturucu()

    sars = FASTA_DOSYASI("SARS", "Cikti_Dosyasi\\SARS.txt")
    sars.dosyaOkuma()

    covid = FASTA_DOSYASI("COVID", "Cikti_Dosyasi\\COVID.txt")
    covid.dosyaOkuma()

    Needleman_Wunsch.farkBul("Cikti_Dosyasi\COVID_Nucleocapsid.txt", "Cikti_Dosyasi\\SARS_Nucleocapsid.txt", "Nucleocapsid")
    Needleman_Wunsch.farkBul("Cikti_Dosyasi\COVID_Polyprotein.txt","Cikti_Dosyasi\\SARS_Polyprotein.txt","Polyprotein")

    Needleman_Wunsch.needlemanWunsch("Cikti_Dosyasi\COVID_Nucleocapsid.txt", "Cikti_Dosyasi\\SARS_Nucleocapsid.txt", "Nucleocapsid")
    Needleman_Wunsch.needlemanWunsch("Cikti_Dosyasi\COVID_Polyprotein.txt","Cikti_Dosyasi\\SARS_Polyprotein.txt","Polyprotein")

if __name__ == "__main__":
    main()
