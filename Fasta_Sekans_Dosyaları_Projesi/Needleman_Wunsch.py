import BLOSUM62
import numpy as np

def needlemanWunsch(covid , sars , protein_turu):
    bosluk_ceza_puani = -4
    covid_protein ="-"+proteinDizisiOku(covid)
    sars_protein ="-"+proteinDizisiOku(sars)
    sat = len(covid_protein)
    sut = len(sars_protein)

    skor_matrisi = np.zeros([sat , sut])
    yol_matrisi = np.zeros([sat,sut])

    for i in range(0,sut,1):
        skor_matrisi[0][i] = i * bosluk_ceza_puani
    for j in range(0,sat,1):
        skor_matrisi[j][0] = j * bosluk_ceza_puani

    # Yol Matrisi : 0 bitis , 1 sol , 2 yukarı , 3 diagonal
    yol_matrisi[0][0] = 0

    for i in range(1,sut,1):
        yol_matrisi[0][i] = 1

    for j in range(1,sat,1):
        yol_matrisi[j][0] = 2

    for satir in range(1,sat,1):
        for sutun in range(1,sut,1):
            if covid_protein[satir] == sars_protein[sutun]:
                skor_matrisi[satir][sutun]=skor_matrisi[satir-1][sutun-1]+BLOSUM62.skorBul(covid_protein[satir],sars_protein[sutun])
                yol_matrisi[satir][sutun]= 3
            else:
                skor_diagonal = skor_matrisi[satir - 1][sutun - 1] + BLOSUM62.skorBul(covid_protein[satir], sars_protein[sutun])
                skor_yukari = skor_matrisi[satir-1][sutun] + BLOSUM62.skorBul(covid_protein[satir],sars_protein[sutun])
                skor_sol = skor_matrisi[satir][sutun-1] + BLOSUM62.skorBul(covid_protein[satir],sars_protein[sutun])
                if skor_yukari >= skor_sol and skor_yukari >= skor_diagonal:
                    skor_matrisi[satir][sutun] = skor_sol
                    yol_matrisi[satir][sutun] = 2
                elif skor_sol >= skor_yukari and skor_sol >= skor_diagonal:
                    skor_matrisi[satir][sutun] = skor_yukari
                    yol_matrisi[satir][sutun] = 1
                else:
                    skor_matrisi[satir][sutun] = skor_diagonal
                    yol_matrisi[satir][sutun] = 3
    hizalama_sonucu_covid = ""
    hizalama_sonucu_sars = ""
    yol_sutun = sut-1
    yol_satir = sat-1
    adim = yol_matrisi[yol_satir][yol_sutun]

    while(adim != 0 ):
        if adim ==1:
            hizalama_sonucu_sars += sars_protein[yol_sutun]
            hizalama_sonucu_covid += "-"
            yol_sutun -= 1
        elif adim == 2:
            hizalama_sonucu_sars += "-"
            hizalama_sonucu_covid += covid_protein[yol_satir]
            yol_satir -= 1
        elif adim == 3:
            hizalama_sonucu_covid += covid_protein[yol_satir]
            hizalama_sonucu_sars += sars_protein[yol_sutun]
            yol_satir -= 1
            yol_sutun -= 1
        adim = yol_matrisi[yol_satir][yol_sutun]

    hizalama_dosyasi = open("Cikti_Dosyasi\\"+protein_turu+"_NWA.txt" , "w")
    hizalama_dosyasi.write(protein_turu+ " Needleman Wunsch Hizalaması \n\n\n" + "nCovid :" + hizalama_sonucu_covid[::-1] + "\n\nSARS   :" + hizalama_sonucu_sars[::-1] )

def proteinDizisiOku(dosya_yolu):
    dosya = open(dosya_yolu, "r")
    protein = dosya.read()
    protein = protein.replace("\n", "")
    protein = protein.replace("\r", "")
    return protein

def farkBul(covid , sars , protein_turu):
    covid_protein = proteinDizisiOku(covid)
    sars_protein = proteinDizisiOku(sars)

    covid_uzunluk = len(covid_protein)-1
    sars_uzunluk = len(sars_protein)-1

    maX_uzunluk = max(covid_uzunluk , sars_uzunluk)

    dosya_metni = ""

    i = 0
    while i<maX_uzunluk:

        fark_metni=""
        fark_sars =""
        fark_nCovid = ""

        if i<covid_uzunluk and i<sars_uzunluk :
            if covid_protein[i] != sars_protein[i]:
                fark_metni += "Konum : " +str(i)

                while covid_protein[i] != sars_protein[i] and i<covid_uzunluk and i<sars_uzunluk:
                    fark_sars += sars_protein[i]
                    fark_nCovid += covid_protein[i]
                    i += 1
                fark_metni += "\nnCovid = " + fark_nCovid + "\nSARS   = " + fark_sars + "\n\n**************\n\n"
            else:
                i += 1

        elif i>=covid_uzunluk and i<=sars_uzunluk:
            fark_metni += "Konum : " +str(i)
            while i<sars_uzunluk:
                fark_sars += sars_protein[i]
                fark_nCovid += "-"
                i += 1
            fark_metni += "\nnCovid = " + fark_nCovid + "\nSARS   = " + fark_sars + "\n\n**************\n\n"

        elif i>=sars_uzunluk and i<=covid_uzunluk :
            fark_metni += "Konum : " +str(i) + "\n"

            while i<covid_uzunluk:
                fark_sars += "-"
                fark_nCovid += covid_protein[i]
                i += 1
            fark_metni += "\nnCovid = " + fark_nCovid + "\nSARS   = " + fark_sars + "\n\n**************\n\n"
        i+=1

        if fark_metni != "":
            print(fark_metni)
            dosya_metni += fark_metni
    dosya = open("Cikti_Dosyasi\\Fark_"+ protein_turu + ".txt" , "w")
    dosya.write(dosya_metni)