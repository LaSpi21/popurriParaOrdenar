# -*- coding: utf-8 -*-
"""
Created on Fri May  6 17:36:47 2022

@author: plasp
"""

data= []
with open("H:/python/cosas labo/data_for_jpreds_disome.txt","r",encoding = "utf-8", newline= "\r\n") as arch_data:
    lineas = arch_data.readlines()
    for linea in lineas:
        data.append(list(linea.split("_")))
        
#%%

dict ={}
for hit in data:
    if hit[1] in dict.keys():
        dict[hit[1]].append([hit[2],int(hit[3].strip("\r\n")),float(hit[0])], hit[1])
    else:
        dict[hit[1]] = [[hit[2],int(hit[3].strip("\r\n")),float(hit[0]), hit[1]]]
        
#%%


with open("H:/python/cosas labo/uniprot-proteome-human-reviewed.fasta","r",encoding = "utf-8", newline= "\n") as proteome_data:
    lines = proteome_data.readlines()
    for key in dict.keys():
        for prot in dict[key]:
            print(key)
            for i in range(len(lines)):
                n = 0
                if "GN=" + prot[0] + " " in lines[i]:
                    if lines[i+1][prot[1]-2:prot[1]+1] == key:
                        if prot[1] >= 50:
                            sec = lines[i+1][prot[1]-50:prot[1]+1]
                        else:
                            sec = lines[i+1][:prot[1]+1]
                        prot.append(sec.strip("\n"))
                    elif lines[i+1][prot[1]-1:prot[1]+2] == key:
                        if prot[1] >= 49:
                            sec = lines[i+1][prot[1]-49:prot[1]+2]
                        else:
                            sec = lines[i+1][:prot[1]+2]
                        prot.append(sec.strip("\n"))
    
                    else: 
                        print("check " + prot[0] + " " + lines[i+1][prot[1]-5:prot[1]-2] + " " + lines[i+1][prot[1]-2:prot[1]+1] + " " +lines[i+1][prot[1]+1:prot[1]+4])


#%%

def Submit_prots(prot, data):
    import jpredapi
    output = jpredapi.submit('single', 'raw', seq = prot)
    return output.headers['Location'][-10:], data

def get_Jpred(job):
    import jpredapi
    import shutil
    jpredapi.get_results(jobid = job, results_dir_path="H:/python/cosas labo/jpred_sspred_disome/results", extract=True)
    try:
        arch_ = job + '.jnet'
        arch_ = 'H:/python/cosas labo/jpred_sspred_disome/results/' + job + "/" +  arch_
        shutil.move(arch_, "H:/python/cosas labo/jpred_sspred_disome")
        shutil.rmtree('H:/python/cosas labo/jpred_sspred_disome/results/' + job, ignore_errors=True)
    except:
        print('revisar ' + job)
    
def get_Jpred2(job):
    arch_ = job + '.jnet'
    arch_ = 'H:/python/cosas labo/jpred_sspred_disome' + arch_
    with open(arch_, 'r',encoding="utf-8") as simple:
        simple.readline()
        linea = simple.readline()
        pred_ss = linea[9:].replace(',', '')
        linea = simple.readline()
        CONF = linea[9:].replace(',', '')
        linea = simple.readline()
        pred_expose = linea[10:]
    return pred_ss, pred_expose, CONF

def get_Jpred3(job):
    arch_ = job + '.jnet'
    arch_ = 'H:/python/cosas labo/jpred_sspred_disome/' + arch_
    with open(arch_, 'r',encoding="utf-8") as simple:
        simple.readline()
        linea = simple.readline()
        pred_ss = linea[9:].replace(',', '')
        linea = simple.readline()
        CONF = linea[9:].replace(',', '')
        linea = simple.readline()
        pred_expose = linea[10:]
    return pred_ss, pred_expose, CONF

#%%
job_IDs = []
for key in dict.keys():
    print(key)
    for prot in dict[key]:
        
        if len(prot) == 5:
            jobs, data = Submit_prots(prot[3], prot)
            job_IDs.append([jobs, data])
            

#%%
job_IDs
    #%%
    
"""ESPERAR UN POCO Y DESPUES RECUPERAR LOS RESULTADOS"""
for i in range(0, len(job_IDs)):
    get_Jpred(job_IDs[i][0])
    
    
#%%

lista_jobs_no_terminados = []
for job in job_IDs:
    arch_ = job[0] + '.jnet'
    arch_ = 'H:/python/cosas labo/jpred_sspred_disome/' + arch_
    try:
        Arch = open(arch_, 'r',encoding="utf-8")
        Arch.close()
    except:
        print(job)
        lista_jobs_no_terminados.append(job)
        
#%%

predicciones_de_ss = []
predicciones_de_exposicion = []
Confianza = []
data = []
for job in job_IDs:
    try:
        pred_ss, pred_exp, conf = get_Jpred3(job[0])
        predicciones_de_ss.append(pred_ss)
        predicciones_de_exposicion.append(pred_exp)
        Confianza.append(conf)
        data.append(job[1])
    except:
        pass
    
    #%%

structure_aling = []

for i in range(len(data)):
    if data[i][1] > 30:
        structure_aling.append(predicciones_de_ss[i][:31])
    else:
        print(data[i][3][:data[i][1]+1])
        structure_aling.append(predicciones_de_ss[i][:data[i][1]+1])

            
#%%


    
['jp_fKmUHJd', ['TMEM38B', 284, 141.8285714, 'EKKSEAKSPSNGVGSLASKPVDVASDNVKKKHTKKNE\n']]
['jp_7qIBKbY', ['DIS3', 951, 137.4462541, 'KIRMSLVEPQIPGISIPTDTSNMDLNGPKKKKMKLGK\n']]
['jp_a2boWX9', ['PWP1', 244, 96.71559633, 'PVIEVWDLDIVDSLEPVFTLGSKLSKKKKKKGKKSSSAEGHTDAVLDLSWNKLIRNVLAS']]
['jp_sE_zvdD', ['SRRM2', 191, 84.0610687, 'QPAPEPPKPYSLVRESSSSRSPTPKQKKKKKKKDRGRRSESSSPRRERKKSSKKKKHRSE']]
['jp_9GDhEg7', ['CHD4', 115, 82.85130112, 'SGEGPEFVEEEEEVALRSDSEGSDYTPGKKKKKKLGPKKEKKSKSKRKEEEEEEDDDDDS']]
['jp_ve3vG2j', ['SAE1', 209, 78.40909091, 'KVAKVSQGVEDGPDTKRAKLDSSETTMVKKKVVFCPVKEALEVDWSSEKAKAALKRTTSD']]
['jp_TpXEgaA', ['ECI2', 160, 77.64705882, 'GTDRKSTGFETLVVTSEDGITKIMFNRPKKKNAINTEMYHEIMRALKAASKDDSIITVLT']]
['jp_7SkkiR6', ['RPS3A', 19, 73.07989691, 'MAVGKNKRLTKGGKKGAKKKVVDPFSKKDWYDVKAPAMFNIRNIGKTLV']]
['jp_18OGwOV', ['SRRM2', 192, 67.24885496, 'PAPEPPKPYSLVRESSSSRSPTPKQKKKKKKKDRGRRSESSSPRRERKKSSKKKKHRSES']]
['jp_VMgGD0x', ['NIFK', 215, 64.3125, 'PSLILQKTESISKTNRQTSTKGQVLRKKKKKVSGTLDTPEKTVDSQGPTPVCTPTFLERR']]
['jp_3rznapo', ['SRRM2', 215, 63.04580153, 'KQKKKKKKKDRGRRSESSSPRRERKKSSKKKKHRSESESKKRKHRSPTPKSKRKSKDKKR']]
['jp_CZmteMr', ['SRRM2', 193, 58.84274809, 'APEPPKPYSLVRESSSSRSPTPKQKKKKKKKDRGRRSESSSPRRERKKSSKKKKHRSESE']]
['jp_eDxi16a', ['ZCRB1', 138, 57.83673469, 'ECGESGHLSYACPKNMLGEREPPKKKEKKKKKKAPEPEEEIEEVEESEDEGEDPALDSLS']]
['jp_S1TUr3L', ['CHD4', 350, 57.63568773, 'DDASINSYSVSDGSTSRSSRSRKKLRTTKKKKKGEEEVTAVDGYETDHQDYCEVCQQGGE']]
['jp_2kTIMxD', ['DIS3', 952, 56.22801303, 'IRMSLVEPQIPGISIPTDTSNMDLNGPKKKKMKLGK\n']]
['jp_WCEmFY1', ['EIF5', 249, 56.0, 'SDHAKVLTLSDDLERTIEERVNILFDFVKKKKEEGVIDSSDKEIVAEAERLDVKAMGPLV']]
['jp_Ih88PpF', ['SF3B2', 334, 56.0, 'GQSASETEEDTVSVSKKEKNRKRRNRKKKKKPQRVRGVSSESSGDREKDSTRSRGSDSPA']]
['jp_0LQLl8v', ['EIF2S2', 83, 55.89387755, 'LEADEEDTRKKDASDDLDDLNFFNQKKKKKKTKKIFDIDEAEEGVKDLKIESDVQEPTEP']]
['jp_Wyn4REt', ['TCOF1', 1479, 52.32907348, 'KDKEKKEKKKKAKKASTKDSESPSQKKKKKKKKTAEQTV\n']]
['jp_lKqhnfh', ['TCOF1', 1481, 52.32907348, 'KEKKEKKKKAKKASTKDSESPSQKKKKKKKKTAEQTV\n']]
['jp_sRvPYqP', ['EIF2S2', 82, 47.71428571, 'DLEADEEDTRKKDASDDLDDLNFFNQKKKKKKTKKIFDIDEAEEGVKDLKIESDVQEPTE']]
['jp_DQA83qO', ['TCOF1', 1363, 47.57188498, 'TKESSRKGWESRKRKLSGDQPAARTPRSKKKKKLGAGEGGEASVSPEKTSTTSKGKAKRD']]
['jp_6IEdInw', ['SF3B2', 333, 46.34482759, 'LGQSASETEEDTVSVSKKEKNRKRRNRKKKKKPQRVRGVSSESSGDREKDSTRSRGSDSP']]
['jp_OWX98zD', ['DDX18', 426, 42.21523179, 'VDGLEQGYVVCPSEKRFLLLFTFLKKNRKKKLMVFFSSCMSVKYHYELLNYIDLPVLAIH']]
['jp_EqAnr4p', ['TMEM263', 111, 41.6, 'LVKGGVSAVAGGVTAVGSAVVNKVPLTGKKKDKSD\n']]
['jp_RYjZQEk', ['MTDH', 446, 40.36153846, 'EPIPDDQKVSDDDKEKGEGALPTGKSKKKKKKKKKQGEDNSTAQDTEELEKEIREDLPVN']]
['jp_13F_c8_', ['TAF7', 154, 38.59315589, 'DKEKKFIWNHGITLPLKNVRKRRFRKTAKKKYIESPDVEKEVKRLLSTDAEAVSTRWEII']]
['jp_eqbAmRs', ['IQCB1', 438, 35.80978261, 'QQRQSLIEYKAAVTLQRAALKFLAKCRKKKKLFAPWRGLQELTDARRVELKKRVDDYVRR']]
['jp_VBHMbzK', ['RPL3', 127, 33.02725021, 'GLRTFKTVFAEHISDECKRRFYKNWHKSKKKAFTKYCKKWQDEDGKKQLEKDFSSMKKYC']]
['jp_GWGum_Q', ['CKAP2', 155, 32.68867925, 'TPHLLLTEDDPQSQHMTLSQAFHLKNNSKKKQMTTEKQKQDANMPKKPVLGSYRGQIVQS']]
['jp_eEAi39o', ['CHD4', 114, 32.42007435, 'SGEGPEFVEEEEEVALRSDSEGSDYTPGKKKKKKLGPKKEKKSKSKRKEEEEEEDDDDDS']]
['jp_LpPTzuh', ['RPS27L', 21, 31.95652174, 'MPLARDLLHPSLEEEKKKHKKKRLVQSPNSYFMDVKCPGCYKITTVFSHAQ']]
['jp_mrF1vVF', ['NKAP', 199, 31.39622642, 'KSTTSASTSEEEKKKKSSRSKERSKKRRKKKSSKRKHKKYSEDSDSDSDSETDSSDEDNK']]
['jp_nPrABOY', ['MTDH', 447, 31.39230769, 'PIPDDQKVSDDDKEKGEGALPTGKSKKKKKKKKKQGEDNSTAQDTEELEKEIREDLPVNT']]
['jp_ercnj8r', ['GTF2F1', 236, 29.32075472, 'IHDLEDDLEMSSDASDASGEEGGRVPKAKKKAPLAKGGRKKKKKKGSDDEAFEDSDDGDF']]
['jp_U8M2yWn', ['RPS27A', 82, 29.19946092, 'GRTLSDYNIQKESTLHLVLRLRGGAKKRKKKSYTTPKKNKHKRKKVKLAVLKYYKVDENG']]
['jp_EB0_D8O', ['EIF5', 250, 28.0, 'DHAKVLTLSDDLERTIEERVNILFDFVKKKKEEGVIDSSDKEIVAEAERLDVKAMGPLVL']]
['jp_Z2hnSi8', ['PWP1', 241, 27.63302752, 'NMTPVIEVWDLDIVDSLEPVFTLGSKLSKKKKKKGKKSSSAEGHTDAVLDLSWNKLIRNV']]
['jp_1vQwgTC', ['NIFK', 213, 25.265625, 'DFPSLILQKTESISKTNRQTSTKGQVLRKKKKKVSGTLDTPEKTVDSQGPTPVCTPTFLE']]
['jp_Ms9jOKO', ['SUPT16H', 423, 24.73156342, 'PEEKTYALFIGDTVLVDEDGPATVLTSVKKKVKNVGIFLKNEDEEEEEEEKDEAEDLLGR']]
['jp_9DWb1qJ', ['PWP1', 242, 23.02752294, 'MTPVIEVWDLDIVDSLEPVFTLGSKLSKKKKKKGKKSSSAEGHTDAVLDLSWNKLIRNVL']]
['jp_4LeO4X1', ['NIFK', 214, 22.96875, 'FPSLILQKTESISKTNRQTSTKGQVLRKKKKKVSGTLDTPEKTVDSQGPTPVCTPTFLER']]
['jp_zEpKMzt', ['EIF2S2', 20, 22.49387755, 'MSGDEMIFDPTMSKKKKKKKKPFMLDEEGDTQTEETQPSETKEVEPEPTE']]
['jp_yVuE3ac', ['SON', 1151, 194.5136612, 'RSMMSMAADSYTDSYTDTYTEAYMVPPLPPEEPPTMPPLPPEEPPMTPPLPPEEPPEGPA']]
['jp_1iIPIOI', ['SDHA', 373, 151.3448276, 'SRSMTLEIREGRGCGPEKDHVYLQLHHLPPEQLATRLPGISETAMIFAGVDVTKEPIPVL']]
['jp_a8c7726', ['PRRC2A', 399, 102.0507497, 'PPEADGKKGNSPNSEPPTPKTAWAETSRPPETEPGPPAPKPPLPPPHRGPAGNWGPPGDY']]
['jp_wtQxqUf', ['CUL4B', 699, 82.74897119, 'NVPGNIELTVNILTMGYWPTYVPMEVHLPPEMVKLQEIFKTFYLGKHSGRKLQWQSTLGH']]
['jp_1PhnURl', ['RPS10', 87, 72.76853612, 'VKEQFAWRHFYWYLTNEGIQYLRDYLHLPPEIVPATLRRSRPETGRPRPKGLEGERPARL']]
['jp_A188ZFw', ['ACOT2', 214, 70.58333333, 'TRHERYFLPPGVRREPVRVGRVRGTLFLPPEPGPFPGIVDMFGTGGGLLEYRASLLAGKG']]
['jp_qSbDgk5', ['NDUFB11', 49, 60.74074074, 'RGLPAARVRWESSFSRTVVAPSAVAGKRPPEPTTPWQEDPEPEDENLYEKNPDSHGYDKD']]
['jp_ahZhBbI', ['DPF2', 88, 59.05454545, 'KRHRGPGLASGQLYSYPARRWRKKRRAHPPEDPRLSFPSIKPDTDQTLKKEGLISQDGSS']]
['jp__T_RlmE', ['GNAS', 787, 58.59677419, 'PPVELANPENQFRVDYILSVMNVPDFDFPPEFYEHAKALWEDEGVRACYERSNEYQLIDC']]
['jp_PQibVAn', ['KCTD12', 263, 45.01904762, 'ARITVCGKTSLAKEVFGDTLNESRDPDRPPERYTSRYYLKFNFLEQAFDKLSESGFHMVA']]
['jp_8Ao4C4E', ['HNRNPUL2', 114, 44.75213675, 'DEEALLEDEDEEPPPAQALGQAAQPPPEPPEAAAMEAAAEPDASEKPAEATAGSGGVNGG']]
['jp_7H9lNQC', ['ADO', 149, 43.70967742, 'MLKVLYGTVRISCMDKLDAGGGQRPRALPPEQQFEPPLQPREREAVRPGVLRSRAEYTEA']]
['jp_4vQIxkx', ['CYB5A', 96, 40.5, 'GHSTDAREMSKTFIIGELHPDDRPKLNKPPETLITTIDSSSSWWTNWVIPAISAVAVALM']]
['jp_bteKOcC', ['HNRNPU', 679, 39.37546867, 'DEITYVELQKEEAQKLLEQYKEESKKALPPEKKQNTGSKKSNKNKSGKNQFNRGGGHRGR']]
['jp_r26rcVX', ['SF3B6', 15, 39.06976744, 'MAMQAAKRANIRLPPEVNRILYIRNLPYKITAEEMYDIFGKYGPI']]
['jp_R0GiaUW', ['AMOT', 285, 36.23987854, 'PQEPGHFYSEHRLNQPGRTEGQLMRYQHPPEYGAARPAQDISLPLSARNSQPHSPTSSLT']]
['jp_948ewVW', ['MRPS18B', 209, 36.07960199, 'SHGAVSATPPAPTLVSGDPWYPWYNWKQPPERELSRLRRLYQGHLQEESGPPPESMPKMP']]
['jp_8E7wzkZ', ['EIF3C', 730, 35.71743487, 'HESDARRRMISKQFHHQLRVGERQPLLGPPESMREHVVAASKAMKMGDWKTCHSFIINEK']]
['jp_XMaV56b', ['EIF3CL', 731, 35.15672396, 'HESDARRRMISKQFHHQLRVGERQPLLGPPESMREHVVAASKAMKMGDWKTCHSFIINEK']]
['jp_RqZ0Llv', ['ZNF706', 69, 31.92682927, 'TCTVCRTQMPDPKTFKQHFESKHPKTPLPPELADVQA\n']]
['jp_ylypXBR', ['RAD21', 506, 31.20987654, 'KRKAGQIDPEPVMPPQQVEQMEIPPVELPPEEPPNICQLIPELELLPEKEKEKEKEKEDD']]
['jp_P_o7f1Z', ['HSPBP1', 64, 31.03448276, 'GSGNSRPPRNLQGLLQMAITAGSEEPDPPPEPMSEERRQWLQEAMSAAFRGQREEVEQMK']]
['jp_Q3sAo1t', ['SLC25A3', 348, 28.39106145, 'IIMIGTLTALQWFIYDSVKVYFRLPRPPPPEMPESLKKKLGLTQ\n']]
['jp_FTzQOQw', ['PSPC1', 53, 27.81176471, 'SAVGESEPAAAAAMALALAGEPAPPAPAPPEDHPDEEMGFTIDIKSFLKPGEKTYTQRCR']]
['jp_gY7LrSs', ['UBE2V2', 138, 27.68965517, 'WQNSYSIKVVLQELRRLMMSKENMKLPQPPEGQTYNN\n']]
['jp_9mYEWSk', ['TAF7', 22, 25.2851711, 'MSKSKDDAPHELESQFILRLPPEYASTVRRAVQSGHVNLKDRLTIELHPDGR']]
['jp_FDlQzZF', ['STC2', 223, 24.73469388, 'QVQCEQNWGSLCSILSFCTSAIQKPPTAPPERQPQVDRTKLSRAHHGEAGHHLPEPSSRE']]
['jp_ltpb6rY', ['HADHA', 232, 24.17721519, 'TGRSIRADRAKKMGLVDQLVEPLGPGLKPPEERTIEYLEEVAITFAKGLADKKISPKRDK']]
['jp_BB5f4iD', ['HMG20A', 94, 151.2719486, 'SSNAAEGNEQRHEDEQRSKRGGWSKGRKRKKPLRDSNAPKSPLTGYVRFMNERREQLRAK']]
['jp_nQSgzIO', ['SSRP1', 533, 135.5454545, 'SSNEGDSDRDEKKRKQLKKAKMAKDRKSRKKPVEVKKGKDPNAPKRPMSAYMLWLNASRE']]
['jp_TOIpgmB', ['STK11', 311, 121.24, 'LSDLLKGMLEYEPAKRFSIRQIRQHSWFRKKHPPAEAPVPIPPSPDTKDRWRSMTVVPYL']]
['jp_ahMLnCI', ['PES1', 25, 121.080292, 'MGGLEKKKYERGSATNYITRNKARKKLQLSLADFRRLCILKGIYPHEPKHKKKVN']]
['jp_2hzwVm1', ['DDX18', 425, 113.3145695, 'TVDGLEQGYVVCPSEKRFLLLFTFLKKNRKKKLMVFFSSCMSVKYHYELLNYIDLPVLAI']]
['jp_XiAfUt1', ['IQCB1', 436, 97.66304348, 'FHQQRQSLIEYKAAVTLQRAALKFLAKCRKKKKLFAPWRGLQELTDARRVELKKRVDDYV']]
['jp_bFMS9aP', ['REEP4', 248, 95.46, 'LIRSQSLRVVKRKPPVREGTSRSLKVRTRKKTVPSDVDS\n']]
['jp_tpiwR3A', ['ITGB3BP', 65, 78.90909091, 'TGTCQMSLFASPTSSEEQKHRNGLSNEKRKKLNHPSLTESKESTTKDNDEFMMLLSKVEK']]
['jp_2nvvKYH', ['DDX56', 512, 65.9895288, 'AVVKPHLGHVPDYLVPPALRGLVRPHKKRKKLSSSCRKAKRAKSQNPLRSFKHKGKKFRP']]
['jp_RKKz8Ji', ['MRPS9', 380, 63.71604938, 'ALCSFVTEDEVEWMRQAGLLTTDPRVRERKKPGQEGARRKFTWKKR\n']]
['jp_6PvCkDu', ['ALKBH2', 242, 58.65671642, 'VVRLPLAHGSLLMMNHPTNTHWYHSLPVRKKVLAPRVNLTFRKILLTKK\n']]
['jp_ZFYGOl_', ['BCLAF1', 31, 56.5681382, 'GRSNSRSHSSRSKSRSQSSSRSRSRSHSRKKRYSSRSRSRTYSRSRSRDRMYSRDYRRDY']]
['jp_RwSJHlF', ['PPIG', 231, 55.76704545, 'SDSDSSSDSQSSSDSSDSESATEEKSKKRKKKHRKNSRKHKKEKKKRKKSKKSASSESEA']]
['jp_SLsnqAp', ['RAB35', 197, 47.27659574, 'VLRAKKDNLAKQQQQQQNDVVKLTKNSKRKKRCC\n']]
['jp_j6zns9L', ['SURF6', 346, 46.62121212, 'KMQQRQDRRRQNLRRKKAARAERRLLRARKKGRILPQDLERAGLV\n']]
['jp_AItV5rG', ['NKAP', 198, 39.24528302, 'KKSTTSASTSEEEKKKKSSRSKERSKKRRKKKSSKRKHKKYSEDSDSDSDSETDSSDEDN']]
['jp_L_XMxzS', ['PPIG', 249, 38.60795455, 'ESATEEKSKKRKKKHRKNSRKHKKEKKKRKKSKKSASSESEAENLEAQPQSTVRPEEIPP']]
['jp_eYeybEF', ['BTF3', 73, 33.55578947, 'QMKETIMNQEKLAKLQAQVRIGGKGTARRKKKVVHRTATADDKKLQFSLKKLGVNNISGI']]
['jp_hjqQ576', ['SURF6', 49, 32.90909091, 'SHSAPEQQARTRAGKTQGSETAGPPKKKRKKTQKKFRKREEKAAEHKAKSLGEKSPAASG']]
['jp_BZ4CZHk', ['PNRC1', 99, 32.43956044, 'PQPRAPAALPNRSLAVAGGTPRAAPKKRRKKKVRASPAGQLPSRFHQYQQHRPSLEGGRS']]
['jp_G1Yj_2R', ['RAP1B', 177, 32.24770642, 'AKSKINVNEIFYDLVRQINRKTPVPGKARKKSSCQLL\n']]
['jp_g4s0FWy', ['U2AF1', 191, 31.99115044, 'ECTRGGFCNFMHLKPISRELRRELYGRRRKKHRSRSRSRERRSRSRDRGRGGGGGGGGGG']]
['jp_veBQzed', ['U2AF1L5', 191, 31.99115044, 'ECTRGGFCNFMHLKPISRELRRELYGRRRKKHRSRSRSRERRSRSRDRGRGGGGGGGGGG']]
['jp_b13mIyP', ['BTF3L4', 24, 29.08536585, 'MNQEKLAKLQAQVRIGGKGTARRKKKVVHRTATADDKKLQSSLKKLAVNNIAGI']]
['jp_bfitrfq', ['NGDN', 299, 28.6, 'HSLTHFSDISALTGGTVHLDEDQNPIKKRKKIPQKGRKKKGFRRRR\n']]
['jp_8D8eSlW', ['RPL7', 35, 26.47594937, 'EKKKEVPAVPETLKKKRRNFAELKIKRLRKKFAQKMLRKARRKLIYEKAKHYHKEYRQMY']]
['jp_53AyQK5', ['EIF1AX', 67, 25.45801527, 'VIKMLGNGRLEAMCFDGVKRLCHIRGKLRKKVWINTSDIILVGLRDYQDNKADVILKYNA']]
['jp_mO0QCh0', ['STX12', 249, 25.18181818, 'EANVESSEVHVERATEQLQRAAYYQKKSRKKMCILVLVLSVIILILGLIIWLVYKTK\n']]
['jp_wMKWWE8', ['RPS27A', 81, 23.61347709, 'DGRTLSDYNIQKESTLHLVLRLRGGAKKRKKKSYTTPKKNKHKRKKVKLAVLKYYKVDEN']]
['jp_ZpXSv8g', ['UBA2', 426, 22.89285714, 'LIVLEGLKILSGKIDQCRTIFLNKQPNPRKKLLVPCALDPPNPNCYVCASKPEVTVRLNV']]
['jp_7sI0yl6', ['U2SURP', 390, 435.5420792, 'PIYIPPSMMEHTLPPPPSGLPFNAQPRERLKNPNAPMLPPPKNKEDFEKTLSQAIVKVVI']]
['jp_IQIR25d', ['SLC25A17', 129, 163.9647059, 'TTGKDLVVGFVAGVVNVLLTTPLWVVNTRLKLQGAKFRNEDIVPTNYKGIIDAFHQIIRD']]
['jp_N9G1Ep9', ['FAM98B', 107, 151.7482517, 'ESFQLEISGFLKEMACPYSVLISGDIKDRLKKKEDCLKLLLFLSTELQASQILQNKKHKN']]
['jp_2jHVjke', ['LHPP', 48, 115.5081967, 'SGVLYDSGAGGGTAIAGSVEAVARLKRSRLKVRFCTNESQKSRAELVGQLQRLGFDISEQ']]
['jp_5FL7amZ', ['TMCO1', 90, 111.1409396, 'LIVFISVCTALLAEGITWVLVYRTDKYKRLKAEVEKQSKKLEKKKETITESAGRQQKKKI']]
['jp_SWIFYjg', ['SPOP', 371, 88.44339623, 'VVSHPHLVAEAYRSLASAQCPFLGPPRKRLKQS\n']]
['jp_guVvsP1', ['C18orf21', 59, 81.67391304, 'TSSHDDKSTFEETCPYCFQLLVLDNSRVRLKPKARLTPKIQKLLNREARNYTLSFKEAKM']]
['jp_RBklIjk', ['CEBPG', 81, 67.31325301, 'SKKSSPMDRNSDEYRQRRERNNMAVKKSRLKSKQKAQDTLQRVNQLKEENERLEAKIKLL']]
['jp_iU2amkY', ['VBP1', 74, 60.26086957, 'TADTVLKKLDEQYQKYKFMELNLAQKKRRLKGQIPEIKQTLEILKYMQKKKESTNSMETR']]
['jp_VKSm89r', ['GNAI3', 91, 55.46875, 'SEDECKQYKVVVYSNTIQSIIAIIRAMGRLKIDFGEAARADDARQLFVLAGSAEEGVMTP']]
['jp_SFelBeX', ['TWF2', 97, 51.12359551, 'YLLYRLDSQNAQGFEWLFLAWSPDNSPVRLKMLYAATRATVKKEFGGGHIKDELFGTVKD']]
['jp_kOH9XNn', ['NSUN5', 102, 47.77777778, 'VLVYELLLGKGFRGGGGRWKALLGRHQARLKAELARLKVHRGVSRNEDLLEVGSRPGPAS']]
['jp_c9CdDEA', ['XRCC6', 259, 47.68961494, 'RVHFEESSKLEDLLRKVRAKETRKRALSRLKLKLNKDIVISVGIYNLVQKALKPPPIKLY']]
['jp_zTGWzcm', ['DAPK3', 302, 44.60784314, 'SWIKAIRRRNVRGEDSGRKPERRRLKTTRLKEYTIKSHSSLPPNNSYADFERFSKVLEEA']]
['jp_FX9aEEX', ['DGUOK', 119, 43.68571429, 'QSLGNLLDMMYREPARWSYTFQTFSFLSRLKVQLEPFPEKLLQARKPVQIFERSVYSDRY']]
['jp_vKh3UUb', ['ZNF593', 89, 40.75471698, 'LHRCLACARYFIDSTNLKTHFRSKDHKKRLKQLSVEPYSQEEAERAAGMGSYVPPRRLAV']]
['jp_nMW2QNl', ['RPS4X', 52, 38.57393337, 'LTGVFAPRPSTGPHKLRECLPLIIFLRNRLKYALTGDEVKKICMQRFIKIDGKVRTDITY']]
['jp_o_6eGkP', ['RPL10A', 55, 36.0097166, 'RKFLETVELQISLKNYDPQKDKRFSGTVRLKSTPRPKFSVCVLGDQQHCDEAKAVDIPHM']]
['jp_ZtZoVpd', ['SF3A3', 88, 32.38709677, 'DLYDDKDGLRKEELNAISGPNEFAEFYNRLKQIKEFHRKHPNEICVPMSVEFEELLKARE']]
['jp_OpbeRlB', ['NUDC', 193, 30.88372093, 'GNGADLPNYRWTQTLSELDLAVPFCVNFRLKGKDMVVDIQRRHLRVGLKGQPAIIDGELY']]
['jp_ul4ExYF', ['HADHA', 733, 29.01265823, 'AVFGLGFPPCLGGPFRFVDLYGAQKIVDRLKKYEAAYGKQFTPCQLLADHANSPNKKFYQ']]
['jp_5gG7Qo2', ['EIF3E', 192, 28.40764331, 'RNALSSLWGKLASEILMQNWDAAMEDLTRLKETIDNNSVSSPLQSLQQRTWLIHWSLFVF']]
['jp_8pcyqj7', ['DDX5', 196, 26.94587629, 'GPICLVLAPTRELAQQVQQVAAEYCRACRLKSTCIYGGAPKGPQIRDLERGVEICIATPG']]
['jp_bAakde9', ['IQCB1', 314, 26.04347826, 'VYQEVEEQKLHQAACLIQAYWKGFQTRKRLKKLPSAVIALQRSFRSKRSKMLLEINRQKE']]
['jp_mYNVFhe', ['RBM5', 199, 25.2046332, 'HIAMHYSNPRPKFEDWLCNKCCLNNFRKRLKCFRCGADKFDSEQEVPPGTTESVQSVDYY']]
['jp_0ecHWkp', ['GADD45GIP1', 190, 24.99137931, 'QELLGYQVDPRSARFQELLQDLEKKERKRLKEEKQKRKKEARAAALAAAVAQDPAASGAP']]
['jp_JcYLoh9', ['LLPH', 105, 23.63636364, 'ETDIKRNKKTLLDQHGQYPIWMNQRQRKRLKAKREKRKGKSKAKAVKVAKGLAW\n']]
['jp_eXbMsXh', ['PRC1', 341, 23.22637795, 'QEQRQAFAPFCAEDYTESLLQLHDAEIVRLKNYYEVHKELFEGVQKWEETWRLFLEFERK']]
['jp_2dEPCta', ['BRD2', 431, 199.006993, 'KMENRDYRDAQEFAADVRLMFSNCYKYNPPDHDVVAMARKLQDVFEFRYAKMPDEPLEPG']]
['jp_FDeHscZ', ['SMARCE1', 66, 179.0319032, 'NNYRLGGNPGTNSRVTASSGITIPKPPKPPDKPLMPYMRYSRKVWDQVKASNPDLKLWEI']]
['jp_uWOvDgF', ['PFAS', 880, 136.9933993, 'LLYVALSPGQHRLGGTALAQCFSQLGEHPPDLDLPENLVRAFSITQGLLKDRLLCSGHDV']]
['jp_Ur7vHZP', ['BPHL', 105, 100.238806, 'DFGPQLKNLNKKLFTVVAWDPRGYGHSRPPDRDFPADFFERDAKDAVDLMKALKFKKVSL']]
['jp_BWk0TBC', ['GRPEL2', 50, 87.31818182, 'AAWESKGWPLPFSTATQRTAGEDCRSEDPPDELGPPLAERALRVKAVKLEKEVQDLTVRY']]
['jp_CY1jHhf', ['TXLNG', 119, 62.72727273, 'TENRNLVSPAYCTQESREEIPGGEARTDPPDGQQDSECNRNKEKTLGKEVLLLMQALNTL']]
['jp_HPtgb4B', ['NONO', 83, 41.54929577, 'GLTIDLKNFRKPGEKTFTQRSRLFVGNLPPDITEEEMRKLFEKYGKAGEVFIHKDKGFGF']]
['jp_6Wl4vnL', ['UBC', 190, 40.72705602, 'TGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLH']]
['jp_f24e3Ao', ['HOXA10', 130, 36.91616766, 'GGLGPGAHGYGPSPIDLWLDAPRSCRMEPPDGPPPPPQQQPPPPPQPPQPAPQATSCSFA']]
['jp_Ildqaub', ['PPIL4', 170, 34.16831683, 'ETFVDKDFVPYQDIRINHTVILDDPFDDPPDLLIPDRSPEPTREQLDSGRIGADEEIDDF']]
['jp_aaUb6dA', ['ABCE1', 138, 34.04255319, 'VGTNGIGKSTALKILAGKQKPNLGKYDDPPDWQEILTYFRGSELQNYFTKILEDDLKAII']]
['jp_s6FZwEh', ['EEF2', 377, 28.5156497, 'DALLQMITIHLPSPVTAQKYRCELLYEGPPDDEAAMGIKSCDPKGPLMMYISKMVPTSDK']]
['jp_ZELRlOm', ['TAF9', 117, 28.39285714, 'PRDFLLDIARQRNQTPLPLIKPYSGPRLPPDRYCLTAPNYRLKSLQKKASTSAGRITVPR']]
['jp_2kAYcaL', ['SCNM1', 225, 28.28571429, 'LRSSGWIPDGRGRWVKDENVEFDSDEEEPPDLPLD\n']]
['jp_2zjYyry', ['KIAA1143', 43, 28.18181818, 'EPAFLARFKERVGYREGPTVETKRIQPQPPDEDGDHSDKEDEQPQVVVLKKGDLSVEEVM']]
['jp_4wXbW2t', ['UBE2S', 28, 27.73333333, 'MNSNVENLPPHIIRLVYKEVTTLTADPPDGIKVFPNEEDLTDLQVTIEGPEGTPYAGG']]
['jp_Vj43MWn', ['NASP', 720, 77.63615561, 'VSMIASRKPTDGASSSNCVTDISHLVRKKRKPEEESPRKDDAKKAKQEPEVNGGSGDAVP']]
['jp_g2pxZY7', ['LMNB1', 416, 68.10497238, 'RLKLSPSPSSRVTVSRASSSRSVRTTRGKRKRVDVEESEASSSVSISHSASATGNVCIEE']]
['jp_2pL90tI', ['SUPT16H', 1046, 61.82890855, 'SRKRKASVHSSGRGSNRGSRHSSAPPKKKRK\n']]
['jp_UR60OxY', ['EIF4A3', 288, 58.05454545, 'EREEWKFDTLCDLYDTLTITQAVIFCNTKRKVDWLTEKMREANFTVSSMHGDMPQKERES']]
['jp_OIfQWI4', ['EMC3', 93, 40.5915493, 'NGKYIPKQSFLTRKYYFNNPEDGFFKKTKRKVVPPSPMTDPTMLTDMMKGNVTNVLPMIL']]
['jp_1uvbH81', ['DDX56', 511, 37.29842932, 'PAVVKPHLGHVPDYLVPPALRGLVRPHKKRKKLSSSCRKAKRAKSQNPLRSFKHKGKKFR']]
['jp_bFbv_ko', ['DDX18', 115, 31.10596026, 'SPQKSTVLTNGEAAMQSSNSESKKKKKKKRKMVNDAEPDTKKAKTENKGKSEEESAETTK']]
['jp_lnf2dJs', ['ZC3H15', 261, 30.0234375, 'LEDLIERERSALGPNVTKITLESFLAWKKRKRQEKIDKLEQDMERRKADFKAGKALVISG']]
['jp_TqFB2fZ', ['CDC5L', 169, 29.27604167, 'DPIDMDEDELEMLSEARARLANTQGKKAKRKAREKQLEEARRLAALQKRRELRAAGIEIQ']]
['jp_k_ezrWN', ['CYBA', 59, 27.22580645, 'AGRFTQWYFGAYSIVAGVFVCLLEYPRGKRKKGSTMERWGQKYMTAVVKLFGPFTRNYYV']]
['jp_TYcIzFR', ['THUMPD3', 146, 26.73684211, 'FEDLAGKLPWSNPLKVWKINASFKKKKAKRKKINQNSSKEKINNGQEVKIDQRNVKKEFT']]
['jp_hctL6F0', ['RPS8', 25, 26.27027027, 'MGISRDNWHKRRKTGGKRKPYHKKRKYELGRPAANTKIGPRRIHTVRVRGGNKKY']]
['jp_1icT921', ['SNRPA', 111, 23.90463215, 'MRIQYAKTDSDIIAKMKGTFVERDRKREKRKPKSQETPATKKAVQGGGATPVVGAVQGPV']]
['jp_GAmie9R', ['LMNB2', 436, 23.70229008, 'SPSSRVTVSRATSSSSGSLSATGRLGRSKRKRLEVEEPLGSGPSVLGTGTGGSGGFHLAQ']]
['jp_V2N6PGV', ['VAT1', 104, 57.99134199, 'GPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRV']]
['jp_700wrWS', ['PREPL', 387, 55.34502924, 'HCVLFLKHSNLLYVNVIGLADDSVRSLKLPPWACGFIMDTNSDPKNCPFQLCSPIRPPKY']]
['jp_rpiKlMQ', ['LRPAP1', 72, 47.73333333, 'KPSPKRESGEEFRMEKLNQLWEKAQRLHLPPVRLAELHADLKIQERDELAWKKLKLDGLD']]
['jp_I1vyflc', ['XRCC5', 506, 41.32107843, 'FPTTKIPNPRFQRLFQCLLHRALHPREPLPPIQQHIWNMLNPPAEVTTKSQIPLSKIKTL']]
['jp_Td3GeFj', ['CUL4B', 698, 37.61316872, 'QNVPGNIELTVNILTMGYWPTYVPMEVHLPPEMVKLQEIFKTFYLGKHSGRKLQWQSTLG']]
['jp_ZfphclS', ['SON', 1150, 35.36612022, 'DRSMMSMAADSYTDSYTDTYTEAYMVPPLPPEEPPTMPPLPPEEPPMTPPLPPEEPPEGP']]
['jp_o4UoJPR', ['PRC1', 519, 34.22834646, 'TTTMSNATANSSIRPIFGGTVYHSPVSRLPPSGSKPVAASTCSGKKTPRTGRHGANKENL']]
['jp_wa_sV0K', ['SPCS1', 132, 31.01351351, 'FIYGYVAEQFGWTVYIVMAGFAFSCLLTLPPWPIYRRHPLKWLPVQESSTDDKKPGERKI']]
['jp_BI6fxex', ['CKMT1A', 176, 30.4, 'SKIRSGYFDERYVLSSRVRTGRSIRGLSLPPACTRAERREVERVVVDALSGLKGDLAGRY']]
['jp_K_aksMY', ['HADHA', 662, 29.01265823, 'GKGFYIYQEGVKRKDLNSDMDSILASLKLPPKSEVSSDEDIQFRLVTRFVNEAVMCLQEG']]
['jp_SZ5Mbke', ['HNRNPU', 678, 28.75736476, 'FDEITYVELQKEEAQKLLEQYKEESKKALPPEKKQNTGSKKSNKNKSGKNQFNRGGGHRG']]
['jp_IszKDDw', ['RPS10', 86, 23.84505703, 'YVKEQFAWRHFYWYLTNEGIQYLRDYLHLPPEIVPATLRRSRPETGRPRPKGLEGERPAR']]
['jp_jWKJrhs', ['KDM4A', 837, 210.0519031, 'VNIAERSPVDVSKIPLPRFKLKCIFCKKRRKRTAGCCVQCSHGRCPTAFHVSCAQAAGVM']]
['jp_uCo7u_1', ['RACGAP1', 247, 64.79527559, 'KTTVTVPNDGGPIEAVSTIETVPYWTRSRRKTGTLQPWNSDSTLNSRQLEPRTETDSVGT']]
['jp_Ih2mncb', ['RPL10A', 26, 56.30931174, 'MSSKVSRDTLYEAVREVLHGNQRKRRKFLETVELQISLKNYDPQKDKRFSGTVRLK']]
['jp_Xj4xt0p', ['SOX4', 75, 50.16163793, 'ADDPSWCKTPSGHIKRPMNAFMVWSQIERRKIMEQSPDMHNAEISKRLGKRWKLLKDSDK']]
['jp_sRD8cju', ['SIX1', 113, 49.56521739, 'QQLWLKAHYVEAEKLRGRPLGAVGKYRVRRKFPLPRTIWDGEETSYCFKEKSRGVLREWY']]
['jp_KuaWacL', ['PNRC1', 98, 46.85714286, 'TPQPRAPAALPNRSLAVAGGTPRAAPKKRRKKKVRASPAGQLPSRFHQYQQHRPSLEGGR']]
['jp_AqZeW94', ['TMEM101', 46, 45.0, 'QLGSVLLTRCPFWGCFSQLMLYAERAEARRKPDIPVPYLYFDMGAAVLCASFMSFGVKRR']]
['jp_N9TwwCH', ['BSG', 349, 39.86557377, 'SHLAALWPFLGIVAEVLVLVTIIFIYEKRRKPEDVLDDDDAGSAPLKSSGQHQNDKGKNV']]
['jp_QxiPNeS', ['RPL23', 90, 30.96031746, 'VMATVKKGKPELRKKVHPAVVIRQRKSYRRKDGVFLYFEDNAGVIVNNKGEMKGSAITGP']]
['jp_ZkF1JTN', ['NKAP', 197, 27.47169811, 'PKKSTTSASTSEEEKKKKSSRSKERSKKRRKKKSSKRKHKKYSEDSDSDSDSETDSSDED']]
['jp_o1rp3ik', ['ABRACL', 56, 24.84848485, 'SVKFGVLFRDDKCANLFEALVGTLKAAKRRKIVTYPGELLLQGVHDDVDIILLQD\n']]
['jp_GxOPmMC', ['C11orf98', 70, 23.27027027, 'IDQGLITRHHLKKRASSARANITLSGKKRRKLLQQIRLAQKEKTAMEVEAPSKPARTSEP']]
['jp_blzFH34', ['ATXN10', 296, 104.0437158, 'VFLRHAELIASTFVDQCKTVLKLASEEPPDDEEALATIRLLDVLCEMTVNTELLGYLQVF']]
['jp_nf3LXmg', ['AAMP', 41, 91.57894737, 'DTPPLETLSFHGDEEIIEVVELDPGPPDPDDLAQEMEDVDFEEEEEEEGNEEGWVLEPQE']]
['jp_20zI2Pj', ['ITM2B', 84, 64.94594595, 'WCMCFGLAFMLAGVILGGAYLYKYFALQPDDVYYCGIKYIKDDVILNEPSADAPAALYQT']]
['jp_bck3sBF', ['DPF2', 161, 59.05454545, 'GAPDPRVDDDSLGEFPVTNSRARKRILEPDDFLDDLDDEDYEEDTPKRRGKGKSKGKGVG']]
['jp_Ms6T714', ['GTF2I', 853, 57.58646617, 'PPRKINSSPNVNTTASGVEDLNIIQVTIPDDDNERLSKVEKARQLREQVNDLFSRKFGEA']]
['jp_RXp1kGo', ['PPP2CA', 204, 26.57142857, 'DTLDHIRALDRLQEVPHEGPMCDLLWSDPDDRGGWGISPRGAGYTFGQDISETFNHANGL']]
['jp_32swZxX', ['SORD', 125, 24.97674419, 'PGAPRENDEFCKMGRYNLSPSIFFCATPPDDGNLCRFYKHNAAFCYKLPDNVTFEEGALI']]
['jp_lB133ii', ['AHSA1', 40, 130.516129, 'WIVEERADATNVNNWHWTERDASNWSTDKLKTLFLAVQVQNEEGKCEVTEVSKLDGEASI']]
['jp_AJpoqH9', ['DDX39A', 268, 94.82634731, 'RKFMQDPMEVFVDDETKLTLHGLQQYYVKLKDSEKNRKLFDLLDVLEFNQVIIFVKSVQR']]
['jp_7CCJBaT', ['SNRPF', 23, 73.66906475, 'MSLPLNPKPFLNGLTGKPVMVKLKWGMEYKGYLVSVDGYMNMQLANTEEYIDG']]
['jp_qylIoM2', ['RALY', 182, 70.31843575, 'PLVRRVKTNVPVKLFARSTAVTTSSAKIKLKSSELQAIKTELTQIKSNIDALLSRLEQIA']]
['jp_a0SUgOF', ['RPS19', 28, 69.00252738, 'MPGVTVKDVNQQEFVRALAAFLKKSGKLKVPEWVDTVKLAKHKELAPYDENWFYTRAA']]
['jp_R0Qt3Sy', ['FASN', 54, 53.33333333, 'NLIGGVDMVTDDDRRWKAGLYGLPRRSGKLKDLSRFDASFFGVHPKQAHTMDPQLRLLLE']]
['jp_50Md3uK', ['POLR3K', 46, 36.75581395, 'GQRCHRFSCNTCPYVHNITRKVTNRKYPKLKEVDDVLGGAAAWENVDSTAESCPKCEHPR']]
['jp_hYNGHKD', ['BTG3', 22, 31.93548387, 'MKNEIAAVVFFFTRLVRKHDKLKKEAVERFAEKLTLILQEKYKNHWYPEKPS']]
['jp_9Qw6Hcl', ['SNRPA1', 152, 29.86666667, 'NPVTNKKHYRLYVIYKVPQVRVLDFQKVKLKERQEAEKMFKGKRGAQLAKDIARRSKTFN']]
['jp_rxwprqq', ['RAE1', 308, 28.9245283, 'IAFHPVHGTLATVGSDGRFSFWDKDARTKLKTSEQLDQPISACCFNHNGNIFAYASSYDW']]
['jp_rb0PWep', ['LBR', 177, 27.17647059, 'NTQEKFSLSQESSYIATQYSLRPRREEVKLKEIDSKEEKYVAKELAVRTFEVTPIRAKDL']]
['jp_6itnheU', ['DDX39B', 269, 26.11658291, 'RKFMQDPMEIFVDDETKLTLHGLQQYYVKLKDNEKNRKLFDLLDVLEFNQVVIFVKSVQR']]
['jp_gz2GBfq', ['CCND2', 147, 25.95330739, 'LTAEKLCIYTDNSIKPQELLEWELVVLGKLKWNLAAVTPHDFIEHILRKLPQQREKLSLI']]
['jp_zDK1uA_', ['AHCYL1', 430, 197.4793388, 'NTEIDVTSLRTPELTWERVRSQVDHVIWPDGKRVVLLAEGRLLNLSCSTVPTFVLSITAT']]
['jp_lOJiedH', ['TXLNG', 120, 108.7272727, 'ENRNLVSPAYCTQESREEIPGGEARTDPPDGQQDSECNRNKEKTLGKEVLLLMQALNTLS']]
['jp_JcglCpL', ['TMED4', 85, 102.2068966, 'YRTQMWDKQKEVFLPSTPGLGMHVEVKDPDGKVVLSRQYGSEGRFTFTSHTPGDHQICLH']]
['jp__nGVsTg', ['TMED7', 81, 73.7704918, 'DIAQGTKCTLEFQVITGGHYDVDCRLEDPDGKVLYKEMKKQYDSFTFTASKNGTYKFCFS']]
['jp_6QGozvx', ['PAWR', 131, 42.09876543, 'APGPRRSEDEPPAASASAAPPPQRDEEEPDGVPEKGKSSGPSARKGKGQIEKRKLREKRR']]
['jp_ypJyz8U', ['GLUD1', 341, 39.0, 'VGLHSMRYLHRFGAKCIAVGESDGSIWNPDGIDPKELEDFKLQHGSILGFPKAKPYEGSI']]
['jp_doci5Gc', ['SHMT2', 168, 33.05799649, 'QPYSGSPANLAVYTALLQPHDRIMGLDLPDGGHLTHGYMSDVKRISATSIFFESMPYKLN']]
['jp_whnRPwt', ['TCOF1', 1242, 23.78594249, 'PANSQASKATPKLDSSPSVSSTLAAKDDPDGKQEAKPQQAAGMLSPKTGGKEAASGTTPQ']]
['jp_02Tzhm3', ['ACBD3', 61, 136.6583333, 'APLLPPPLPPPSPPGSGRGPGASGEQPEPGEAAAGGAAEEARRLEQRWGFGLEELYGLAL']]
['jp_C5HUUsK', ['TOMM40', 242, 119.1666667, 'GSGILVAHYLQSITPCLALGGELVYHRRPGEEGTVMSLAGKYTLNNWLATVTLGQAGMHA']]
['jp_oClHUoP', ['GNAS', 972, 92.08064516, 'KVLAGKSKIEDYFPEFARYTTPEDATPEPGEDPRVTRAKYFIRDEFLRISTASGDGRHYC']]
['jp_Xsphh8a', ['TACC3', 270, 50.70879121, 'VCAPAAVATSPPGAIPKEACGGAPLQGLPGEALGCPAGVGTPVPADGTQTLTCAHTSAPE']]
['jp_Zcy1_PJ', ['RAD23A', 289, 50.05, 'QQLGQENPQLLQQISRHQEQFIQMLNEPPGELADISDVEGEVGAIGEEAPQMNYIQVTPQ']]
['jp_aCIVZpx', ['ABCE1', 103, 34.04255319, 'NLPSNLEKETTHRYCANAFKLHRLPIPRPGEVLGLVGTNGIGKSTALKILAGKQKPNLGK']]
['jp_xYn93UA', ['RFXANK', 28, 30.70588235, 'MELTQPAEDLIQTQQTPASELGDPEDPGEEAADGSDTVVLSLFPCTPEPVNPEPDASV']]
['jp_iv8iEni', ['NONO', 66, 26.5915493, 'PPPPIPANGQQASSQNEGLTIDLKNFRKPGEKTFTQRSRLFVGNLPPDITEEEMRKLFEK']]
['jp_JCy3420', ['VRK1', 297, 23.97986577, 'PKYVRDSKIRYRENIASLMDKCFPEKNKPGEIAKYMETVKLLDYTEKPLYENLRDILLQG']]
['jp__4bHKX0', ['PSMB6', 108, 114.5945946, 'RSGSAADTQAVADAVTYQLGFHSIELNEPPLVHTAASLFKEMCYRYREDLMAGIIIAGWD']]
['jp_snzoMmr', ['VAT1', 105, 86.98701299, 'PGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIAVGEGVSDRKAGDRVM']]
['jp_mnTL3nQ', ['PUF60', 58, 74.66666667, 'AGDKWKPPQGTDSIKMENGQSTAAKLGLPPLTPEQQEALQKAKKYAMEQSIKSVLVKQTI']]
['jp_0dNr08i', ['MRPS23', 40, 46.08433735, 'SIFSRTRDLVRAGVLKEKPLWFDVYDAFPPLREPVFQRPRVRYGKAKAPIQDIWYHEDRI']]
['jp_8gHAIDr', ['NDUFB9', 165, 40.64516129, 'EVKQLQEETPPGGPLTEALPPARKEGDLPPLWWYIVTRPRERPM\n']]
['jp_qPAl4un', ['SON', 1148, 39.78688525, 'TADRSMMSMAADSYTDSYTDTYTEAYMVPPLPPEEPPTMPPLPPEEPPMTPPLPPEEPPE']]
['jp_BcJv6Rb', ['PARP1', 750, 29.13875598, 'SQGSSDSQILDLSNRFYTLIPHDFGMKKPPLLNNADSVQAKVEMLDNLLDIEVAYSLLRG']]
['jp_HSpRoLv', ['SON', 1159, 26.52459016, 'DSYTDSYTDTYTEAYMVPPLPPEEPPTMPPLPPEEPPMTPPLPPEEPPEGPALPTEQSAL']]
['jp_wSzms1R', ['MAD2L2', 106, 26.47058824, 'KNDVEKVVVVILDKEHRPVEKFVFEITQPPLLSISSDSLLSHVEQLLRAFILKISVCDAV']]
['jp_YRKdXLa', ['DHFR', 27, 26.4375, 'MVGSLNCIVAVSQNMGIGKNGDLPWPPLRNEFRYFQRMTTTSSVEGKQNLVIMGKKT']]
['jp__Ayqpm7', ['TACC3', 366, 23.04945055, 'KLEFDVSDGATSKRAPPPRRLGERSGLKPPLRKAAVRQQKAPQEVEEDDGRSGAGEDPPM']]
['jp_yi2VoFw', ['ILF2', 49, 78.71617162, 'GGGFRPFVPHIPFDFYLCEMAFPRVKPAPDETSFSEALLKRNQDLAPNSAEQASILSLVT']]
['jp_kwpFOrr', ['IPO7', 597, 73.90376569, 'EYSEEVTPIAVEMTQHLAMTFNQVIQTGPDEEGSDDKAVTAMGILNTIDTLLSVVEDHKE']]
['jp_KlLti4c', ['GSTA4', 209, 38.93650794, 'FLQEYTVKLSNIPTIKRFLEPGSKKKPPPDEIYVRTVYNIFRP\n']]
['jp_LgskHBe', ['CDC34', 200, 38.33823529, 'GTKVDAERDGVKVPTTLAEYCVKTKAPAPDEGSDLFYDDYYEDGEVEEEADSCFGDDEDD']]
['jp_GRfhN2h', ['MRPL13', 130, 31.37113402, 'VAIVKLAIYGMLPKNLHRRTMMERLHLFPDEYIPEDILKNLVEELPQPRKIPKRLDEYTQ']]
['jp_UiaLKgD', ['GRPEL2', 51, 28.25, 'AWESKGWPLPFSTATQRTAGEDCRSEDPPDELGPPLAERALRVKAVKLEKEVQDLTVRYQ']]
['jp_vKvrRy0', ['PREPL', 166, 25.54385965, 'ADNDNYEVLFNLEELKLDQPFIDCIRVAPDEKYVAAKIRTEDSEASTCVIIKLSDQPVME']]
['jp_0SF9b84', ['HSP90B1', 752, 167.8383838, 'LRSGYLLPDTKAYGDRIERMLRLSLNIDPDAKVEEEPEEEPEETAEDTTEDTEQDEDEEM']]
['jp_JNyVdOU', ['THRAP3', 225, 102.819222, 'NQGDEAKEQTFSGGTSQDTKASESSKPWPDATYGTGSASRASAVSELSPRERSPALKSPL']]
['jp_pYQ7Dg7', ['VDAC1', 230, 75.20160481, 'KLETAVNLAWTAGNSNTRFGIAAKYQIDPDACFSAKVNNSSLIGLGYTQTLKPGIKLTLS']]
['jp_1NBlcyH', ['CALR', 212, 68.52459016, 'EVKIDNSQVESGSLEDDWDFLPPKKIKDPDASKPEDWDERAKIDDPTDSKPEDWDKPEHI']]
['jp_8g4vhJX', ['PEBP1', 72, 58.35756385, 'PTQVKNRPTSISWDGLDSGKLYTLVLTDPDAPSRKDPKYREWHHFLVVNMKGNDISSGTV']]
['jp_TOSFAI1', ['ANXA5', 164, 41.68831169, 'SSLEDDVVGDTSGYYQRMLVVLLQANRDPDAGIDEAQVEQDAQALFQAGELKWGTDEEKF']]
['jp_fYqLCFY', ['MAT2A', 51, 33.55932203, 'SESVGEGHPDKICDQISDAVLDAHLQQDPDAKVACETVAKTGMILLAGEITSRAAVDYQK']]
['jp_yow6q0V', ['XPO1', 58, 32.62608696, 'DNVVNCLYHGEGAQQRMAQEVLTHLKEHPDAWTRVDTILEFSQNMNTKYYGLQILENVIK']]
['jp_WOPC7qU', ['PPA1', 165, 28.43925234, 'IGVKVLGILAMIDEGETDWKVIAINVDDPDAANYNDINDVKRLKPGYLEATVDWFRRYKV']]
['jp_jtQByxc', ['NDUFA9', 257, 22.5, 'LGSLGWKTVKQPVYVVDVSKGIVNAVKDPDANGKSFAFVGPSRYLLFHLVKYIFAVAHRL']]
['jp_Ob5fzt2', ['PRPF8', 1299, 117.1550152, 'LMTYFREAVVNTQELLDLLVKCENKIQTRIKIGLNSKMPSRFPPVVFYTPKELGGLGMLS']]
['jp_Ega37dZ', ['DRAP1', 17, 109.097561, 'MPSKKKKYNARFPPARIKKIMQTDEEIGKVAAAVPVIISRALELFLE']]
['jp_S6ygXXI', ['DDX21', 319, 101.8181818, 'ACFYGGTPYGGQFERMRNGIDILVGTPGRIKDHIQNGKLDLTKLKHVVLDEVDQMLDMGF']]
['jp_YrwY6Ro', ['MRPL47', 219, 36.73170732, 'LNKRYNRKRFFALPYVDHFLRLEREKRARIKARKENLERKKAKILLKKFPHLAEAQKSSL']]
['jp_jEBb1bB', ['TOP2A', 435, 34.63652174, 'GIVESILNWVKFKAQVQLNKKCSAVKHNRIKGIPKLDDANDAGGRNSTECTLILTEGDSA']]
['jp_VobOTo1', ['MRPL30', 72, 31.29545455, 'EKVFQASPEDHEKYGGDPQNPHKLHIVTRIKSTRRRPYWEKDIIKMLGLEKAHTPQVHKN']]
['jp_Dvhsy2T', ['MRPL23', 87, 28.0, 'EGIYNVPVAAVRTRVQHGSNKRRDHRNVRIKKPDYKVAYVQLAHGQTFTFPDLFPEKDES']]
['jp_AR2_OiR', ['BAG6', 743, 55.8739726, 'GSLGARAGSSESIAAFIQRLSGSSNIFEPGADGALGFFGALLSLLCQNFSMVDVVMLLHG']]
['jp__txLXnP', ['MRPL27', 78, 54.18181818, 'RRQGIKKMEGHYVHAGNIIATQRHFRWHPGAHVGVGKNKCLYALEEGIVRYTKEVYVPHP']]
['jp_H67ruT0', ['PSMA5', 131, 46.67142857, 'TYNETMTVESVTQAVSNLALQFGEEDADPGAMSRPFGVALLFGGVDEKGPQLFHMDPSGT']]
['jp_r1r9Uip', ['MRPL36', 45, 46.0, 'LSRHTVKPRALSTFLFGSIRGAAPVAVEPGAAVRSLLSPGLLPHLLPALGFKNKTVLKKR']]
['jp_XUTIt0k', ['SAE1', 77, 34.09090909, 'GAEIAKNLILAGVKGLTMLDHEQVTPEDPGAQFLIRTGSVGRNRAEASLERAQNLNPMVD']]
['jp_Cdm4OrS', ['SORD', 97, 24.97674419, 'LGHEASGTVEKVGSSVKHLKPGDRVAIEPGAPRENDEFCKMGRYNLSPSIFFCATPPDDG']]
['jp_bS29278', ['NRGN', 25, 24.6875, 'MDCCTENACSKPDDDILDIPLDDPGANAAAAKIQASFRGHMARKKIKSGERGRKG']]
['jp_d8Z8nR9', ['RPL32', 31, 60.37592745, 'AALRPLVKPKIVKKRTKKFIRHQSDRYVKIKRNWRKPRGIDNRVRRRFKGQILMPNIGYG']]
['jp_GdUc1hy', ['SSRP1', 565, 60.24242424, 'VEVKKGKDPNAPKRPMSAYMLWLNASREKIKSDHPGISITDLSKKAGEIWKGMSKEKKEE']]
['jp_TVFtA1Z', ['RPL27', 68, 35.28129952, 'SHALVAGIDRYPRKVTAAMGKKKIAKRSKIKSFVKVYNYNHLMPTRYSVDIPLDKTVVNK']]
['jp_EuTLvU5', ['VDAC1', 114, 33.89769308, 'ITVEDQLARGLKLTFDSSFSPNTGKKNAKIKTGYKREHINLGCDMDFDIAGPSIRGALVL']]
['jp_rnlW4iI', ['GTF2I', 801, 22.53383459, 'EGLPEGVPFRRPSTFGIPRLEKILRNKAKIKFIIKKPEMFETAIKESTSSKSPPRKINSS']]
['jp_XkED5vk', ['CCDC28B', 53, 64.76056338, 'LRRVPVPTSHSGSLALGLPHLPSPKQRAKFKRVGKEKCRPVLAGGGSGSAGTPLQHSFLT']]
['jp_qWrptKL', ['PHLDA3', 105, 52.48920863, 'VTEGGGEIDFRCPLEDPGWNAQITLGLVKFKNQQAIQTVRARQSLGTGTLVS\n']]
['jp_prUAaYf', ['RPL10', 157, 47.8960396, 'RVHIGQVIMSIRTKLQNKEHVIEALRRAKFKFPGRQKIHISKKWGFTKFNADEFEDMVAE']]
['jp_wwM5msM', ['DLX2', 207, 46.36912752, 'ALPERAELAASLGLTQTQVKIWFQNRRSKFKKMWKSGEIPSEQHPGASASPPCASPPVSA']]
['jp_jGrDrLl', ['TSPYL1', 321, 44.33742331, 'MIRGQDAEMLRYITNLEVKELRHPRTGCKFKFFFRRNPYFRNKLIVKEYEVRSSGRVVSL']]
['jp_iUlsE0M', ['KDELR1', 82, 39.17875648, 'TNYISLYNTCMKVVYIACSFTTVWLIYSKFKATYDGNHDTFRVEFLVVPTAILAFLVNHD']]
['jp_yQ7D4hn', ['IPO7', 488, 26.08368201, 'NHVFPLFSSELGYMRARACWVLHYFCEVKFKSDQNLQTALELTRRCLIDDREMPVKVEAA']]
['jp_ZeQLPFo', ['NACA', 2006, 24.3442623, 'SPASDTYIVFGEAKIEDLSQQAQLAAAEKFKVQGEAVSNIQENTQTPTVQEESEEEEVDE']]
['jp_TejxOW9', ['TIPIN', 90, 23.05343511, 'KRNIPKLDAQRLISERGLPALRHVFDKAKFKGKGHEAEDLKMLIRHMEHWAHRLFPKLQF']]
['jp_xKJbC68', ['IRS4', 369, 312.5197481, 'CRSYSISIGAHLLTLLSARRHLGLVPLEPGGWLRRSRFEQFCHLRAIGDGEDEMLFTRRF']]
['jp_9eVk6km', ['DXO', 358, 287.9340659, 'NFCAAFLSFAQSTVVQDDPRLVHLFSWEPGGPVTVSVHQDAPYAFLPIWYVEAMTQDLPS']]
['jp_mghnsj5', ['GRN', 35, 240.2696629, 'SWVALTAGLVAGTRCPDGQFCPVACCLDPGGASYSCCRPLLDKWPTTLSRHLGGPCQVDA']]
['jp_z0XhXR0', ['SMARCA4', 278, 115.2983871, 'YSRPHGMGGPNMPPPGPSGVPPGMPGQPPGGPPKPWPEGPMANAAAPTSTPQKLIPPQPT']]
['jp_egSZdZO', ['RNF149', 102, 76.14503817, 'AHGLVGVPWAPGGDLEGCAPDTRFFVPEPGGRGAAPWVALVARGGCTFKDKVLVAARRNA']]
['jp_9zuxKU9', ['ZNF503', 189, 36.90114068, 'GPLKLSDIGVEDKSSFKPYSKPGSDKKEPGGGGGGGGGGGGGGGGVSSEKSGFRVPSATC']]
['jp_09PVmDr', ['HSPA1L', 272, 309.6, 'FVEEFKRKHKKDISQNKRAVRRLRTACERAKRTLSSSTQANLEIDSLYEGIDFYTSITRA']]
['jp_JxFokpO', ['NACA', 2056, 102.2459016, 'EESEEEEVDETGVEVKDIELVMSQANVSRAKAVRALKNNSNDIVNAIMELTM\n']]
['jp_lEwVrbE', ['RDX', 82, 88.73333333, 'KGYSTWLKLNKKVTQQDVKKENPLQFKFRAKFFPEDVSEELIQEITQRLFFLQVKEAILN']]
['jp_g2L3hyg', ['EZR', 82, 38.61842105, 'KGFPTWLKLDKKVSAQEVRKENPLQFKFRAKFYPEDVAEELIQDITQKLFFLQVKEGILS']]
['jp_O0hwwJB', ['ARPC2', 268, 34.3245614, 'NTINLIHTFRDYLHYHIKCSKAYIHTRMRAKTSDFLKVLNRARPDAEKKEMKTITGKTFS']]
['jp_xEYVipE', ['XRCC6', 248, 27.04784131, 'DIISIAEDEDLRVHFEESSKLEDLLRKVRAKETRKRALSRLKLKLNKDIVISVGIYNLVQ']]
['jp_RlwniUC', ['NOP56', 155, 25.16169154, 'RLHFHNLVKGLTDLSACKAQLGLGHSYSRAKVKFNVNRVDNMIIQSISLLDQLDKDINTF']]
['jp_T4Vrdip', ['HSPA1B', 270, 22.7193807, 'FVEEFKRKHKKDISQNKRAVRRLRTACERAKRTLSSSTQASLEIDSLFEGIDFYTSITRA']]
['jp_Jk1Wcr7', ['C16orf72', 135, 114.2068966, 'DTHQRSFDIGIQIGYQRRNKDVLAWVKKRRRTIRREDLISFLCGKVPPPRNSRAPPRLTV']]
['jp_S4pkU1x', ['CLK2', 121, 50.93181818, 'DTDYRHSYEYQRENSSYRSQRSSRRKHRRRRRRSRTFSRSSSQHSSRRAKSVEDDAEGHL']]
['jp_xSUHZHI', ['TRABD', 54, 38.29457364, 'PVPRVLSGDPQNLSDVDAFNLLLEMKLKRRRQRPNLPRTVTQLVAEDGSRVYVVGTAHFS']]
['jp_8aCCoGG', ['SNRNP27', 27, 37.03977273, 'MGRSRSRSPRRERRRSRSTSRERERRRRERSRSRERDRRRSRSRSPHRRRSRSPRRH']]
['jp_XzDSsUl', ['RAB24', 132, 30.38297872, 'VKELRSLEEGCQIYLCGTKSDLLEEDRRRRRVDFHDVQDYADNIKAQLFETSSKTGQSVD']]
['jp_BcvPEDu', ['SURF6', 147, 27.42424242, 'RQRLHEKIQEARGQGSAKELSPAALEKRRRRKQERDRKKRKRKELRAKEKARKAEEATEA']]
['jp_3xiJMN2', ['HNRNPA3', 20, 24.4516129, 'MEVKPPPGRPQPDSGRRRRRRGEEGHDPKEPEQLRKLFIGGLSFETTDDS']]
['jp_PB5nJX6', ['ATF3', 92, 23.73913043, 'ESVTVSDRPLGVSITKAEVAPEEDERKKRRRERNKIAAAKCRNKKKEKTECLQKESEKLE']]
['jp_lQ0uMqc', ['CLK2', 122, 22.63636364, 'TDYRHSYEYQRENSSYRSQRSSRRKHRRRRRRSRTFSRSSSQHSSRRAKSVEDDAEGHLI']]
['jp_915q_0_', ['TCEA1', 218, 82.36363636, 'NRVRSRISNLKDAKNPNLRKNVLCGNIPPDLFARMTAEEMASDELKEMRKNLTKEAIREH']]
['jp_D3wnpYB', ['GTF2B', 298, 55.01652893, 'RTQKEIGDIAGVADVTIRQSYRLIYPRAPDLFPTDFKFDTPVDKLPQL\n']]
['jp_VJXds4o', ['CYC1', 196, 52.73529412, 'RPGKLFDYFPKPYPNSEAARAANNGALPPDLSYIVRARHGGEDYVFSLLTGYCEPPTGVS']]
['jp_NuYj7IF', ['PRKCSH', 296, 43.83333333, 'SFYDRVWAAIRDKYRSEALPTDLPAPSAPDLTEPKEEQPPVPSSPTEEEEEEEEEEEEEA']]
['jp_PfaeavZ', ['BIRC5', 53, 37.454039, 'NWPFLEGCACTPERMAEAGFIHCPTENEPDLAQCFFCFKELEGWEPDDDPIEEHKKHSSG']]
['jp_FYCmN9d', ['FKBP8', 208, 31.54861111, 'YGPQGRSPYIPPHAALCLEVTLKTAVDGPDLEMLTGQERVALANRKRECGNAHYQRADFV']]
['jp_pdmRz0M', ['PFAS', 881, 30.9339934, 'LYVALSPGQHRLGGTALAQCFSQLGEHPPDLDLPENLVRAFSITQGLLKDRLLCSGHDVS']]
['jp_IO8PNj2', ['GADD45A', 121, 29.88, 'VSNPGRLAELLLLETDAGPAASEGAEQPPDLHCVLVTNPHSSQWKDPALSQLICFCRESR']]
['jp_wa8CGje', ['CSE1L', 129, 25.92, 'MLSSPEQIQKQLSDAISIIGREDFPQKWPDLLTEMVNRFQSGDFHVINGVLRTAHSLFKR']]
['jp_WE62SWo', ['PPIL4', 171, 24.40594059, 'TFVDKDFVPYQDIRINHTVILDDPFDDPPDLLIPDRSPEPTREQLDSGRIGADEEIDDFK']]
['jp_gu_IYnt', ['SDHB', 49, 67.84965831, 'GACLQASRGAQTAAATAPRIKKFAIYRWDPDKAGDKPHMQTYEVDLNKCGPMVLDALIKI']]
['jp_ALAZD8m', ['LTA4H', 375, 57.30337079, 'GWGELQNSVKTFGETHPFTKLVVDLTDIDPDVAYSSVPYEKGFALLFYLEQLLGGPEIFL']]
['jp_xzJwefZ', ['NGRN', 39, 44.0754717, 'GGRVCAAVTRCGFATRGVAGPGPIGREPDPDSDWEPEERELQEVESTLKRQKQAIRFQKI']]
['jp_GJdhIML', ['HSP90B1', 751, 35.19191919, 'TLRSGYLLPDTKAYGDRIERMLRLSLNIDPDAKVEEEPEEEPEETAEDTTEDTEQDEDEE']]
['jp_nfqezln', ['GOT1', 29, 31.57627119, 'MAPPSVFAEVPQAQPVLVFKLTADFREDPDPRKVNLGVGAYRTDDCHPWVLPVVKKVEQ']]
['jp_fdhfiAo', ['METTL3', 559, 31.19463087, 'KIELFGRPHNVQPNWITLGNQLDGIHLLDPDVVARFKQRYPDGIISKPKNL\n']]
['jp_khXQWUJ', ['AAMP', 40, 25.18421053, 'ADTPPLETLSFHGDEEIIEVVELDPGPPDPDDLAQEMEDVDFEEEEEEEGNEEGWVLEPQ']]
['jp_h5pKeil', ['MAT2A', 50, 25.16949153, 'TSESVGEGHPDKICDQISDAVLDAHLQQDPDAKVACETVAKTGMILLAGEITSRAAVDYQ']]