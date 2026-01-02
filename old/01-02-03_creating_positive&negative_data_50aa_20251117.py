# Spike proteinを開裂するproteasesの基質たちを学習して，新たなSpike protein を開裂するproteases を探す機械学習モデルです．

# merops からcopy&pasteで作ったproteaseの基質たちの情報が入ったファイルの名前を配列に格納します．
#ここもいずれは自動化したい．
#おそらく，merops databaseを構築して，そこから情報を検索するようにすればいいのかと．

#2022/07/10現在
#MEROPSのDBは構築済みなので，Uniprot IDのリストをMEROPS DBから作成すればよくなりました．
#Web上のMEROPSでは存在する基質データが構築したDBにはまだ反映されていなくて存在しない例もあるようなので注意です．


#1
import re
import numpy as np
import glob
import random
import os
import matplotlib.pyplot as plt
import numpy as np

#ライブラリのインポート
import pandas as pd
from urllib.request import urlopen
import numpy as np
from lxml import etree

#MEROPSはmysqlでもう動くので本当はここは使いたくない．
#2022/07/26現在訂正します．MySQLによるMROPSデータベースでは，更新が間に合っていないデータがあるので，
# ウェブからコピペで取得したものを使うべきです．たった４つしかポジティブデータがないので，
# 手動でできる．この時，Dataframeによってコピペしたファイルを扱います．

# mysqlに接続します．
#ここからはじまる．
#import mysql.connector as mydb
import MySQLdb as mydb

import csv
import pickle

def main():
    cur = con_db()
    merops_code_mece = create_merops_code_table()
    print("merops_code_mece.")
    print(merops_code_mece)
    print(len(merops_code_mece))

    np.random.seed(0)
    subs_count_list = create_subs_count_list()

    df_cleave_pattern = pd.DataFrame(
    data = {'uniprot_id':[], 'cleave_pattern':[]}
    )
    print("df_cleave_pattern")
    print(df_cleave_pattern) 
    df_full_aa_less_than_eight = pd.DataFrame({
        'protease_num': [], 'merops_id': [], 'substrate_num': [], 'uniprot_id': [], 'full_aa': [], 'full_aa_length': []
    }) 

    #count = 0
    #protease_turn = 0
    #substrate_turn = 0
    # error_list = ["DAA158"]
    error_list = []
    
    create_posi_nega_dataset(merops_code_mece, cur, error_list, df_full_aa_less_than_eight)


def con_db():
    # コネクションの作成
    conn = mydb.connect(
        host='localhost',
        port=3306,
        user='root',
        password='miyazakilab',
        #database='meropsrefs01'
        database='meropsweb12_1'
    )

    # DB操作用にカーソルを作成
    cur = conn.cursor()
    return cur

# MySQLで使うMEROPS code辞書作成
#HATL1は開裂する基質情報が登録されていなかった．（2022年7月24日現在）
#TMPRSS13はMSPLとも呼ぶ． MSPLは4つの基質情報があったが，Uniprot IDが記載されていなかった．(2022年7月24日現在)

#import csv

def create_merops_code_table():
    csv_file = open("./learning_data/merops_code_mece.csv", "r", encoding="ms932", errors="", newline="" )
    #リスト形式
    f = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

    print(f)
    #print(len(f))

    i = 0
    for row in f:
        print(row)
        merops_code_mece = row
        i = i + 1
        if i == 1:
            break

    print("merops_code_mece")
    print(merops_code_mece)
    len(merops_code_mece)
    return merops_code_mece


def create_subs_count_list():
    csv_file = open("./leaning_data/subs_count_list.csv", "r", encoding="ms932", errors="", newline="" )
    #リスト形式
    f2 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

    for row in f2:
        #print(row)
        subs_count_list = row

    print("subs_count_list.")
    print(subs_count_list)
    return subs_count_list


# Uniprot APIを用いてアミノ酸配列を取得するスクリプトを用意します．
#https://seiyakukenkyusya.com/programming/collecting-uniprot-information/
#上記のリンク先を参考にしました．
#uid_eg = 'P78365'

def aaseq_from_uid(uid, protease_turn, substrate_turn):
    df = pd.DataFrame(np.arange(3).reshape(1, 3), columns=['uniprotKB_accession', 'function', 'sequence'], index=['protease'+str(protease_turn)+'_substrate'+str(substrate_turn)])

    for column_name in df:
        df[column_name] = df[column_name].astype(str)

    df['uniprotKB_accession'][0] = uid

    #display(df)    

    url = "https://www.uniprot.org/uniprot/" + uid + ".xml"
    f = urlopen(url)
    xml = f.read()
    root = etree.fromstring(xml)
    
    #以下のコードは下の説明を参照
    function = root.find('./entry/comment[@type="function"]', root.nsmap)
    if function==None:
        print("function was not detected.")
        pass
    else:
        df["function"][0] = function[0].text
        #print(function[0].text+"¥n")

    sequence = root.find('./entry/sequence', root.nsmap) 
    if sequence==None: 
        print("sequence was not detected.")
        pass 
    else: 
        df["sequence"][0] = sequence.text 
        #print(sequence.text+"¥n")

    #display(df) 
    #df0=pd.concat([df0, df], axis=0)
    #display(df0)
    print(df)
    print(df["sequence"][0]) 
    return df["sequence"][0]
    
    


#print("aaseq_from_uid(uid_eg, 0, 0): {}".format(aaseq_from_uid(uid_eg, 0, 0)))


# 取得したアミノ酸からP1, P1'を特定する．今回はP4部位を戻り値にした．

# ファイルを各々読み込んで，Uniprot IDからアミノ酸配列を取得します．
#&nbsp;取得したアミノ酸配列の開裂部位を特定して，開裂部位のP1, P1' を中心に，開裂部位左に10アミノ酸，開裂部位右に10アミノ酸だけ取り出して，1or0の値，電荷値，親水性値を階数3のテンソルに代入しています．値が入っている最奥のベクトルの次元が３次元で，cahannel数が3と呼びます．
#行数がアミノ酸の種類である20になり，列数が取得したデータのアミノ酸である20になり，20×20型行列ができることになる．
#例えると，カラー画像の入力がR, G, Bであるのと同じになっています．<br>
#A = [<br>
#        &emsp;[<br>
#          &emsp;&emsp;[R, G, B], [R, G, B], [R, G, B]<br>
#           &emsp;], <br>
#        &emsp;[<br>
#          &emsp;&emsp;[R, G, B], [R, G, B], [R, G, B]<br>
#            &emsp;],<br> 
#        &emsp;[<br>
#          &emsp;&emsp;[R, G, B], [R, G, B], [R, G, B]<br>
#            &emsp;], <br>
#                                                   ]<br>
#A[0][0][0] == [R, G, B]<br>
#というような形の階数3のテンソルです．階数3とは，各値に基底が3種類あって，その基底3つのすべての組み合わせでできる基底でできた空間を階数3のテンソルと呼びます．<br>
#基底3種類をx, y, zとすると，(x or y or z) × (x or y or z) × (x or y or z)で各括弧内のグループから1つずつ取り出した基底ができて，27種類できることになる．例えば，ある階数3のテンソルの基底はAxxz, Ayzx, Ayyxなどが考えられます．<br>
#<br>
#https://www.youtube.com/watch?v=f5liqUk0ZTw<br>
#の10:30を参照ください．<br>
#<br>
#よって，今回の3階テンソルの場合，今回の20×20型行列であり，各ベクトルには3つの値が3channelとして入っているので，次のように43種類の基底があって，(20種類の基底)×(20種類の基底)×(3種類)=1200パターンの基底を持つ3階テンソルと考えられます．<br>
#<br>
#&emsp;&emsp;x_train[基質i][20種類のアミノ酸][取得したアミノ酸1つずつが20個][1or0, 電荷値，親水性値の各値一個ずつ]<br>
#<br>
#以上のような3階テンソル（半変テンソルか，共変テンソルか，どちらなのかは後ほど調べます．）を訓練データとして作成します．<br>

## 注意!! Uniprot IDは正規表現で[A-Z]{1}[0-9A-Z]{5}のように表す．
#だけど，文字列"MERNUM"とかがヒットするのはどう対処するのか．
#https://trade-and-develop.hatenablog.com/entry/2017/02/23/021119

#cur = conn.cursor(buffered=True)
#cur = conn.cursor()

# creating positive data 
#print("merops_code_mece.")
#print(merops_code_mece)
#print(len(merops_code_mece))


#main function

def create_posi_nega_dataset(merops_code_mece, cur, error_list, df_full_aa_less_than_eight):
    #count = 0
    protease_turn = 0
    #substrate_turn = 0
    #error_list = ["DAA158"]
    for protease in merops_code_mece:
        print("==================================protease"+str(protease_turn)+"=======================================")
        print("protease.")
        print(protease)
        print("protease_turn.")
        print(protease_turn)

        #if protease_turn < 858:
        #    protease_turn = protease_turn + 1
        #    continue
        protease_turn = protease_turn + 1
            
        #if protease in already_read:
        #    continue
        
        substrate_turn = 0
        
        df_cleave_pattern=pd.DataFrame(
            data={'', 'uniprot_id':[], 'p1':[], 'cleave_pattern':[]}
        )
        #display(df_cleave_pattern)
        print("df_cleave_pattern.")
        print(df_cleave_pattern)

        df_negative_pattern=pd.DataFrame(
        data={'uniprot_id':[], 'p1':[], 'negative_pattern':[]}
        )
        #display(df_negative_pattern)    
        print("df_negative_pattern.")
        print(df_negative_pattern)
        
        #func1
        #cleave_pattern = ""
        #cleave_pattern = cptn()       

        merops_code = [protease]
        print("this is ok.")

        #proteaseに紐づけられている基質を取得する
        print("merops_code:" + protease)
        cur.execute("SELECT uniprot_acc, p1 FROM cleavage where code=(%s);", merops_code)
        subs = cur.fetchall()
        print("subs.")
        print(subs)
        for sub in subs:
            print("sub.")
            print(sub)

        print("subs.")
        print(len(subs))

        """
        if len(subs) > 1000:
            len_subs = 1000
            
            np.random.seed(0)
            a = np.arange(len(subs))
            np.random.shuffle(a)
            b = a[0:1000]
            print("b.")
            print(b)
            print(len(b))
            
            list_subs = []
            for i in b:
                list_subs.append(subs[i])
        else:
            len_subs = len(subs)
            list_subs = subs
        
        if len_subs == 0:
            continue
        """    
        list_subs = len(subs)
        
            
        #positive_dataの元となる配列データを作成する
        df_cleave_pattern = create_posi_dataset(protease, list_subs, df_cleave_pattern, protease_turn, error_list, df_full_aa_less_than_eight)

        #negative_dataの元となる配列データを作成する
        df_negative_pattern = create_nega_dataset(protease, list_subs, df_cleave_pattern, protease_turn, error_list, df_full_aa_less_than_eight) 
    print("END")
    #return df_cleave_pattern, df_negative_pattern

def create_posi_dataset(protease, list_subs, df_cleave_pattern, protease_turn, error_list, df_full_aa_less_than_eight):
    for i in range(len(list_subs)):
        print(f"----------------------protease: {protease_turn}, positive_data: {i}/{len(list_subs)}--------------------------")
        substrate_turn = i
        if protease == 'S01.300':
            break
        #if protease_turn == 5:
        #    break
        
        
        # rを付けることを推奨。
        # バックスラッシュをそのままで分かりやすいため。
        uniprot_id = list_subs[i][0].strip()
        """
        print(f"uid: {uniprot_id}")
        content = fr'{uniprot_id}' 
        pattern = '.*?([A-Z]{1})([A-Z0-9]{2})([A-Z0-9]{3})'

        result = re.match(pattern, content)

        if result: #none以外の場合
            print("result.")
            print(result) 
            # output:<_sre.SRE_Match object; span=(0, 3), match='hel'>
            print("result.span().")
            print(result.span()) 
            # output:(0, 3)
            print("result.group().")
            print(result.group()) 
            # output:hel
        else:   
            print("result.")
            print(result) 
            continue
            
        uniprot_id = result.group(1) + result.group(2) + result.group(3)
        print(f"uniprot_id: {uniprot_id}")

        """

        #if uniprot_id in error_list:
        #    print("error_list.")
        #    print(error_list)
        #    break
        #uniprot_id = uniprot_id.strip()
        print(f"uniprot_id: {uniprot_id}")
        content = fr'{uniprot_id}' 

        #pattern = '.*?([A-Z]{1})([A-Z0-9]{2})([A-Z0-9]{3})'
        pattern01 = '.*?([A-Z0-9]{6})$'
        pattern02 = '(^[A-Z0-9]{6}.*?)'
        result01 = re.match(pattern01, content)
        print("result01: ")
        print(result01)
        #result01 = result01[0].strip()
        result02 = re.match(pattern02, content)
        print("result02")
        print(result02)
        #result02 = result02[0].strip()
        if result01 is not None:
            if len(result01[0]) == 6:
                result = result01
        elif result02 is not None:
            if len(result02[0]) == 6:
                result = result02
        else:
            ("re error!")
            pass
        #else:
        #    pass
        uniprot_id = result[0]
        print("uniprot_id after re: {}".format(uniprot_id))
        print("result.")
        print(result) 

        if result: #none以外の場合
            #print("result.")
            #print(result) 
            # output:<_sre.SRE_Match object; span=(0, 3), match='hel'>
            #print("result.span().")
            #print(result.span()) 
            # output:(0, 3)
            #print("result.group().")
            #print(result.group()) 
            # output:hel
            pass
        else:   
            #print("result.")
            #print(result) 
            continue
        #p1 = sub[j][1]
        #print(f"p1: {p1}") 
        
        p1 = list_subs[i][1]
        print(f"p1: {p1}")

        #アミノ酸配列全長を取得する
        full_aa = aaseq_from_uid(uniprot_id, protease_turn, substrate_turn)
        print("full_aa.")
        print(full_aa)
        print(type(full_aa))
        #print(f"length of full_aa: "+str(len(full_aa)))
        #if len(full_aa) < 3:
        #    print("full_aa len < 3.")
        #    continue
        
        
        #例外である端っこも取得できるように工夫する
        #if p1 - 4 >= 0 and len(full_aa) - p1 >= 4:
        #    cleave_pattern = full_aa[p1-4:p1+4]
        #elif p1 - 4 < 0:
        #    term = 4 - p1
        #    cleave_pattern = "-"*term + full_aa[0:8-term]
        #elif len(full_aa) - p1 < 4:
        #    term = 4 - (len(full_aa) - p1)
        #    cleave_pattern = full_aa[len(full_aa) - (8-term):len(full_aa)] + "-"*term
        #else:
        #    pass
        
        #開裂パターンを表示する
        #print("cleave_pattern.")
        #print(cleave_pattern)
        
        #df_cleave_pattern.loc[f'{i}'] = [uniprot_id, p1, cleave_pattern]
        df_cleave_pattern.loc[f'{i}'] = [uniprot_id, p1, full_aa]
        #display(df_cleave_pattern)
        print("df_cleave_pattern.")
        print(df_cleave_pattern)

        substrate_turn = substrate_turn + 1
        
    filename = f'./proteases/cleave_pattern_one_letter_aa_{protease}.csv'
    df_cleave_pattern.to_csv(filename)
    print("df_cleave_pattern.")
    print(df_cleave_pattern)   

    filename_short_aa = f'./proteases/positive_pattern_aa_less_than_eight_{protease}.csv'
    df_full_aa_less_than_eight.to_csv(filename_short_aa)

    substrate_turn = 0
    return df_cleave_pattern

def create_nega_dataset(protease, list_subs, df_negative_pattern, protease_turn, error_list, df_full_aa_less_than_eight):
    for i in range(len(list_subs)):
        print(f"----------------------protease: {protease_turn}, negative_data: {i}/{len(list_subs)}--------------------------")
        substrate_turn = i

        # rを付けることを推奨。
        # バックスラッシュをそのままで分かりやすいため。
        uniprot_id = list_subs[i][0]
        print(f"uid: {uniprot_id}")
        content = fr'{uniprot_id}' 
        pattern = '.*?([A-Z]{1})([A-Z0-9]{2})([A-Z0-9]{3})'

        result = re.match(pattern, content)

        if result: #none以外の場合
            print(result) 
            # output:<_sre.SRE_Match object; span=(0, 3), match='hel'>
            print(result.span()) 
            # output:(0, 3)
            print(result.group()) 
            # output:hel
        else:
            continue
            
        uniprot_id = result.group(1) + result.group(2) + result.group(3)
        print(f"uniprot_id: {uniprot_id}")
        
        if uniprot_id in error_list:
            break
        
        p1 = list_subs[i][1]
        print(f"p1: {p1}")
        
        full_aa = aaseq_from_uid(uniprot_id, protease_turn, substrate_turn)
        #cleave_pattern = full_aa[p1-3:p1+4]
        #print(cleave_pattern)
        
        print(f"len_full_aa: {len(full_aa)}")
        if len(full_aa) < 8:
            #'protease_num': [], 'merops_id': [], 'substrate_num': [], 'uniprot_id': [], 'full_aa': [], 'full_aa_length': []
            df_full_aa_less_than_eight.loc['f{substrate_turn}'] = [protease_turn, protease, substrate_turn, uniprot_id, full_aa, len(full_aa)]
            continue
        
        np.random.seed(i)
        a = np.arange(len(full_aa)-1)

        print(f"len_a :{len(a)}")
        
        for j in range(8):
            print(f"j: {j}")
            if p1- 4 + j < 0:
                continue
            elif p1 - 4 + j > len(full_aa)-1:
                continue
            a = a[a != p1 - 4 + j]

        print("a.")
        print(a)   
        np.random.shuffle(a)

        print("a.")
        print(a)
        print(a[0])

        
        #例外である端っこも取得できるように工夫する
        if len(full_aa)-1 - a[0] >= 8-1:
            negative_pattern = full_aa[a[0]:a[0]+8]
        else:
            term = 8 - (len(full_aa)-1 - a[0])
            negative_pattern = full_aa[a[0]:len(full_aa)] + "-"*term
        
        #negative_pattern = full_aa[a[0]:a[0]+8]
        
        print("negative_pattern.")
        print(negative_pattern)
            
        df_negative_pattern.loc[f'{i}'] = [uniprot_id, a[0]+4, negative_pattern]
        #display(df_negative_pattern)
        print("df_negative_pattern.")
        print(df_negative_pattern)

        substrate_turn = substrate_turn + 1
        
    filename = f'./proteases/negative_pattern_one_letter_aa_{protease}.csv'
    df_negative_pattern.to_csv(filename)

    filename_short_aa = f'./proteases/negative_pattern_aa_less_than_eight_{protease}.csv'
    df_full_aa_less_than_eight.to_csv(filename_short_aa)

    count = count + 1
    return  df_negative_pattern    


#if __name__ == '__name__':
#    main()

if __name__ == "__main__":
   main()
   print("END.")


















