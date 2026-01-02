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

import MySQLdb as mydb

import csv
import pickle

def main():
    cur = con_db()
    merops_code_mece = create_merops_code_table()
    print("merops_code_mece.")
    print(merops_code_mece)
    print(len(merops_code_mece))

    trim_len = 160 #50 # 108

    np.random.seed(42)
    subs_count_list = create_subs_count_list()


    df_full_aa_shorter_than_trim_len = pd.DataFrame({
        'protease_num': [], 'merops_id': [], 'substrate_num': [], 'uniprot_id': [], 'full_aa': [], 'full_aa_length':[], 'p1': []
    }) 

    error_list = []
    
    create_posi_nega_dataset(merops_code_mece, cur, error_list, df_full_aa_shorter_than_trim_len, trim_len)


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
    csv_file = open("./learning_data/subs_count_list.csv", "r", encoding="ms932", errors="", newline="" )
    #リスト形式
    f2 = csv.reader(csv_file, delimiter=",", doublequote=True, lineterminator="\r\n", quotechar='"', skipinitialspace=True)

    for row in f2:
        #print(row)
        subs_count_list = row

    print("subs_count_list.")
    print(subs_count_list)
    return subs_count_list

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

    print(df)
    print(df["sequence"][0]) 
    return df["sequence"][0]

def create_posi_nega_dataset(merops_code_mece, cur, error_list, df_full_aa_shorter_than_trim_len, trim_len):
    #count = 0
    # protease_turn = 0
    protease_num_dic = {'S01.247': 566, 'S08.071': 741}
    protease_turn = protease_num_dic['S08.071']

    #substrate_turn = 0
    #error_list = ["DAA158"]
    # for protease in merops_code_mece:
    # for protease in ['S01.247']:
    

    for protease in ['S08.071']: # furin protease
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
        
        # substrate_turn = 0
        
        df_cleave_pattern = pd.DataFrame(
            data={'protease_turn': [], 'merops_id': [], 'substrate_turn': [], 'uniprot_id':[], 'p1':[], 'len_full_aa': [], 'cleave_pattern':[], 'full_aa': []}
        )
        #display(df_cleave_pattern)
        print("df_cleave_pattern.")
        print(df_cleave_pattern)

        df_negative_pattern = pd.DataFrame(data={
            'protease_turn': [], 'merops_id': [], 'substrate_turn': [], 'uniprot_id':[], 'p1':[], 'negative_pattern':[], 'full_aa': [], 'len_full_aa': [], 'start': [], 'end': [], 'a[k]': []
        })

        #display(df_negative_pattern)    
        print("df_negative_pattern.")
        print(df_negative_pattern)

        df_uniprot_id_error = pd.DataFrame({
            'protease_turn': [], 'merops_id': [], 'substrate_turn': [], 'uniprot_id':[]
        })
        
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

        #list_subs = len(subs)
        list_subs = subs
        #list_subs = (('A3R530', 320), ('A3R534', 338), ('F2P1E0', 338), ('Q96T73', 255), ('P51170', 135), ('P51170', 136), ('P51170', 137), ('P51170', 138), ('P51170', 153), ('P51170', 168), ('P51170', 170), ('P51170', 172), ('P51170', 178), ('P51170', 179), ('P51170', 180), ('P51170', 181), ('P51170', 189), ('K9N5Q8', 887), ('P0DTC2', 815), ('P59594', 797), ('O15393', 255))
            
        #positive_dataの元となる配列データを作成する
        df_cleave_pattern = create_posi_dataset(protease, list_subs, df_cleave_pattern, protease_turn, error_list, df_full_aa_shorter_than_trim_len, df_uniprot_id_error, trim_len)

        #negative_dataの元となる配列データを作成する
        df_negative_pattern = create_nega_dataset(protease, list_subs, df_negative_pattern, protease_turn, error_list, df_full_aa_shorter_than_trim_len, df_uniprot_id_error, trim_len) 
    print("END")
    #return df_cleave_pattern, df_negative_pattern

def create_posi_dataset(protease, list_subs, df_cleave_pattern, protease_turn, error_list, df_full_aa_shorter_than_trim_len, df_uniprot_id_error, trim_len):
    print("list_subs: ")
    print(list_subs)
    for i in range(len(list_subs)):
        
        print(f"----------------------protease: {protease_turn}, {protease}; positive_data: {i+1}/{len(list_subs)}, {list_subs[i]}--------------------------")
        substrate_turn = i
        # if protease == 'S01.300':
        #     break
        #if protease_turn == 5:
        #    break
        
        
        # rを付けることを推奨。
        # バックスラッシュをそのままで分かりやすいため。
        # uniprot_id = list_subs[i][0].strip()
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
        uniprot_id = list_subs[i][0].strip()
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
            df_uniprot_id_error[f'{i}'] = [protease_turn, protease, substrate_turn, uniprot_id]
            continue
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
            df_uniprot_id_error[f'{i}'] = [protease_turn, protease, substrate_turn, uniprot_id]
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
        if len(full_aa) < trim_len:
            print("full_aa len < trim_len: {}.".format(trim_len))
            df_full_aa_shorter_than_trim_len.loc['f{substrate_turn}'] = [protease_turn, protease, substrate_turn, uniprot_id, full_aa, len(full_aa), p1]
            continue
        
        
        #例外である端っこも取得できるように工夫する
        if p1 - trim_len/2 >= 0 and len(full_aa) - p1 >= trim_len/2:
            cleave_pattern = full_aa[int(p1-trim_len/2) :int(p1+trim_len/2)]
        elif p1 - trim_len/2 < 0:
            term = int(trim_len/2 - p1)
            print("term: {}".format(term))
            print("type of term: {}".format(type(term)))
            cleave_pattern = "-"*term + full_aa[0:trim_len-term]
            print("cleave_pattern in p1 - trim_len/2 < 0: ")
            print(cleave_pattern)
        elif len(full_aa) - p1 < trim_len/2:
            term = int(trim_len/2 - (len(full_aa) - p1))
            cleave_pattern = full_aa[len(full_aa) - (trim_len-term):len(full_aa)] + "-"*term
            print("cleave_pattern in len(full_aa) - p1 < trim_len/2: ")
            print(cleave_pattern)
        else:
            pass
        
        #開裂パターンを表示する
        #print("cleave_pattern.")
        #print(cleave_pattern)
        
        #df_cleave_pattern.loc[f'{i}'] = [uniprot_id, p1, cleave_pattern]
        # data={'protease_turn': [], 'merops_id': [], 'substrate_turn': [], 'uniprot_id':[], 'p1':[], 'cleave_pattern':[], 'full_aa': []}
        df_cleave_pattern.loc[f'{i}'] = [protease_turn, protease, substrate_turn, uniprot_id, p1, len(full_aa), cleave_pattern, full_aa]

        #display(df_cleave_pattern)
        print("df_cleave_pattern.")
        print(df_cleave_pattern)

        substrate_turn = substrate_turn + 1

        filename = f'./proteases/cleave_pattern_one_letter_aa_{protease}.csv'
        df_cleave_pattern.to_csv(filename)
        print("df_cleave_pattern.")
        print(df_cleave_pattern)   

        filename_short_aa = f'./proteases/positive_pattern_aa_less_than_eight_{protease}.csv'
        df_full_aa_shorter_than_trim_len.to_csv(filename_short_aa)

        filename_uniprot_id_error = f'./proteases/positive_pattern_uniprot_id_error_{protease}.csv'
        df_uniprot_id_error.to_csv(filename_uniprot_id_error)
        
    filename = f'./proteases/cleave_pattern_one_letter_aa_{protease}.csv'
    df_cleave_pattern.to_csv(filename)
    print("df_cleave_pattern.")
    print(df_cleave_pattern)   

    filename_short_aa = f'./proteases/positive_pattern_aa_less_than_eight_{protease}.csv'
    df_full_aa_shorter_than_trim_len.to_csv(filename_short_aa)

    filename_uniprot_id_error = f'./proteases/positive_pattern_uniprot_id_error_{protease}.csv'
    df_uniprot_id_error.to_csv(filename_uniprot_id_error)

    substrate_turn = 0
    return df_cleave_pattern

def create_nega_dataset(protease, list_subs, df_negative_pattern, protease_turn, error_list, df_full_aa_shorter_than_trim_len, df_uniprot_id_error, trim_len):
    subs_memo = []
    p1_memo = []
    df_no_negative_data = pd.DataFrame({
        'protease_turn': [], 'merops_id': [], 'substrate_tun': [], 'uniprot_id': [], 'full_aa': [], 'full_aa_length': [] #, 'p1': []
    })
    num_no_negative_data = 0
    temp_negative_pattern_list = []
    num_negative_data = 0

    for i in range(len(list_subs)):
        print(f"----------------------protease: {protease_turn}, {protease}; negative_data: {i+1}/{len(list_subs)}, {list_subs[i]}--------------------------")
        substrate_turn = i

        uniprot_id = list_subs[i][0].strip()
        def uniprot_id_re(uniprot_id):
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
                df_uniprot_id_error[f'{i}'] = [protease_turn, protease, substrate_turn, uniprot_id]
                # continue
            #else:
            #    pass
            uniprot_id_re = result[0]
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
                df_uniprot_id_error[f'{i}'] = [protease_turn, protease, substrate_turn, uniprot_id]
                #continue
            return uniprot_id_re
        uniprot_id = uniprot_id_re(uniprot_id)

        p1 = list_subs[i][1]
        print(f"p1: {p1}")


        subs_memo.append(list_subs[i][0])
        p1_memo.append(list_subs[i][1])
        if i > 0:
            print(f"list_subs[i][0]: {list_subs[i][0]}")
            print(f"list_subs[i-1][0]: {list_subs[i-1][0]}")
            print(f"list_subs[i][1]: {list_subs[i][1]}")
            print(f"list_subs[i-1][1]: {list_subs[i-1][1]}")
            if list_subs[i][0] == list_subs[i-1][0]:
                continue
            elif list_subs[i][0] != list_subs[i-1][0] and len(subs_memo) >= 3: 
                print("#"*10+"i: {}, ".format(i)+"make negative patterns for {}, {}".format(protease, uniprot_id))
                print("subs_memo: ")
                last_subs_memo = subs_memo.pop(-1)
                print(subs_memo)
                print(len(subs_memo))
                print(last_subs_memo)
                print("p1_memo: ")
                last_p1_memo = p1_memo.pop(-1)
                print(p1_memo)
                print(len(p1_memo))
                print(last_p1_memo)

                uniprot_id = list_subs[i-1][0].strip()
                uniprot_id = uniprot_id_re(uniprot_id)
                full_aa = aaseq_from_uid(uniprot_id, protease_turn-1, substrate_turn-1)
                df_full_aa = pd.DataFrame(list(full_aa), columns=["char"])
                print("df_full_aa: ")
                print(df_full_aa)

                index_list = list(range(0, len(full_aa)))
                df_index = pd.DataFrame({
                    'index_list': index_list
                })
                cond_list = []
                for j in range(len(p1_memo)):
                    cond = (df_index['index_list'] <  p1_memo[j]-trim_len/2) | (df_index['index_list'] > p1_memo[j]+trim_len/2)
                    cond_list.append(cond)
                #and_cond_list = []
                #for j in range(len(cond_list)):
                #    temp_cond_list = cond_list[j]
                #    for k in range():
                #result = [all(values) for values in zip(*lists)]
                and_cond_list = [all(values) for values in zip(*cond_list)]
                print("and_cond_list: ")
                print(and_cond_list)
                print(len(and_cond_list))

                not_positive_regions = df_full_aa[and_cond_list]
                print("not_positive_regions: ")
                print(not_positive_regions)
                print(len(not_positive_regions))

                remain_index_list = df_index[and_cond_list]
                print("remain_index_list: ")
                print(remain_index_list)
                remain_index_list = df_index['index_list'].tolist()

                """
                nums = [1, 2, 5, 6, 7, 8, 9, 14, 19, 31, 33, 35, 36, 37, 38, 39, 41, 45]
                ranges = []
                start = 0
                for i in range(1, len(nums) + 1):
                    # 連続が途切れた時
                    if i == len(nums) or nums[i] != nums[i-1] + 1:
                        length = i - start
                        if length >= 5:
                            ranges.append((start, i - 1))  # (開始インデックス, 終了インデックス)
                        start = i
                print(ranges)
                """
                ranges = []
                start = 0
                for j in range(1, len(remain_index_list) + 1):
                    # 連続が途切れた時
                    if j == len(remain_index_list) or remain_index_list[j] != remain_index_list[j-1] + 1:
                        length = j - start
                        if length >= trim_len:
                            # ranges.append((start, j - 1))  # (開始インデックス, 終了インデックス)
                            ranges.append((remain_index_list[start], remain_index_list[j - 1]))  # (開始インデックス, 終了インデックス)
                        start = j
                print("ranges: ")
                print(ranges)
                print(len(ranges))

                if len(ranges) == 0:
                    filename_no_negative_data = './proteases/df_no_nagative_data.csv'
                    df_no_negative_data.loc[f'{num_no_negative_data}'] = [protease_turn, protease, substrate_turn, uniprot_id, full_aa, len(full_aa)] #, p1_memo]
                    df_no_negative_data.to_csv(filename_no_negative_data)
                    num_no_negative_data += 1
                    subs_memo = []
                    p1_memo = []
                    continue
                else:
                    base_num = i-1 - len(ranges) 
                    for j in range(len(ranges)):
                        print("ranges[j]: {}".format(ranges[j]))
                        temp_negative_pattern = ranges[j]
                        print('temp_negative_pattern: ')
                        print(temp_negative_pattern)
                        temp_negative_pattern_list.append(temp_negative_pattern)
                        pd.DataFrame({
                            'temp_negative_data': temp_negative_pattern_list, 
                            'uniprot_id': uniprot_id
                        }).to_csv(f'./proteases/temp_negative_data_{protease}.csv')
                        filename_negatiev_data = f'./proteases/negative_pattern_one_letter_aa_{protease}.csv'
                        print('full_aa: ')
                        print(full_aa)
                        
                        random.seed(42)
                        a = np.arange(temp_negative_pattern[0], temp_negative_pattern[1]-trim_len)
                        np.random.shuffle(a)
                        repeat_num = (temp_negative_pattern[1] - temp_negative_pattern[0])//trim_len
                        for k in range(repeat_num): 
                            #df_negative_pattern.loc[f'{base_num + j}'] = [protease_turn, protease, substrate_turn, uniprot_id, p1, full_aa[a[0]: a[0]+50], full_aa, ranges[j][0], ranges[j][1]]
                            df_negative_pattern.loc[f'{num_negative_data}'] = [protease_turn, protease, substrate_turn, uniprot_id, p1, full_aa[a[k]: a[k]+trim_len], full_aa, len(full_aa), ranges[j][0], ranges[j][1], a[k]]
                            num_negative_data += 1
                            df_negative_pattern.to_csv(filename_negatiev_data)
                    subs_memo = []
                    p1_memo = []
                    continue
                #subs_memo = []
                #p1_memo = []
            else:
                subs_memo = []
                p1_memo = []
                subs_memo.append(list_subs[i][0])
                p1_memo.append(list_subs[i][1])
                pass
        else:
            pass
        full_aa = aaseq_from_uid(uniprot_id, protease_turn, substrate_turn)
        #cleave_pattern = full_aa[p1-3:p1+4]
        #print(cleave_pattern)
        
        print(f"len_full_aa: {len(full_aa)}")
        if len(full_aa) < trim_len:
            #'protease_num': [], 'merops_id': [], 'substrate_num': [], 'uniprot_id': [], 'full_aa': [], 'full_aa_length': []
            #df_full_aa_shorter_than_trim_len = pd.DataFrame({
            #    'protease_num': [], 'merops_id': [], 'substrate_num': [], 'uniprot_id': [], 'full_aa': [], 'full_aa_length':[], 'p1': []
            #}) 
            df_full_aa_shorter_than_trim_len.loc['f{substrate_turn}'] = [protease_turn, protease, substrate_turn, uniprot_id, full_aa, len(full_aa), p1]
            continue
        
        # np.random.seed(i)
        random.seed(42)
        #a = np.arange(len(full_aa)-1)
        a = np.arange(len(full_aa))

        print(f"len_a :{len(a)}")
        
        for j in range(trim_len):
            # print(f"j: {j}")
            if p1- trim_len/2 + j < 0:
                continue
            elif p1 - trim_len/2 + j > len(full_aa)-1:
                continue
            a = a[a != p1 - trim_len/2 + j]
        print("p1: {}".format(p1))
        print("a after eliminate around p1: ")
        print(a)
        print(len(a))

        ranges = []
        start = 0
        remain_index_list = a
        for j in range(1, len(remain_index_list) + 1):
            # 連続が途切れた時
            if j == len(remain_index_list) or remain_index_list[j] != remain_index_list[j-1] + 1:
                length = j - start
                if length >= trim_len:
                    # ranges.append((start, j - 1))  # (開始インデックス, 終了インデックス)
                    ranges.append((remain_index_list[start], remain_index_list[j - 1]))  # (開始インデックス, 終了インデックス)
                start = j
        print("ranges: ")
        print(ranges)
        print(len(ranges))
        print("len of full aa: ")
        print(len(full_aa))
        print("p1: ")
        print(p1)

        for j in range(len(ranges)):
            print("ranges[j]: {}".format(ranges[j]))
            temp_negative_pattern = ranges[j]
            print('temp_negative_pattern: ')
            print(temp_negative_pattern)
            temp_negative_pattern_list.append(temp_negative_pattern)
            pd.DataFrame({
                'temp_negative_data': temp_negative_pattern_list, 
                'uniprot_id': uniprot_id
            }).to_csv(f'./proteases/temp_negative_data_{protease}.csv')
            filename_negatiev_data = f'./proteases/negative_pattern_one_letter_aa_{protease}.csv'
            print('full_aa: ')
            print(full_aa)
            
            random.seed(42)
            a = np.arange(temp_negative_pattern[0], temp_negative_pattern[1]-trim_len)
            np.random.shuffle(a)
            repeat_num = (temp_negative_pattern[1] - temp_negative_pattern[0])//trim_len
            print("repeat_num in only one uniprot_id: ")
            print(repeat_num)
            for k in range(repeat_num): 
                #df_negative_pattern.loc[f'{base_num + j}'] = [protease_turn, protease, substrate_turn, uniprot_id, p1, full_aa[a[0]: a[0]+50], full_aa, ranges[j][0], ranges[j][1]]
                df_negative_pattern.loc[f'{num_negative_data}'] = [protease_turn, protease, substrate_turn, uniprot_id, p1, full_aa[a[k]: a[k]+trim_len], full_aa, len(full_aa), ranges[j][0], ranges[j][1], a[k]]
                num_negative_data += 1
                df_negative_pattern.to_csv(filename_negatiev_data)
        """
        print("a.")
        print(a)   
        random.seed(42)
        np.random.shuffle(a)

        print("a.")
        print(a)
        print(a[0])

        
        #例外である端っこも取得できるように工夫する
        if len(full_aa)-1 - a[0] >= trim_len-1:
            negative_pattern = full_aa[a[0]:a[0]+trim_len]
        else:
            term = int(trim_len - (len(full_aa)-1 - a[0]))
            negative_pattern = full_aa[a[0]:len(full_aa)] + "-"*term
        
        #negative_pattern = full_aa[a[0]:a[0]+8]
        
        print("negative_pattern.")
        print(negative_pattern)
        
        print("a[0]+trim_len/2: {}".format(a[0]+trim_len/2))
        print("p1: {}".format(p1))
        #data={'protease_turn': [], 'merops_id': [], 'substrate_turn': [], 'uniprot_id':[], 'p1':[], 'negative_pattern':[]}
        #df_negative_pattern.loc[f'{i}'] = [protease_turn, protease, substrate_turn, uniprot_id, a[0]+trim_len/2, negative_pattern]
        df_negative_pattern.loc[f'{num_negative_data}'] = [protease_turn, protease, substrate_turn, uniprot_id, p1, negative_pattern, full_aa, p1-trim_len/2, p1+trim_len/2]
        num_negative_data += 1
        print("df_negative_pattern.")
        print(df_negative_pattern)
        """

        substrate_turn = substrate_turn + 1

        filename = f'./proteases/negative_pattern_one_letter_aa_{protease}.csv'
        df_negative_pattern.to_csv(filename)

        filename_short_aa = f'./proteases/negative_pattern_aa_shorter_than_trim_len_{protease}.csv'
        df_full_aa_shorter_than_trim_len.to_csv(filename_short_aa)

        filename_uniprot_id_error = f'./proteases/negative_pattern_uniprot_id_error_{protease}.csv'
        df_uniprot_id_error.to_csv(filename_uniprot_id_error)
        
    filename = f'./proteases/negative_pattern_one_letter_aa_{protease}.csv'
    df_negative_pattern.to_csv(filename)

    filename_short_aa = f'./proteases/negative_pattern_aa_shorter_than_trim_len_{protease}.csv'
    df_full_aa_shorter_than_trim_len.to_csv(filename_short_aa)

    filename_uniprot_id_error = f'./proteases/negative_pattern_uniprot_id_error_{protease}.csv'
    df_uniprot_id_error.to_csv(filename_uniprot_id_error)

    #count = count + 1
    return  df_negative_pattern    


#if __name__ == '__name__':
#    main()

if __name__ == "__main__":
   main()
   print("END.")


















