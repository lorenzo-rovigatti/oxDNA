import os
from datetime import datetime
import random
import sys
# ランダムな長さ10のsequenceを作る

def make_filename():
    now = datetime.now()
    date_time_str = now.strftime("%Y-%m-%d-%H%M%S")
    path = f'seq/seq_{date_time_str}.dat'
    return path

def make_file(path):
    f = open(path, 'w')
    return f

def make_random_seq(len):
    seq = ''
    amide = ['A', 'T', 'G', 'C']
    for i in range(len):
        tmp = random.random()*4
        seq += amide[int(tmp)]
    return seq

def main():
    # ファイルを作成
    path = make_filename()
    print(path)
    f = make_file(path)
    # ランダムな文字列の生成
    sequence = make_random_seq(10)
    # ファイルに書き込む
    f.write(sequence)
    f.close()

if __name__ == '__main__':
    # print(sys.argv[1])
    main()