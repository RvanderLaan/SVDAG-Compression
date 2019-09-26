import os
import numpy as np

# def int_arr_from_bytes(bytes):

class Node():
  def __init__(self):


def parse_int(f):
  bts = f.read(4)
  return int.from_bytes(bts, byteorder='little')

def parse_int3(f):
  res = []
  res.append(parse_int(f))
  res.append(parse_int(f))
  res.append(parse_int(f))
  return res

def load_svo(filename: str):
  with open(filename, 'rb') as f:
    bbox_min = parse_int3(f)
    bbox_max = parse_int3(f)
    root_side = parse_int(f)
    levels = parse_int(f)
    count = parse_int(f)

    print('File {} has {} levels and {} nodes'.format(filename, levels, count))

    

    np.array()

load_svo('./models/sponza-10-lossless.svdag')
