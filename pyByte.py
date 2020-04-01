import binascii
import struct
from struct import *
#def test_write_bin(self):  
fname = 'tt.bin'
#f = open(fname, 'wb')
## ff d8 ff e0 00 10 4a
#f.write(bytearray([255, 216, 255, 224, 0, 16, 74]))
#f.close()
#00 00 '00 00 '4F 81 81 08 00 00 00 00 '33 00 '33 00 '00 40 00 00
with open(fname, "rb") as binfile :
  line = binfile.readline()
  for line in binfile:
    kk = bytearray(line[1:5])
    board, channel = struct.unpack('>HH',kk) 
    print board, channel
    tt = bytearray(line[5:13])
# #    #tt = tt[::-1] # invert
    timestamp = struct.unpack('>q',tt)# big-endian, q: 8-byte long long
    print hex(timestamp[0])
#    ck = struct.unpack('cccccccc',tt)
#    time = timestamp[::-1]
#   bt = tt.decode('hex')
#    for i in time:
#      ss = hex(i)
#      print ss
#    #format = '%dB' % len(ss)
#    #d = struct.unpack_from(format,ss)
#    print timestamp#int('0x'+ss,16)
#    ee = line[13:15]
#    elong = struct.unpack('BB',ee) 
#    print elong
#    print len(line)





#with open(fname, "rb") as f1:
#  while True:
#    current_byte = f1.read(1)
#    if (not current_byte):
#      break
#    val = ord(current_byte)
#    print hex(val),
#print
