import xml.etree.ElementTree
from argparse import ArgumentParser
def init_parser():
  parser = ArgumentParser()
  parser.add_argument('--config', type=str, help='Detector config XML file (full path)', default='/mnt/research/NuInt/DUNEPrismTools/jake-contrib/try.xml')
  parser.add_argument('--DIR', type=str, help='Directory for plots', default='')
  parser.add_argument('--setting', type=str, help='Horn Current (FD?)', default='FHC')

  return parser

#print a.config

