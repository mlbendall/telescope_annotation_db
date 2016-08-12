#! /usr/bin/env python

"""
IGV.py

Provides python interface for controlling IGV through a port

"""

from socket import socket

class IGV:
  def __init__(self, host='0.0.0.0', port=60151):
    self.host = host
    self.port = port
  
  def send(self, command, confirm=True):
    s = socket()
    s.connect((self.host,self.port))
    s.send('%s\n' % command.strip('\n'))
    response = s.recv(100)
    if confirm:
      assert response.startswith('OK\n'), "ERROR: Response is %s" % response 
    s.close()
  
  def new(self):
    ''' Create a new session. Unloads all tracks except the default genome annotations
    '''
    self.send("new")
  
  def load(self, path):
    ''' Loads data or session files. Specify a comma-delimited list of full paths or URLs
    '''
    if not isinstance(path,str):
      self.send('load %s' % ','.join(path))
    else:
      self.send("load %s" % path)    

  def exit(self):
    ''' Exit (close) the IGV application
    '''
    self.send("exit",confirm=False)
  
  def genome(self,genome):
    ''' Selects a genome
    '''
    self.send("genome %s" % genome)
  
  def goto(self,locus):
    ''' Scrolls to a single locus or space-delimited list of loci. If a list is provided,
        these loci will be displayed in a split screen view. Use any syntax that is valid
        in the IGV search box.
    '''
    self.send("goto %s" % locus)
  
  def region(self, chr, start, end):
    ''' Not implemented '''
    pass

  def maxPanelHeight(self, height):
    ''' Sets the number of vertical pixels (height) of each panel to include in image
    '''
    self.send("maxPanelHeight %d" % height)
        
  def snapshot(self,filename=None):
    ''' Saves a snapshot of the IGV window to an image file. If filename is omitted,
        writes a PNG file with a filename generated based on the locus. If filename is
        specified, the filename extension determines the image file format, which must be
        .png, .jpg, or .svg
    '''
    extensions = set(['png','svg','jpg'])
    if filename is None: self.send('snapshot')
    else:
      assert filename.split('.')[-1] in extensions , 'ERROR: Filename "%s" is invalid' % filename
      self.send('snapshot %s' % filename)

  def snapshotDirectory(self,dname='.'):
    ''' Sets the directory in which to write images
    '''
    import os
    assert os.path.isdir(dname)
    self.send('snapshotDirectory %s' % os.path.abspath(dname))
  
  def viewaspairs(self, trackName=None):
    ''' Set the display mode for an alignment track to "View as pairs". trackName is
        optional.
    '''
    if trackName is None: self.send('viewaspairs')
    else: self.send('viewaspairs %s' % trackName)

  def squish(self, trackName=None):
    ''' Squish a given trackName. trackName is optional, and if it is not supplied all
        annotation tracks are squished.
    '''
    if trackName is None: self.send('squish')
    else: self.send('squish %s' % trackName)

  def collapse(self, trackName=None):
    ''' Collapse a given trackName. trackName is optional, and if it is not supplied all
        annotation tracks are collapsed.
    '''
    if trackName is None: self.send('collapse')
    else: self.send('collapse %s' % trackName)

  def expand(self, trackName=None):
    ''' Expand a given trackName. trackName is optional, and if it is not supplied all
        annotation tracks are expanded.
    '''
    if trackName is None: self.send('expand')
    else: self.send('expand %s' % trackName)

  def sort(self, option, locus):
    ''' Not implemented yet '''
    pass
    
  def preference(self, trackName=None):
    ''' Not implemented '''
    pass

  def setSleepInterval(self, ms=1000):
      ''' Sets a delay (sleep) time in milliseconds.  The sleep interval is invoked 
          between successive commands.
      '''
      self.send('setSleepInterval %d' % int(ms))

