import datetime,os

def GetDatetimeFromMyFormat(l):
  # l = "2018-5-17 19:18:45"
  l_words = l.split()
  l_days = l_words[0].split('-')
  l_times = l_words[1].split(':')
  out = datetime.datetime(int(l_days[0]),int(l_days[1]),int(l_days[2]),int(l_times[0]),int(l_times[1]),int(l_times[2]))
  return out
