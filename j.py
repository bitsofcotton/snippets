import sys

a = ""
for line in sys.stdin:
  for l in line:
    # N.B. thanks to https://arc-tech.hatenablog.com/entry/2021/01/20/105620
    try:
      if(('ぁ' <= l and l <= 'ん') or l == 'ー' or \
         ('ァ' <= l and l <= 'ヶ') or \
         ('一' <= l and l <= '龠') ):
        a += str(l)
    except:
      pass
print(a)

