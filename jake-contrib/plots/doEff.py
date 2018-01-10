
def doEff(total, acc, name):
  eff = acc.Clone(name)

  eff.Divide(total)
  return eff

