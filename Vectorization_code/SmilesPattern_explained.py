from pyparsing import Literal, Word, White, alphas, nestedExpr, quotedString,\
    cStyleComment, alphanums, nums, StringStart, StringEnd

class SmilesPattern:
    def __init__(self):
        pass

    def addRawStr(toks):
        if 'branch' in toks:
            toks['rawStr'] = toks['branch']
        else:
            toks['rawStr'] = ''.join(toks[:])

    whitespace = White().leaveWhitespace()

    ### ATOM SECTION

    # Organic subset section
    _aliphatic_organic = ( Literal('Cl').setResultsName('symbol') \
                          | Literal('Br').setResultsName('symbol') \
                          | Word('BCNOSPFI',exact=1).setResultsName('symbol') \
                          ).setResultsName('organic')
    _aromatic_organic = ( Literal('c').setResultsName('symbol') \
                         | Word('bnosp',exact=1).setResultsName('symbol') \
                         ).setResultsName('organic'))
