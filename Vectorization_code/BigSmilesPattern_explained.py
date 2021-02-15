from SmilesPattern import SmilesPattern
from pyparsing import Word, Literal, nestedExpr, nums, StringStart, StringEnd,\
    printables, alphanums
from utility import flatten_list

class BigSmilesPattern(SmilesPattern):
    
    ### DEFINITION OF PATTERNS INVOLVED IN BigSmiles_Bond
    _BigSmilesBondChar = '$<>'
    _BondDesc = Word(_BigSmilesBondChar, exact=1).setResultsName('BigSmiles_BondType') + \
        ( Word(nums, exact=1).setResultsName('BigSmiles_Bondid') | \
           Literal('%') + Word(nums, exact=2).setResultsName('BigSmiles_Bondid') ) \
         * (0,1)
    _bigSmilesBond = (Literall('[') + _BondDesc.setResultsName('BigSmiles_Bond') +\
                      Literal(']'))

    ### DEFINITIONS OF PATTERNS INVOLVED IN Argument_Smiles
    # Redefinition for the elements used in parsing of Augmented Smiles strings
    _AugmentedSmilesChar = SmilesPattern._SmilesChar | _BigSmilesBond
