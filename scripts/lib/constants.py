roleName_proteinBackbone = "CA"
roleName_nucleicBackbone = "C1'"

roleNamesForBackbone = [roleName_proteinBackbone, roleName_nucleicBackbone]

roleName_proteinSideChain = "CB"
# Purines bond to the C1' of the sugar at their N9 atoms. Pyrimidines bond to the sugar C1' atom at their N1 atoms.
roleName_nucleicPurines = "N1"
roleName_nucleicPyrimidines = "N9"

roleNamesForSideChain = [
    roleName_proteinSideChain,
    roleName_nucleicPurines,
    roleName_nucleicPyrimidines,
]

metalsTempFactor = 50.00
donorsTempFactor = 40.00
ligandsTempFactor = 30.00
neighborDonorsTempFactor = 20.00
neighborLigandsTempFactor = 10.00
